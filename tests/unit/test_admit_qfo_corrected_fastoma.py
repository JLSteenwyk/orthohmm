from copy import deepcopy

import pytest

from benchmark_tools.admit_qfo_corrected_fastoma import validate_execution, tree_clades, bind_task_outputs, admit
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def execution():
    expected = {"output_root": "/out", "native_argv": ["nextflow", "run"], "cwd": "/out/run", "environment": {"NXF_OFFLINE": "true"}}
    scheduler = {"State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon", "AllocCPUS": "180", "ReqMem": "720G", "JobIDRaw": "123"}
    source = {"path": "/runner", "sha256": "source", "bytes": 10}
    result = {**deepcopy(expected), "status": "process_succeeded_pending_native_admission", "exit_code": 0,
              "job_id": "123", "node": "bizon", "source": source, "started_epoch": 1.0, "finished_epoch": 2.0,
              "accuracy_admitted": False, "native_outputs_validated": False, "publication_ready": False}
    for key, path in (("log", "/out/native.log"), ("timing", "/out/time.txt"),
                      ("trace", "/out/run/trace.txt"), ("nextflow_log", "/out/run/.nextflow.log")):
        result[key] = {"path": path, "bytes": 1, "sha256": "artifact"}
    return result, expected, scheduler, source


def test_matching_execution():
    validate_execution(*execution())


@pytest.mark.parametrize("key,value", [("status", "failed"), ("exit_code", 1), ("job_id", "other"),
    ("node", "other"), ("started_epoch", 0), ("finished_epoch", float("nan")),
    ("finished_epoch", float("inf")), ("finished_epoch", .5), ("accuracy_admitted", True),
    ("publication_ready", True), ("native_outputs_validated", True), ("source", {}),
    ("native_argv", ["changed"]), ("cwd", "/wrong"), ("environment", {})])
def test_execution_drift_rejected(key, value):
    result, expected, scheduler, source = execution()
    result[key] = value
    with pytest.raises(ValueError):
        validate_execution(result, expected, scheduler, source)


@pytest.mark.parametrize("key,value", [("State", "RUNNING"), ("ExitCode", "1:0"),
                                      ("AllocCPUS", "32"), ("ReqMem", "700G"), ("NodeList", "other")])
def test_scheduler_drift_rejected(key, value):
    result, expected, scheduler, source = execution()
    scheduler[key] = value
    with pytest.raises(ValueError):
        validate_execution(result, expected, scheduler, source)


def test_empty_or_wrong_trace_rejected():
    for field, value in (("path", "/old/trace.txt"), ("bytes", 0)):
        result, expected, scheduler, source = execution()
        result["trace"][field] = value
        with pytest.raises(ValueError):
            validate_execution(result, expected, scheduler, source)


def test_species_tree_renaming_not_topology_change(tmp_path):
    first, second, different = [tmp_path / name for name in ("a", "b", "c")]
    first.write_text("((s1:1,s2:1)N1:1,s3:1)N0;")
    second.write_text("(s3:1,(s2:1,s1:1)renamed:1)root;")
    different.write_text("((s1:1,s3:1):1,s2:1);")
    species = ["s1", "s2", "s3"]
    assert tree_clades(first, species) == tree_clades(second, species)
    assert tree_clades(first, species) != tree_clades(different, species)
    with pytest.raises(ValueError):
        tree_clades(first, ["s1", "s2"])


def bindings(tmp_path):
    output = tmp_path / "out"
    output.mkdir()
    tasks = []
    names = {"collect_subhogs": ["FastOMA_HOGs.orthoxml", "RootHOGs.tsv", "OrthologousGroups.tsv"],
             "extract_pairwise_ortholog_relations": ["orthologs.tsv.gz", "FastOMA_HOGs.orthoxml"],
             "check_input": ["species_tree_checked.nwk"], "fastoma_report": ["report.ipynb", "report.html"]}
    for process, files in names.items():
        directory = tmp_path / process
        directory.mkdir()
        tasks.append({"trace": {"name": process}, "directory": str(directory)})
        for name in files:
            (directory / name).write_text(name)
            (output / name).write_text(name)
    return {"tasks": tasks}, output


def test_published_output_bindings(tmp_path):
    tasks, output = bindings(tmp_path)
    assert len(bind_task_outputs(tasks, output, [])) == 8
    (output / "orthologs.tsv.gz").write_text("different")
    with pytest.raises(ValueError, match="differs"):
        bind_task_outputs(tasks, output, [])


@pytest.mark.parametrize("changed", [None, "query", "published", "consumed"])
def test_omamer_input_and_mapping_bindings(tmp_path, changed):
    tasks, output = bindings(tmp_path)
    omamer, infer = tmp_path / "omamer", tmp_path / "infer"
    omamer.mkdir()
    (infer / "hogmaps").mkdir(parents=True)
    (output / "hogmap").mkdir()
    tasks["tasks"].extend([{"trace": {"name": "omamer_run (sp.fa)"}, "directory": str(omamer)},
                           {"trace": {"name": "infer_roothogs"}, "directory": str(infer)}])
    (omamer / ".command.sh").write_text("omamer search --query sp.fa")
    original = tmp_path / "sp.fa"
    original.write_text(">protein\nACDE\n")
    query = omamer / "sp.fa"
    query.write_bytes(original.read_bytes())
    paths = {"query": query, "published": output / "hogmap/sp.fa.hogmap",
             "consumed": infer / "hogmaps/sp.fa.hogmap"}
    for path in [omamer / "sp.fa.hogmap", paths["published"], paths["consumed"]]:
        path.write_text("mapping")
    inputs = [record(original)]
    if changed:
        paths[changed].write_text("changed")
        with pytest.raises(ValueError):
            bind_task_outputs(tasks, output, inputs)
    else:
        assert len(bind_task_outputs(tasks, output, inputs)) == 11


def test_existing_admission_never_overwritten(tmp_path):
    path = tmp_path / "report"
    path.write_text("preserve")
    with pytest.raises(FileExistsError):
        admit(tmp_path, 1, 2, path)
    assert path.read_text() == "preserve"


def test_reviewed_trace_gate_precedes_task_access(tmp_path, monkeypatch):
    from benchmark_tools import admit_qfo_corrected_fastoma as module
    path = tmp_path / "trace.txt"
    path.write_text("unreviewed trace\n")
    monkeypatch.setattr(module, "audit_tasks", lambda *a, **k: pytest.fail("Audited unreviewed trace"))
    with pytest.raises(ValueError, match="explicitly reviewed"):
        module.reviewed_tasks(path, tmp_path / "work", ["a.fa"])
