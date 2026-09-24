import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from benchmark_tools import admit_recovered_orthomcl as module


@pytest.fixture
def execution(tmp_path):
    executor, base = tmp_path / "executor", tmp_path / "base"
    source = executor / "benchmark_tools/run_recovered_orthomcl.py"
    source.parent.mkdir(parents=True)
    shutil.copyfile(Path(module.__file__).with_name("run_recovered_orthomcl.py"), source)
    (base / "inputs").mkdir(parents=True)
    for name in ("native.log", "native.time.txt", "inputs/staging.json"):
        (base / name).write_text("fixture")
    scheduler = dict(JobIDRaw="125", State="COMPLETED", ExitCode="0:0", NodeList="bizon", AllocCPUS="180", ReqMem="900G")
    runtime = {"fixture": "runtime"}
    tool, inputs = base / "native_tool", base / "inputs"
    report = dict(status="recovered_native_exited_zero_pending_admission", job_id="125", node="bizon",
        native_exit_code=0, accuracy_admitted=False, publication_ready=False, downstream_execution_authorized=False,
        source=module.record(source), cwd=str(base), environment={**module.environment(), "ORTHOMCL_PAIR_WORKERS": "64"},
        started_epoch=1, native_started_epoch=2, native_finished_epoch=3, finished_epoch=4,
        command=[str(module.TOOL / "venv_orthomcl/bin/perl"), "-I" + str(tool),
            str(executor / "benchmark_tools/run_orthomcl_perl_script.pl"), str(tool / "orthomcl.pl"),
            "--mode", "4", "--bpo_file", str(inputs / "all.bpo"), "--gg_file", str(inputs / "all.gg")])
    for key, name in (("native_log", "native.log"), ("native_timing", "native.time.txt"), ("staging", "inputs/staging.json")):
        report[key] = module.record(base / name)
    for key in ("runtime_before", "runtime_after"):
        report[key] = dict(status="dedicated_bpo_python_runtime_verified", manifest=runtime, mapped_files=[])
    return report, scheduler, executor, base, runtime


def test_execution(execution):
    module.validate_execution(*execution)


@pytest.mark.parametrize("key,value", [
    ("status", "corrected_orthomcl_native_exited_zero_pending_admission"), ("job_id", "123"),
    ("node", "other"), ("native_exit_code", False), ("native_exit_code", 1),
    ("accuracy_admitted", True), ("publication_ready", True), ("downstream_execution_authorized", True),
    ("source", {}), ("command", []), ("cwd", "/wrong"), ("environment", {}), ("staging", {}),
    ("native_log", {}), ("native_timing", {}), ("started_epoch", 0), ("finished_epoch", 2),
    ("native_finished_epoch", float("nan")), ("native_started_epoch", True)])
def test_execution_rejection(execution, key, value):
    execution[0][key] = value
    with pytest.raises(ValueError):
        module.validate_execution(*execution)


@pytest.mark.parametrize("key,value", [("State", "RUNNING"), ("ExitCode", "1:0"),
    ("NodeList", "other"), ("AllocCPUS", "2"), ("ReqMem", "64G")])
def test_scheduler_rejection(execution, key, value):
    execution[1][key] = value
    with pytest.raises(ValueError):
        module.validate_execution(*execution)


@pytest.mark.parametrize("key", ["runtime_before", "runtime_after"])
def test_runtime_rejection(execution, key):
    execution[0][key]["manifest"] = {}
    with pytest.raises(ValueError):
        module.validate_execution(*execution)


@pytest.mark.parametrize("problem", [None, "source", "record", "threads", "path", "patch"])
def test_source_reconstruction(tmp_path, monkeypatch, problem):
    tool = tmp_path / "native_tool"
    tool.mkdir()
    script, config = tmp_path / "original.pl", tmp_path / "original.pm"
    script.write_text("original script")
    config.write_text('our $PATH_TO_ORTHOMCL = "old";\nour $ORTHOMCL_DATA_DIR = "old";\nour $BLAST_NOCPU = 1;\n')
    originals = [module.record(script), module.record(config)]
    monkeypatch.setattr(module, "parallelize_source", lambda text: text + "\npatched")
    (tool / "orthomcl.pl").write_text("original script\npatched")
    (tool / "orthomcl_module.pm").write_text(
        f'our $PATH_TO_ORTHOMCL = "{tool}/";\nour $ORTHOMCL_DATA_DIR = "{tmp_path}/inputs/";\nour $BLAST_NOCPU = 180;\n')
    sources = dict(originals=originals, threads=180, pair_parallel_patch=True,
                   tool_directory=str(tool), data_directory=str(tmp_path / "inputs"),
                   configured_sources=[module.record(tool / n) for n in ("orthomcl.pl", "orthomcl_module.pm")])
    if problem == "source":
        (tool / "orthomcl.pl").write_text("changed science")
        sources["configured_sources"][0] = module.record(tool / "orthomcl.pl")
    elif problem == "record":
        sources["configured_sources"] = []
    elif problem == "threads":
        sources["threads"] = 1
    elif problem == "path":
        sources["data_directory"] = "/wrong"
    elif problem == "patch":
        sources["pair_parallel_patch"] = False
    if problem:
        with pytest.raises(ValueError):
            module.validate_sources(dict(native_sources=sources), originals, tmp_path)
    else:
        module.validate_sources(dict(native_sources=sources), originals, tmp_path)


@pytest.mark.parametrize("problem", [None, "groups", "inventory", "missing_matrix", "symlink", "cache"])
def test_outputs(tmp_path, monkeypatch, problem):
    native = tmp_path / "native_tool/run"
    (native / "tmp").mkdir(parents=True)
    for name in ("all_orthomcl.out", "tmp/all_ortho.idx", "tmp/all_ortho.mtx", "tmp/all_ortho.mcl"):
        (native / name).write_text("fixture")
    if problem == "missing_matrix":
        (native / "tmp/all_ortho.mtx").unlink()
    elif problem == "symlink":
        (native / "linked").symlink_to(native / "all_orthomcl.out")
    report = dict(native_groups=module.record(native / "all_orthomcl.out"),
                  native_outputs=[module.record(p) for p in sorted(native.rglob("*")) if p.is_file()], pair_caches=[])
    monkeypatch.setattr(module, "pair_cache_records", lambda *a: ["changed"] if problem == "cache" else [])
    if problem == "groups":
        (native / "all_orthomcl.out").write_text("changed")
    elif problem == "inventory":
        (native / "new").write_text("new")
    if problem:
        with pytest.raises(ValueError):
            module.validate_outputs(report, tmp_path)
    else:
        assert module.validate_outputs(report, tmp_path) == native


def test_isolated_cli(tmp_path):
    subprocess.run([sys.executable, "-I", "-B", module.__file__, "--help"], cwd=tmp_path,
                   check=True, capture_output=True, text=True)


@pytest.mark.parametrize("problem", [None, "search_coverage", "audit", "counts", "indexes", "changed_input", "runtime_after"])
def test_admission_flow(tmp_path, monkeypatch, problem):
    root = tmp_path
    base = root / "benchmarks/results/qfo_blast_recovery_native_v1"
    (base / "inputs").mkdir(parents=True)
    native = base / "native_tool/run"
    native.mkdir(parents=True)
    (native / "all_orthomcl.out").write_text("fixture")
    parent = root / "bpo.json"
    parent.write_text("{}")
    parent_record = module.record(parent)
    for name in ("inputs/staging.json", "native.log", "native.time.txt"):
        (base / name).write_text("{}")
    results = root / "benchmark_tools/results"
    results.mkdir(parents=True)
    source_path = results / "qfo_corrected_orthomcl_native_sources_20260918.json"
    source_path.write_text("{}")
    report = dict(admission=parent_record, admission_job="124", admission_scheduler={},
        query_coverage={"failed_queries": 2}, checked_records=[parent_record],
        staging=module.record(base / "inputs/staging.json"), native_log=module.record(base / "native.log"),
        native_timing=module.record(base / "native.time.txt"), staged_inputs={},
        native_sources={"configured_sources": []}, native_groups=module.record(native / "all_orthomcl.out"),
        content=dict(groups=1, grouped_proteins=2, ungrouped_proteins=984135))
    (base / "report.json").write_text(json.dumps(report))
    executor = root / "benchmarks/work/executor"
    monkeypatch.setattr(module, "execution_identity", lambda: dict(admission_job_id="126", node="bizon"))
    monkeypatch.setattr(module, "completed", lambda *a: ({}, ""))
    monkeypatch.setattr(module, "frozen", lambda *a: None)
    monkeypatch.setattr(module, "validate_execution", lambda *a: None)
    monkeypatch.setattr(module, "validate_staging", lambda *a: None)
    monkeypatch.setattr(module, "validate_sources", lambda *a: None)
    monkeypatch.setattr(module, "validate_outputs", lambda *a: native)
    monkeypatch.setattr(module, "verify", lambda *a: None)
    runtime_calls = []
    def runtime(*args):
        runtime_calls.append(1)
        if problem == "runtime_after" and len(runtime_calls) == 2:
            raise ValueError("runtime changed")
        return {"manifest": {}}
    monkeypatch.setattr(module, "verify_runtime", runtime)
    original_read = module.read_frozen
    def read(path, digest):
        if path == source_path:
            return dict(originals=[], checked_records=[])
        if path.parent == results:
            return {}
        return original_read(path, digest)
    monkeypatch.setattr(module, "read_frozen", read)
    evidence = dict(scheduler={}, query_coverage={"failed_queries": 2}, admission=parent_record,
                    inputs={}, index_validation={"records": 3}, checked_records=[parent_record])
    if problem == "search_coverage":
        evidence["query_coverage"]["failed_queries"] = 0
    monkeypatch.setattr(module, "verify_inputs", lambda *a: evidence)
    def audit(*args):
        if problem == "audit":
            raise ValueError("groups differ from partition")
        args[-1].write_text("{}")
        return dict(content=dict(input_proteins=984137, input_species=78,
            final_groups=2 if problem == "counts" else 1, grouped_proteins=2, ungrouped_input_proteins=984135))
    monkeypatch.setattr(module, "audit", audit)
    def indexes(argv, output, env, label):
        (output / "indexes.stdout").write_text(json.dumps({"records": 2 if problem == "indexes" else 3}))
        (output / "indexes.stderr").write_text("")
        if problem == "changed_input":
            parent.write_text("changed")
    monkeypatch.setattr(module, "native_step", indexes)
    destination = root / "admission"
    args = (root, 125, executor, "commit", executor, "bpo_commit", destination)
    if problem:
        with pytest.raises(ValueError):
            module.admit(*args)
        if problem == "search_coverage":
            assert not destination.exists()
            return
        result = json.loads((destination / "report.json").read_text())
        assert result["status"] == "recovered_native_validation_failed"
        assert result["conversion_authorized"] is False
    else:
        result = module.admit(*args)
        assert result["status"] == "recovered_orthomcl_native_outputs_admitted"
        assert result["conversion_authorized"] is True
        assert result["pair_semantics"] == "cross_species_final_group_cliques"
        with pytest.raises(FileExistsError):
            module.admit(*args)
    assert result["admission_job_id"] == "126"
    assert result["accuracy_admitted"] is False
    assert result["publication_ready"] is False
    assert result["query_coverage"]["failed_queries"] == 2
