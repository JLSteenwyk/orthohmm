import json
from pathlib import Path
import subprocess
import sys
import shutil
from types import SimpleNamespace

import pytest

from benchmark_tools import run_qfo_recovered_orthomcl_assessment as module


def fixture():
    stage = dict(status="recovered_orthomcl_pairs_prepared_unscored", accuracy_evaluated=False,
        publication_ready=False, method="orthomcl", participant="qfo_corrected_orthomcl_recovered",
        semantics=module.SEMANTICS, job_id="127", total_pairs=1, retained_pairs=1,
        removed_mapping_pairs=0, native_duplicate_relations=0,
        pairs=dict(bytes=4, sha256="same"), filtered_pairs=dict(bytes=4, sha256="same"),
        content=dict(total_pairs=1, final_groups=1, grouped_proteins=2, ungrouped_input_proteins=984135))
    scheduler = dict(JobIDRaw="127", State="COMPLETED", ExitCode="0:0", NodeList="bizon", AllocCPUS="2", ReqMem="64G")
    return stage, scheduler


def test_stage():
    module.validate_stage(*fixture())


@pytest.mark.parametrize("key,value", [("status", "corrected_orthomcl_pairs_prepared_unscored"),
    ("accuracy_evaluated", True), ("publication_ready", True), ("method", "other"),
    ("participant", "qfo_corrected_orthomcl"), ("semantics", "graph_edges"), ("job_id", "1"),
    ("total_pairs", 0), ("retained_pairs", 0), ("removed_mapping_pairs", 1),
    ("native_duplicate_relations", 1), ("total_pairs", True)])
def test_stage_rejection(key, value):
    stage, scheduler = fixture()
    stage[key] = value
    with pytest.raises(ValueError):
        module.validate_stage(stage, scheduler)


@pytest.mark.parametrize("key", ["State", "ExitCode", "NodeList", "AllocCPUS", "ReqMem"])
def test_scheduler_rejection(key):
    stage, scheduler = fixture()
    scheduler[key] = "wrong"
    with pytest.raises(ValueError):
        module.validate_stage(stage, scheduler)


@pytest.mark.parametrize("problem", ["sha", "bytes", "scope", "counts", "keys", "bool"])
def test_pair_content_rejection(problem):
    stage, scheduler = fixture()
    if problem in {"sha", "bytes"}:
        stage["filtered_pairs"]["sha256" if problem == "sha" else "bytes"] = "changed"
    elif problem == "scope":
        stage["content"]["ungrouped_input_proteins"] = 1
    elif problem == "counts":
        stage["content"]["total_pairs"] = 2
    elif problem == "keys":
        stage["content"]["extra"] = 1
    elif problem == "bool":
        stage["content"]["final_groups"] = True
    with pytest.raises(ValueError):
        module.validate_stage(stage, scheduler)


@pytest.mark.parametrize("problem", [None, "check_only", "allocation", "exit", "mutation"])
def test_run(tmp_path, monkeypatch, problem):
    source = tmp_path / "source.py"
    source.write_text("fixture")
    output, results = tmp_path / "output", tmp_path / "results"
    report = dict(cwd=str(output), results=str(results), command=["fixture"],
                  environment_overrides={"FROZEN": "1"}, source=module.record(source), verified_records=[],
                  status="prepared_unrun", accuracy_admitted=False, publication_ready=False)
    monkeypatch.setattr(module, "prepare", lambda *a: report)
    monkeypatch.setenv("SLURM_JOB_ID", "128")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "1" if problem == "allocation" else "8")
    monkeypatch.setenv("SLURM_MEM_PER_NODE", "65536")
    monkeypatch.setattr(module.os, "uname", lambda: SimpleNamespace(nodename="bizon"))
    calls = []
    def score(argv, **kwargs):
        calls.append(argv)
        assert kwargs["env"]["FROZEN"] == "1"
        results.mkdir()
        (results / "fixture.json").write_text("{}")
        if problem == "mutation":
            source.write_text("changed")
        return SimpleNamespace(returncode=1 if problem == "exit" else 0)
    monkeypatch.setattr(module.subprocess, "run", score)
    args = (tmp_path, "digest", 127, tmp_path / "executor", "commit")
    if problem in {"allocation", "exit", "mutation"}:
        with pytest.raises((ValueError, RuntimeError)):
            module.run(*args)
        if problem == "allocation":
            assert not output.exists()
            assert not calls
            return
        result = json.loads((output / "results.json").read_text())
        assert result["status"] == "failed"
    else:
        result = module.run(*args, check_only=problem == "check_only")
        if problem == "check_only":
            assert not calls
            assert not output.exists()
        else:
            assert result["status"] == "process_succeeded_pending_independent_admission"
            assert len(result["outputs"]) == 1
    assert result["accuracy_admitted"] is False
    assert result["publication_ready"] is False


def test_pin_and_isolated_cli(tmp_path):
    assert module.record(Path(module.__file__).with_name("prepare_recovered_orthomcl_pairs.py"))["sha256"] == module.CONVERTER_SHA
    subprocess.run([sys.executable, "-I", "-B", module.__file__, "--help"], cwd=tmp_path,
                   check=True, capture_output=True, text=True)


@pytest.mark.parametrize("problem", [None, "coverage", "groups", "mapping", "existing", "unbound", "changed_pair"])
def test_prepare(tmp_path, monkeypatch, problem):
    root = tmp_path
    stage, scheduler = fixture()
    executor = root / "benchmarks/work/executor"
    source = executor / "benchmark_tools/prepare_recovered_orthomcl_pairs.py"
    source.parent.mkdir(parents=True)
    shutil.copyfile(Path(module.__file__).with_name(source.name), source)
    directory = root / "benchmarks/results/qfo_blast_recovery_pairs_v1"
    directory.mkdir(parents=True)
    for key, name in (("pairs", "pairs.tsv"), ("filtered_pairs", "pairs.qfo.tsv")):
        (directory / name).write_text("A\tB\n")
        stage[key] = module.record(directory / name)
    content = dict(input_proteins=984137, input_species=78, final_groups=1,
                   grouped_proteins=2, ungrouped_input_proteins=984135, cross_species_clique_pairs=1)
    native = dict(status="recovered_orthomcl_native_outputs_admitted", conversion_authorized=True,
        admission_job_id="126", node="bizon", allocated_cpus=2, memory_mib=65536,
        accuracy_admitted=False, publication_ready=False, pair_semantics=module.SEMANTICS,
        scheduler=dict(State="COMPLETED", ExitCode="0:0", NodeList="bizon", AllocCPUS="180", ReqMem="900G"),
        content=content, query_coverage={"failed_queries": 2},
        source=module.record(Path(module.__file__).with_name("admit_recovered_orthomcl.py")), checked_records=[], outputs=[])
    prior = root / "benchmarks/results/qfo_blast_recovery_native_admission_v1/report.json"
    prior.parent.mkdir()
    prior.write_text(json.dumps(native))
    group = dict(status="native_final_groups_match_mcl_partition", content=content,
                 accuracy_admitted=False, publication_ready=False, checked_records=[])
    if problem == "groups":
        group["status"] = "failed"
    (directory / "groups.json").write_text(json.dumps(group))
    mapping = root / "mapping.json.gz"
    mapping.write_text("fixture")
    env = root / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    env.parent.mkdir(parents=True)
    manifest = dict(reference_files=[module.record(mapping)], environment_overrides={})
    env.write_text(json.dumps(manifest))
    monkeypatch.setattr(module, "ENV_SHA", module.record(env)["sha256"])
    stage.update(source=module.record(source), admission=module.record(prior), admission_scheduler={},
                 group_audit=module.record(directory / "groups.json"), query_coverage={"failed_queries": 2},
                 mapping=module.record(mapping), checked_records=[module.record(source), module.record(prior)])
    if problem == "coverage":
        stage["query_coverage"]["failed_queries"] = 0
    elif problem == "mapping":
        stage["mapping"] = {**stage["mapping"], "sha256": "changed"}
    elif problem == "unbound":
        stage["checked_records"].pop()
    (directory / "results.json").write_text(json.dumps(stage))
    digest = module.record(directory / "results.json")["sha256"]
    if problem == "changed_pair":
        (directory / "pairs.tsv").write_text("mutated")
    output = root / "benchmarks/results/qfo_blast_recovery_assessment_v1"
    if problem == "existing":
        output.mkdir()
    monkeypatch.setattr(module, "completed", lambda job, *a: (scheduler if job == 127 else {}, ""))
    monkeypatch.setattr(module, "frozen", lambda *a: None)
    monkeypatch.setattr(module, "environment_records", lambda m: m["reference_files"])
    monkeypatch.setattr(module, "command_for", lambda r, s, m, w, o: [s["filtered_pairs"]["path"], s["participant"], str(w), str(o)])
    if problem:
        with pytest.raises((ValueError, FileExistsError)):
            module.prepare(root, digest, 127, executor, "commit")
    else:
        result = module.prepare(root, digest, 127, executor, "commit")
        assert result["command"][1] == "qfo_corrected_orthomcl_recovered"
        assert result["work"] == str(root / "qfo_benchmark/w/qc_mcr")
        assert result["accuracy_admitted"] is False
    if problem != "existing":
        assert not output.exists()
