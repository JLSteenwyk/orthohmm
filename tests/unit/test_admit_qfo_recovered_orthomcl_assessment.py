import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from benchmark_tools import admit_qfo_recovered_orthomcl_assessment as module


@pytest.mark.parametrize("problem", [None, "pending", "exit", "bool_exit", "publication", "preflight",
    "provenance", "outputs", "trace", "metrics", "mutation", "existing"])
def test_admission(tmp_path, monkeypatch, problem):
    root = tmp_path
    executor = root / "benchmarks/work/executor"
    helpers = executor / "benchmark_tools"
    helpers.mkdir(parents=True)
    source = helpers / "run_qfo_recovered_orthomcl_assessment.py"
    shutil.copyfile(Path(module.__file__).with_name(source.name), source)
    directory = root / "benchmarks/results/qfo_blast_recovery_assessment_v1"
    directory.mkdir(parents=True)
    results = root / "scoring"
    (results / "stats").mkdir(parents=True)
    trace = results / "stats/trace_1.txt"
    trace.write_text("fixture trace")
    metric = results / "metric.json"
    metric.write_text("{}")
    (directory / "scoring.log").write_text("fixture log")
    env = root / "environment.json"
    env.write_text(json.dumps(dict(pipeline=str(root / "pipeline"))))
    monkeypatch.setattr(module, "ENV_SHA", module.record(env)["sha256"])
    expected = dict(status="prepared_unrun", source=module.record(source), verified_records=[module.record(env)],
        results=str(results), environment_manifest=module.record(env), accuracy_admitted=False, publication_ready=False,
        conversion_scheduler={}, pairs_manifest={}, command=["frozen"], stage=dict(
            participant="qfo_corrected_orthomcl_recovered", query_coverage={"failed_queries": 2},
            content={"ungrouped_input_proteins": 4}, semantics="cross_species_final_group_cliques", group_audit={}))
    preflight = {**expected, "status": "running", "job_id": "128"}
    report = {**preflight, "status": "process_succeeded_pending_independent_admission", "exit_code": 0,
              "outputs": [module.record(p) for p in sorted(results.rglob("*")) if p.is_file()],
              "log": module.record(directory / "scoring.log")}
    if problem == "exit":
        report["exit_code"] = 1
    elif problem == "bool_exit":
        report["exit_code"] = False
    elif problem == "publication":
        report["publication_ready"] = preflight["publication_ready"] = True
    elif problem == "preflight":
        preflight["command"] = ["changed"]
    elif problem == "provenance":
        report["command"] = preflight["command"] = ["changed"]
    elif problem == "outputs":
        metric.write_text("changed")
    (directory / "results.json").write_text(json.dumps(report))
    (directory / "preflight.json").write_text(json.dumps(preflight))
    scheduler = dict(State="RUNNING" if problem == "pending" else "COMPLETED", ExitCode="0:0",
                     NodeList="bizon", AllocCPUS="8", JobIDRaw="128", ReqMem="64G")
    monkeypatch.setattr(module, "completed", lambda *a: (scheduler, "fixture"))
    monkeypatch.setattr(module, "frozen", lambda *a: None)
    def prepare(*args, **kwargs):
        assert kwargs == dict(require_fresh=False, helpers=helpers)
        return expected
    monkeypatch.setattr(module, "prepare", prepare)
    def validate_trace(text):
        if problem == "trace":
            raise ValueError("incomplete trace")
        return ["fixture"]
    monkeypatch.setattr(module, "validate_trace", validate_trace)
    def validate_directory(*args):
        if problem == "metrics":
            raise ValueError("invalid endpoint")
        if problem == "mutation":
            env.write_text("changed")
        return {"fixture": "scores"}, [metric]
    monkeypatch.setattr(module, "validate_directory", validate_directory)
    output = root / "admission.json"
    if problem == "existing":
        output.write_text("existing")
    args = (root, 128, 127, "digest", executor, "commit", executor, "converter", output)
    if problem:
        with pytest.raises((ValueError, FileExistsError)):
            module.admit(*args)
        if problem == "existing":
            assert output.read_text() == "existing"
        else:
            assert not output.exists()
    else:
        result = module.admit(*args)
        assert result["status"] == "recovered_orthomcl_assessment_admitted"
        assert result["accuracy_admitted"] is True
        assert result["publication_ready"] is False
        assert result["query_coverage"]["failed_queries"] == 2
        assert result["pair_semantics"] == "cross_species_final_group_cliques"


def test_cli_and_pin(tmp_path):
    assert module.record(Path(module.__file__).with_name("run_qfo_recovered_orthomcl_assessment.py"))["sha256"] == module.RUNNER_SHA
    subprocess.run([sys.executable, "-I", "-B", module.__file__, "--help"], cwd=tmp_path,
                   check=True, capture_output=True, text=True)
