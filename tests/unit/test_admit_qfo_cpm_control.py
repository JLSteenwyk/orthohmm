import copy
import json
import sys

import pytest

from benchmark_tools import admit_qfo_cpm_control as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_admit_qfo_corrected_replay import parent_fixture


def fixture(root):
    parent, worker, replay, context, _, scheduler, _, executor, _ = parent_fixture(root)
    source = {"path": str(executor / "benchmark_tools/run_qfo_cpm_control.py"), "sha256": "fixture", "bytes": 1}
    helpers = [{"path": str(executor / "benchmark_tools" / name)} for name in module.HELPERS]
    context["output_root"] = str(root / "output")
    parent.update(status="cpm_control_reproduced_pending_independent_admission", job_id=module.JOB,
        executor_commit=module.EXECUTOR, context=copy.deepcopy(context), source=source, helpers=helpers,
        checked_inputs=[], changed_arms_authorized=False,
        worker_command=[sys.executable, source["path"], "--root", str(root), "--worker"])
    worker.update(context=copy.deepcopy(context), source=source)
    scheduler["JobIDRaw"] = module.JOB
    for stage, filename in zip(replay["stages"], module.FILENAMES):
        stage["output"] = {"path": str(root / "output/replay" / filename)}
    return parent, worker, replay, context, scheduler, root, executor, source, helpers, []


@pytest.mark.parametrize("problem", [None, "job", "cpus", "node", "scheduler", "parent", "exit",
    "commit", "time", "source", "context", "inputs", "helpers", "authorization", "accuracy",
    "command", "worker_command", "cwd", "parameters", "species", "profiles", "path", "scoring"])
def test_parent_gate(tmp_path, problem):
    args = copy.deepcopy(fixture(tmp_path))
    parent, worker, replay, context, scheduler, *_ = args
    mutations = {
        "job": (scheduler, "JobIDRaw", "other"), "cpus": (scheduler, "AllocCPUS", "16"),
        "node": (scheduler, "NodeList", "other"), "scheduler": (scheduler, "State", "RUNNING"),
        "parent": (parent, "status", "failed"), "exit": (parent, "exit_code", False),
        "commit": (parent, "executor_commit", "other"), "time": (parent, "finished_epoch", float("nan")),
        "source": (worker, "source", {}), "context": (worker, "context", {}),
        "inputs": (parent, "checked_inputs", [{}]), "helpers": (parent, "helpers", []),
        "authorization": (parent, "changed_arms_authorized", True), "accuracy": (worker, "accuracy_evaluated", True),
        "command": (replay, "command", []), "worker_command": (parent, "worker_command", []),
        "cwd": (replay, "cwd", "other"), "parameters": (replay["parameters"], "cpm_resolution", .12),
        "species": (replay["counts"], "species", 77), "profiles": (replay["counts"], "profiles_built", 0),
        "path": (replay["stages"][0]["output"], "path", "other"),
        "scoring": (replay["stages"][0], "official_orthobench", {})}
    if problem:
        target, key, value = mutations[problem]
        target[key] = value
        with pytest.raises(ValueError):
            module.validate_parent(*args)
    else:
        module.validate_parent(*args)


@pytest.mark.parametrize("state", ["RUNNING", "PENDING", "FAILED", "CANCELLED"])
def test_scheduler_gate_precedes_output_reads(tmp_path, monkeypatch, state):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
        f"JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n{module.JOB}|{state}|0:0|00:01:00|bizon|32\n")
    output = tmp_path / "admission.json"
    with pytest.raises(ValueError, match="COMPLETED"):
        module.admit(tmp_path, output)
    assert not output.exists()


@pytest.mark.parametrize("problem", [None, "worker_hash", "runtime", "comparison", "unequal", "count", "checkpoint"])
def test_admission_orchestration(tmp_path, monkeypatch, problem):
    def write(path, value):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value))
        return record(path)
    root = tmp_path
    executor = root / "benchmarks/work/publication_qfo_cpm_control_v1"
    launcher = root / "launcher"
    directory = root / "output"
    for name in (*module.HELPERS, "run_qfo_cpm_control.py"):
        write(executor / "benchmark_tools" / name, {})
    scientific = write(launcher / "benchmark_tools/replay_high_sensitivity.py", {})
    auditor = write(launcher / "benchmark_tools/audit_accuracy_checkpoint.py", {})
    names = root / "names.txt"
    names.write_text("".join(f"g{i}\n" for i in range(984137)))
    numeric = {"status": "verified", "summary": {"genes": 984137}, "manifest": {"fixture": True}}
    admission_record = write(root / "native.json", {"content": {"numeric_checkpoint": numeric}})
    plan = {"runtime": {"fixture": True}}
    plan_record = write(root / "plan.json", plan)
    baseline = {"status": "corrected_checked_replay_admitted", "plan": plan_record, "checked_records": []}
    baseline_record = write(root / "benchmark_tools/results/qfo_corrected_replay_admission_21757.json", baseline)
    monkeypatch.setattr(module, "BASELINE_SHA", baseline_record["sha256"])
    context = {"output_root": str(directory), "cwd": str(launcher), "checked_records": []}
    monkeypatch.setattr(module, "corrected_evidence", lambda *a: (plan, plan_record, admission_record, record(names)))
    monkeypatch.setattr(module, "evidence", lambda *a: context)
    monkeypatch.setattr(module, "validate_parent", lambda *a: None)
    def audit(*args, **kwargs):
        assert kwargs == {"cpm_arm": "control"}
        return {"clustering": ["fixture"], "checked_records": []}
    monkeypatch.setattr(module, "audit", audit)
    monkeypatch.setattr(module, "verify", lambda *a: {"fixture": problem != "runtime"})
    partition = write(directory / "partition.txt", {})
    comparisons = [{"partition_equal": problem != "unequal", "groups": 1, "output": partition}]
    monkeypatch.setattr(module, "compare_partitions", lambda *a: comparisons)
    worker_record = write(directory / "checked_worker.json", {})
    replay_input = {**numeric, "accuracy_evaluated": False, "auditor": auditor}
    if problem == "checkpoint":
        replay_input["summary"] = {}
    replay_record = write(directory / "replay.json", {"source": scientific, "input": replay_input,
        "stages": [{"output": partition, "clusters": 2 if problem == "count" else 1}]})
    parent = {"worker": worker_record, "replay": replay_record, "runtime_before": plan["runtime"],
        "runtime_after": plan["runtime"], "stage_comparisons": [] if problem == "comparison" else comparisons}
    if problem == "worker_hash":
        parent["worker"] = {**worker_record, "sha256": "changed"}
    write(directory / "results.json", parent)
    write(directory / "time.txt", {})
    write(directory / "replay.log", {})
    def check_output(command, **kwargs):
        if command[0] == "sacct":
            return f"JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n{module.JOB}|COMPLETED|0:0|00:01:00|bizon|32\n"
        return module.EXECUTOR + "\n"
    monkeypatch.setattr(module.subprocess, "check_output", check_output)
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k: None)
    output = root / "admitted.json"
    if problem:
        with pytest.raises(ValueError):
            module.admit(root, output)
        assert not output.exists()
    else:
        result = module.admit(root, output)
        assert json.loads(output.read_text()) == result
        assert result["status"] == "cpm_control_reproduced_and_admitted"
        assert result["changed_arms_authorized"] == ["cpm_low", "cpm_high"]
        assert result["accuracy_evaluated"] is False and result["publication_ready"] is False
        with pytest.raises(FileExistsError):
            module.admit(root, output)
