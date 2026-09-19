import copy
import json
import sys

import pytest

from benchmark_tools import admit_qfo_cpm_variant as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_admit_qfo_corrected_replay import parent_fixture


@pytest.mark.parametrize("problem", [None, "missing", "duplicate", "running", "failed", "exit", "node", "cpus", "index"])
def test_completed_task(problem):
    row = ["21960_0", "21961", "COMPLETED", "0:0", "00:10:00", "bizon", "32"]
    if problem in ("running", "failed"):
        row[2] = problem.upper()
    elif problem == "exit":
        row[3] = "1:0"
    elif problem == "node":
        row[5] = "other"
    elif problem == "cpus":
        row[6] = "2"
    accounting = "JobID|JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n"
    accounting += ("|".join(row) + "\n") * (0 if problem == "missing" else 2 if problem == "duplicate" else 1)
    if problem:
        with pytest.raises(ValueError):
            module.completed_task(accounting, True if problem == "index" else 0)
    else:
        assert module.completed_task(accounting, 0)["JobIDRaw"] == "21961"


@pytest.mark.parametrize("index", [0, 1])
@pytest.mark.parametrize("problem", [None, "identity", "status", "exit", "context", "worker_source", "helpers",
    "time", "command", "cwd", "resolution", "profiles", "species", "output", "scoring", "accuracy"])
def test_parent(tmp_path, index, problem):
    parent, worker, replay, context, _, scheduler, root, executor, _ = parent_fixture(tmp_path)
    executor.joinpath("benchmark_tools").mkdir(parents=True)
    for name in ("run_qfo_cpm_variant.py", "run_qfo_cpm_control.py"):
        (executor / "benchmark_tools" / name).write_text(name)
    context.update(resolution=[.08, .12][index], output_root=str(tmp_path / "output"))
    source = record(executor / "benchmark_tools/run_qfo_cpm_variant.py")
    shared = record(executor / "benchmark_tools/run_qfo_cpm_control.py")
    helpers = [source, shared]
    parent.update(status="cpm_variant_replay_complete_pending_admission", arm=module.ARMS[index], index=index,
        executor_commit=module.EXECUTOR, context=copy.deepcopy(context), source=source, helpers=helpers,
        publication_ready=False, worker_command=[sys.executable, "-B", source["path"], "--root", str(root),
                                               "--index", str(index), "--worker"])
    worker.update(context=copy.deepcopy(context), source=shared)
    replay["parameters"]["cpm_resolution"] = context["resolution"]
    for stage, filename in zip(replay["stages"], module.FILENAMES):
        stage["output"] = {"path": str(tmp_path / "output/replay" / filename)}
    mutations = {
        "identity": (parent, "arm", "control"), "status": (parent, "status", "failed"),
        "exit": (parent, "exit_code", False), "context": (worker, "context", {}),
        "worker_source": (worker, "source", source), "helpers": (parent, "helpers", []),
        "time": (parent, "finished_epoch", float("nan")), "command": (parent, "worker_command", []),
        "cwd": (replay, "cwd", "other"), "resolution": (replay["parameters"], "cpm_resolution", .1),
        "profiles": (replay["counts"], "profiles_built", 0), "species": (replay["counts"], "species", 77),
        "output": (replay["stages"][0]["output"], "path", "other"),
        "scoring": (replay["stages"][0], "official_orthobench", {}), "accuracy": (parent, "accuracy_evaluated", True)}
    args = (parent, worker, replay, context, scheduler, root, executor, index, helpers)
    if problem:
        target, key, value = mutations[problem]
        target[key] = value
        with pytest.raises(ValueError):
            module.validate_parent(*args)
    else:
        module.validate_parent(*args)


@pytest.mark.parametrize("problem", [None, "different", "missing", "duplicate", "count", "stage_order", "hash"])
def test_partition_comparison(tmp_path, problem):
    stages, control = [], []
    for index, label in enumerate(module.STAGES):
        old, new = tmp_path / f"old{index}", tmp_path / f"new{index}"
        old.write_text("a b\nc\n")
        new.write_text({"different": "a\nb c\n", "missing": "a b\n", "duplicate": "a b\nb c\n"}.get(problem, "c\nb a\n"))
        stages.append({"label": label, "clusters": 3 if problem == "count" else 2, "output": record(new)})
        control.append({"label": label, "output": record(old)})
    if problem == "stage_order":
        control.reverse()
    if problem == "hash":
        new.write_text("changed")
    if problem in (None, "different"):
        coverage, comparisons = module.partitions({"stages": stages}, control, set("abc"))
        assert len(coverage) == len(comparisons) == 4
        assert all(r["partition_equal"] is (problem is None) for r in comparisons)
        assert all(r["variant_only_groups"] == (0 if problem is None else 2) for r in comparisons)
    else:
        with pytest.raises(ValueError):
            module.partitions({"stages": stages}, control, set("abc"))


@pytest.mark.parametrize("state", ["PENDING", "RUNNING", "FAILED"])
def test_live_scheduler_precedes_files(tmp_path, monkeypatch, state):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
        f"JobID|JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n21960_0|21961|{state}|0:0|00:01:00|bizon|32\n")
    output = tmp_path / "admission.json"
    with pytest.raises(ValueError, match="completed"):
        module.admit(tmp_path, 0, output)
    assert not output.exists()


@pytest.mark.parametrize("problem", [None, "fresh", "runtime", "coverage", "checkpoint", "command"])
def test_admission_orchestration(tmp_path, monkeypatch, problem):
    def write(path, value):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value))
        return record(path)
    root, directory, launcher = tmp_path, tmp_path / "output", tmp_path / "launcher"
    executor = root / "benchmarks/work/publication_qfo_cpm_variant_v1"
    write(executor / "benchmark_tools/run_qfo_cpm_variant.py", {})
    scientific = write(launcher / "benchmark_tools/replay_high_sensitivity.py", {})
    auditor = write(launcher / "benchmark_tools/audit_accuracy_checkpoint.py", {})
    names = root / "names.txt"
    names.write_text("".join(f"g{i}\n" for i in range(984137)))
    numeric = {"status": "verified", "summary": {"genes": 984137}, "manifest": {"fixture": True}}
    admission_record = write(root / "native.json", {"content": {"numeric_checkpoint": numeric}})
    plan = {"runtime": {}}
    plan_record = write(root / "plan.json", plan)
    context = {"output_root": str(directory), "cwd": str(launcher), "checked_records": []}
    monkeypatch.setattr(module, "corrected_evidence", lambda *a: (plan, plan_record, admission_record, record(names)))
    monkeypatch.setattr(module, "evidence", lambda *a: context)
    control_report = {"stage_comparisons": []}
    authorization = {"executor": str(root / "admitter"), "report": control_report,
        "report_record": write(root / "control.json", control_report), "checked_records": []}
    monkeypatch.setattr(module, "control_evidence", lambda *a: authorization)
    monkeypatch.setattr(module, "validate_parent", lambda *a: None)
    def audit(*args, **kwargs):
        assert kwargs == {"cpm_arm": "cpm_low"}
        return {"clustering": ["fixture"], "checked_records": []}
    monkeypatch.setattr(module, "audit", audit)
    monkeypatch.setattr(module, "verify", lambda *a: {"changed": True} if problem == "runtime" else {})
    partition = write(directory / "partition.txt", {})
    coverage = [{"label": "fixture", "groups": 1, "output": partition}]
    def partitions(replay, control, universe):
        assert len(universe) == 984137
        return coverage, [{"partition_equal": False}]
    monkeypatch.setattr(module, "partitions", partitions)
    worker_record = write(directory / "checked_worker.json", {})
    replay_input = {**numeric, "accuracy_evaluated": False, "auditor": auditor}
    if problem == "checkpoint":
        replay_input["summary"] = {}
    replay_record = write(directory / "replay.json", {"source": scientific, "input": replay_input})
    fresh_record = write(directory / "fresh_control_admission.json", {} if problem == "fresh" else control_report)
    expected = [sys.executable, "-B", str(root / "admitter/benchmark_tools/admit_qfo_cpm_control.py"),
                "--root", str(root), "--output", fresh_record["path"]]
    parent = {"worker": worker_record, "replay": replay_record, "fresh_control_admission": fresh_record,
        "authorization": authorization, "admission_command": [] if problem == "command" else expected,
        "runtime_before": {}, "runtime_after": {}, "coverage": [] if problem == "coverage" else coverage}
    write(directory / "results.json", parent)
    for name in ("control_validation.log", "replay.log", "time.txt"):
        write(directory / name, {})
    def check_output(command, **kwargs):
        if command[0] == "sacct":
            return "JobID|JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n21960_0|21961|COMPLETED|0:0|00:01:00|bizon|32\n"
        return module.EXECUTOR
    monkeypatch.setattr(module.subprocess, "check_output", check_output)
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k: None)
    output = root / "admitted.json"
    if problem:
        with pytest.raises(ValueError):
            module.admit(root, 0, output)
        assert not output.exists()
    else:
        result = module.admit(root, 0, output)
        assert json.loads(output.read_text()) == result
        assert result["status"] == "cpm_variant_replay_admitted_unscored"
        assert result["control_comparison"] == [{"partition_equal": False}]
        assert result["accuracy_evaluated"] is False and result["publication_ready"] is False
        with pytest.raises(FileExistsError):
            module.admit(root, 0, output)
