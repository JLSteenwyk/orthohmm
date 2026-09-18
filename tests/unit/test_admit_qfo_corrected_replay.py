import copy
import json
import sys

import pytest

from benchmark_tools import admit_qfo_corrected_replay as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def parent_fixture(root):
    executor = root / "executor"
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    helpers = [{"path": str(executor / "benchmark_tools" / name), "sha256": name, "bytes": 1} for name in module.HELPERS]
    plan_record = {"path": str(root / "plan.json"), "sha256": "fixture", "bytes": 1}
    plan = {"native_command": ["fixture"], "expected_stages": module.STAGES, "cwd": str(launcher)}
    worker = {"status": "corrected_checked_replay_worker_returned", "accuracy_evaluated": False,
              "source": helpers[0], "plan": plan_record, "replay_source": {"fixture": True},
              "calls": [{"stage": s, "status": "checked", "exit_code": 0}
                        for s in ("initial", "multipass", "profile_base", "profile_expanded")]}
    replay = {"command": plan["native_command"], "source": worker["replay_source"], "cwd": str(launcher),
              "parameters": {"accuracy_profile": "high_sensitivity", "cpm_resolution": .1, "profile_expansion": True,
                  "profile_iterations": 1, "jackknife_profile_thresholds": False,
                  "jackknife_single_copy_profiles": False, "profile_min_species": 1,
                  "matrix": "BLOSUM62", "leiden_seed": 4},
              "counts": {"genes": 984137, "species": 78, "profiles_built": 1},
              "stages": [{"label": s} for s in module.STAGES]}
    parent = {"status": "corrected_checked_replay_complete_pending_admission", "job_id": "17",
              "exit_code": 0, "accuracy_evaluated": False, "executor_commit": module.EXECUTOR,
              "started_epoch": 1., "finished_epoch": 2., "plan": plan_record,
              "source": helpers[0], "helpers": helpers, "worker_command": [sys.executable,
                  str(executor / "benchmark_tools/run_qfo_corrected_replay.py"), "--root", str(root),
                  "--plan", plan_record["path"], "--plan-sha256", plan_record["sha256"], "--replay-worker"]}
    scheduler = {"State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon", "AllocCPUS": "32", "JobIDRaw": "17"}
    return parent, worker, replay, plan, plan_record, scheduler, root, executor, helpers


@pytest.mark.parametrize("problem", [None, "job", "node", "cpus", "scheduler_failed", "parent_failed",
    "exit", "commit", "time_nan", "time_order", "plan", "helpers", "worker_command", "command", "cwd",
    "settings", "stages", "species", "profiles", "scoring"])
def test_parent_gate(tmp_path, problem):
    args = copy.deepcopy(parent_fixture(tmp_path))
    parent, worker, replay, plan, _, scheduler, *_ = args
    if problem == "job":
        parent["job_id"] = "18"
    elif problem == "node":
        scheduler["NodeList"] = "other"
    elif problem == "cpus":
        scheduler["AllocCPUS"] = "16"
    elif problem == "scheduler_failed":
        scheduler["State"] = "FAILED"
    elif problem == "parent_failed":
        parent["status"] = "failed"
    elif problem == "exit":
        parent["exit_code"] = False
    elif problem == "commit":
        parent["executor_commit"] = "other"
    elif problem == "time_nan":
        parent["finished_epoch"] = float("nan")
    elif problem == "time_order":
        parent["finished_epoch"] = 0
    elif problem == "plan":
        worker["plan"] = {}
    elif problem == "helpers":
        parent["helpers"] = []
    elif problem == "worker_command":
        parent["worker_command"] = []
    elif problem == "command":
        replay["command"] = []
    elif problem == "cwd":
        replay["cwd"] = "other"
    elif problem == "settings":
        replay["parameters"]["cpm_resolution"] = .2
    elif problem == "stages":
        replay["stages"] = []
    elif problem == "species":
        replay["counts"]["species"] = 77
    elif problem == "profiles":
        replay["counts"]["profiles_built"] = 0
    elif problem == "scoring":
        replay["stages"][0]["official_orthobench"] = {}
    if problem:
        with pytest.raises(ValueError):
            module.validate_parent(*args)
    else:
        module.validate_parent(*args)


@pytest.mark.parametrize("problem", [None, "different", "missing", "duplicate", "foreign", "count",
                                    "path", "changed_hash", "native_missing", "native_duplicate", "stage_order"])
def test_partition_gate_and_nonequivalence(tmp_path, problem):
    output = tmp_path / "replay"
    output.mkdir()
    stages = []
    for label, filename in zip(module.STAGES, module.FILENAMES):
        path = output / filename
        path.write_text("a b\nc\n")
        stages.append({"label": label, "clusters": 2, "output": record(path)})
    native = tmp_path / "native.txt"
    native.write_text("OG0: a b\nOG1: c\n")
    final = output / module.FILENAMES[-1]
    if problem in ("different", "missing", "duplicate", "foreign"):
        final.write_text({"different": "a\nb c\n", "missing": "a b\n", "duplicate": "a b\nb c\n",
                          "foreign": "a b\nc d\n"}[problem])
        stages[-1]["output"] = record(final)
    elif problem == "count":
        stages[-1]["clusters"] = 3
    elif problem == "path":
        stages[-1]["output"] = stages[0]["output"]
    elif problem == "changed_hash":
        final.write_text("a\nb c\n")
    elif problem == "native_missing":
        native.write_text("OG0: a b\n")
    elif problem == "native_duplicate":
        native.write_text("OG0: a b\nOG1: b c\n")
    elif problem == "stage_order":
        stages.reverse()
    if problem in (None, "different"):
        coverage, comparison = module.validate_partitions(tmp_path, {"stages": stages}, set("abc"), record(native))
        assert len(coverage) == 4
        assert comparison["partition_equal"] is (problem is None)
        assert comparison["native_only_groups"] == (2 if problem else 0)
        assert comparison["replay_only_groups"] == (2 if problem else 0)
    else:
        with pytest.raises(ValueError):
            module.validate_partitions(tmp_path, {"stages": stages}, set("abc"), record(native))


@pytest.mark.parametrize("state", ["RUNNING", "PENDING", "FAILED", "CANCELLED"])
def test_scheduler_gate_precedes_missing_files(tmp_path, monkeypatch, state):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
        f"JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n17|{state}|0:0|00:01:00|bizon|32\n")
    output = tmp_path / "admission.json"
    with pytest.raises(ValueError, match="COMPLETED"):
        module.admit(tmp_path, tmp_path / "absent.json", "missing", "missing", "17", output)
    assert not output.exists()


@pytest.mark.parametrize("problem", [None, "report_hash", "runtime", "worker_hash", "comparison"])
def test_admission_orchestration(tmp_path, monkeypatch, problem):
    # Exercise orchestration with real artifact hashing. Scheduler/runtime,
    # native stage audit and partition parsing are independently tested gates.
    def write(path, value):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value))
        return record(path)
    root = tmp_path
    executor = root / "benchmarks/work/publication_qfo_corrected_replay_v1"
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    directory = root / "benchmarks/results/qfo_corrected_checked_replay_v1"
    directory.mkdir(parents=True)
    for name in module.HELPERS:
        write(executor / "benchmark_tools" / name, {})
    scientific = write(launcher / "benchmark_tools/replay_high_sensitivity.py", {})
    auditor = write(launcher / "benchmark_tools/audit_accuracy_checkpoint.py", {})
    primary = {"input_directory": str(root / "inputs")}
    primary_path = root / "benchmark_tools/results/qfo_corrected_primary_commands_20260918.json"
    primary_record = write(primary_path, primary)
    monkeypatch.setattr(module, "PLAN_SHA", primary_record["sha256"])
    checkpoint = root / "checkpoint"
    checkpoint.mkdir()
    names = checkpoint / "gene_names.txt"
    names.write_text("".join(f"g{i}\n" for i in range(984137)))
    native_record = write(root / "native.txt", {})
    numeric = {"status": "numeric_checkpoint_verified", "summary": {"genes": 984137}, "manifest": {"fixture": True}}
    admission_record = write(root / "native_admission.json", {"content": {"numeric_checkpoint": numeric,
                                                                          "native_groups": native_record}})
    plan = {"output_root": str(directory), "input_fastas": [], "primary_plan": primary_record,
            "source": record(executor / "benchmark_tools/prepare_qfo_corrected_replay.py"),
            "native_command": module.command_for(sys.executable, launcher, directory, root / "inputs", checkpoint, "checkpoint"),
            "environment_overrides": {"PYTHONPATH": str(launcher), "PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1",
                                      "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"},
            "runtime": {"fixture": True}, "checked_records": [record(names)]}
    plan_path = root / "plan.json"
    plan_record = write(plan_path, plan)
    monkeypatch.setattr(module, "corrected_evidence", lambda *a: (plan, plan_record, admission_record, record(names)))
    monkeypatch.setattr(module, "validate_admission", lambda *a: (checkpoint, "checkpoint"))
    monkeypatch.setattr(module, "validate_parent", lambda *a: None)
    monkeypatch.setattr(module, "audit_stages", lambda *a: {"clustering": ["fixture"], "checked_records": []})
    coverage, comparison = [], {"partition_equal": False, "native_only_groups": 2}
    monkeypatch.setattr(module, "validate_partitions", lambda *a: (coverage, comparison))
    monkeypatch.setattr(module, "verify", lambda *a: {"fixture": problem != "runtime"})
    worker_record = write(directory / "checked_worker.json", {})
    replay_record = write(directory / "replay.json", {"source": scientific, "stages": [],
        "input": {**numeric, "accuracy_evaluated": False, "auditor": auditor}})
    parent = {"worker": worker_record, "replay": replay_record, "runtime_before": {"fixture": True},
              "runtime_after": {"fixture": True}, "coverage": coverage, "native_partition_comparison": comparison}
    if problem == "worker_hash":
        parent["worker"] = {**worker_record, "sha256": "changed"}
    elif problem == "comparison":
        parent["native_partition_comparison"] = {"partition_equal": True}
    parent_record = write(directory / "results.json", parent)
    write(directory / "replay.log", {})
    write(directory / "time.txt", {})
    def check_output(command, **kwargs):
        if command[0] == "sacct":
            return "JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n17|COMPLETED|0:0|00:01:00|bizon|32\n"
        return module.EXECUTOR + "\n"
    monkeypatch.setattr(module.subprocess, "check_output", check_output)
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k: None)
    output = root / "admission.json"
    digest = "0" * 64 if problem == "report_hash" else parent_record["sha256"]
    if problem:
        with pytest.raises(ValueError):
            module.admit(root, plan_path, plan_record["sha256"], digest, "17", output)
        assert not output.exists()
    else:
        result = module.admit(root, plan_path, plan_record["sha256"], digest, "17", output)
        assert json.loads(output.read_text()) == result
        assert result["status"] == "corrected_checked_replay_admitted"
        assert result["native_partition_comparison"]["partition_equal"] is False
        assert result["accuracy_evaluated"] is False
        assert result["publication_ready"] is False
