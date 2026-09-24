import copy
import sys

import pytest

from benchmark_tools import validate_cpm_refinement_reconstruction as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record


@pytest.mark.parametrize("problem", [None, "running", "failed", "exit", "cpu", "node", "missing", "duplicate"])
def test_scheduler(problem):
    fields = ["22153", "COMPLETED", "0:0", "bizon", "1"]
    for key, index, value in (("running", 1, "RUNNING"), ("failed", 1, "FAILED"),
                             ("exit", 2, "1:0"), ("cpu", 4, "2"), ("node", 3, "other")):
        if problem == key:
            fields[index] = value
    rows = ["|".join(fields)]
    if problem == "missing":
        rows = []
    elif problem == "duplicate":
        rows *= 2
    accounting = "JobID|State|ExitCode|NodeList|AllocCPUS\n" + "\n".join(rows)
    if problem:
        with pytest.raises(ValueError):
            module.completed(accounting)
    else:
        assert module.completed(accounting)["JobID"] == "22153"


@pytest.mark.parametrize("problem", [None, "groups", "membership", "missing", "duplicate", "names"])
def test_independent_memberships(tmp_path, problem):
    left, right = tmp_path / "left", tmp_path / "right"
    left.write_text("a b\nc\n")
    right.write_text({"membership": "a\nb c\n", "missing": "a b\n",
                      "duplicate": "a b\na c\n"}.get(problem, "c\nb a\n"))
    names = ["a", "b", "c"] + (["a"] if problem == "names" else [])
    groups = 1 if problem == "groups" else 2
    if problem:
        with pytest.raises(ValueError):
            module.compare_partitions(left, right, names, groups)
    else:
        assert module.compare_partitions(left, right, names, groups)["partition_equal"] is True


@pytest.mark.parametrize("problem", [None, "source", "status", "job", "exit", "flags", "runtime", "worker",
    "command", "module", "output", "genes", "hits", "numeric", "auditor", "bool_exit"])
def test_report_bindings(tmp_path, monkeypatch, problem):
    root, executor, directory = tmp_path, tmp_path / "executor", tmp_path / "result"
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    def file(path):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("fixture")
        return record(path)
    source = file(executor / "benchmark_tools/reconstruct_cpm_refinement.py")
    monkeypatch.setattr(module, "SOURCE_SHA", source["sha256"])
    modules = [file(launcher / name) for name in (
        "benchmark_tools/replay_high_sensitivity.py", "orthohmm/accuracy.py", "orthohmm/refinement.py",
        "benchmark_tools/audit_historical_profile_ablation.py", "benchmark_tools/audit_accuracy_checkpoint.py")]
    numeric = {"status": "numeric_checkpoint_verified", "auditor": {**modules[-1], "path": "/original/auditor.py"}}
    child = {"modules": modules, "output": file(directory / "reconstructed.txt"), "accuracy_evaluated": False,
        "partition_equal": True, "genes": 984137, "refinement_directed_hits": 0,
        "numeric_checkpoint": {**numeric, "auditor": modules[-1]}}
    runtime = {"status": "verified"}
    parent = {"status": "original_cpm_refinement_reproduced_unscored", "job_id": "22153", "returncode": 0,
        "source": source, "expected": file(root / "benchmarks/results/qfo_parameter_cpm_replay_v3/cpm_high/replay/orthogroups_multipass_refined.txt"),
        "recovery_authorized": False, "accuracy_evaluated": False, "publication_ready": False,
        "runtime_before": runtime, "runtime_after": runtime, "worker_report": file(directory / "worker.json"),
        "worker_log": file(directory / "worker.log"), "result": copy.deepcopy(child),
        "command": [sys.executable, "-B", source["path"], "--root", str(root), "--output", str(directory), "--worker"]}
    if problem == "source":
        monkeypatch.setattr(module, "SOURCE_SHA", "changed")
    elif problem in ("status", "job", "exit", "flags", "runtime", "worker", "bool_exit"):
        key, value = {"status": ("status", "failed"), "job": ("job_id", "other"), "exit": ("returncode", 1),
            "flags": ("recovery_authorized", True), "runtime": ("runtime_after", {}),
            "worker": ("worker_report", {}), "bool_exit": ("returncode", False)}[problem]
        parent[key] = value
    elif problem == "command":
        parent["command"][-1] = "other"
    elif problem:
        if problem == "module":
            child["modules"] = []
        elif problem == "output":
            child["output"] = {}
        elif problem == "genes":
            child["genes"] = 976504
        elif problem == "hits":
            child["refinement_directed_hits"] = 1
        elif problem == "numeric":
            child["numeric_checkpoint"] = {}
        elif problem == "auditor":
            child["numeric_checkpoint"]["auditor"] = numeric["auditor"]
        parent["result"] = copy.deepcopy(child)
    if problem:
        with pytest.raises(ValueError):
            module.check_report(parent, child, root, executor, directory, runtime, numeric)
    else:
        assert len(module.check_report(parent, child, root, executor, directory, runtime, numeric)) == 10


def test_pending_job_rejected_before_outputs(tmp_path, monkeypatch):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
        "JobID|State|ExitCode|NodeList|AllocCPUS\n22153|PENDING|0:0|None assigned|0\n")
    with pytest.raises(ValueError, match="completed single-CPU"):
        module.validate(tmp_path)
