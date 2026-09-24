import copy
from pathlib import Path
import subprocess
import sys

import pytest

from benchmark_tools import admit_cpm_checkpoint_recovery as module


@pytest.mark.parametrize("problem", [None, "running", "failed", "exit", "cpu", "memory", "node", "missing", "duplicate", "job"])
def test_scheduler(problem):
    fields = ["123", "COMPLETED", "0:0", "1", "64G", "bizon"]
    for label, index, value in (("running", 1, "RUNNING"), ("failed", 1, "FAILED"), ("exit", 2, "1:0"),
        ("cpu", 3, "2"), ("memory", 4, "32G"), ("node", 5, "other")):
        if problem == label:
            fields[index] = value
    rows = ["|".join(fields)]
    if problem == "missing":
        rows = []
    elif problem == "duplicate":
        rows *= 2
    text = "JobID|State|ExitCode|AllocCPUS|ReqMem|NodeList\n" + "\n".join(rows)
    if problem:
        with pytest.raises(ValueError):
            module.completed(text, "123_0" if problem == "job" else "123")
    else:
        assert module.completed(text, "123")["State"] == "COMPLETED"


@pytest.mark.parametrize("problem", [None, "status", "job", "commit", "source", "attempt", "accuracy", "admitted",
    "invented_stats", "phase_order", "extra_optimizer", "phase_exit", "bool_exit", "command", "nan", "negative",
    "log", "stage_origin", "stage_output"])
def test_parent_contract(tmp_path, monkeypatch, problem):
    def record(path):
        return dict(path=str(path), bytes=7, sha256="fixture")
    monkeypatch.setattr(module, "record", record)
    root, directory = tmp_path, tmp_path / "recovery"
    source = record(tmp_path / "runner.py")
    phases = [dict(mode=mode, command=[sys.executable, "-B", source["path"], "--root", str(root),
        "--output", str(directory), "--mode", mode], status="completed", returncode=0, wall_s=1.,
        log=record(directory / f"{mode}.log")) for mode in ("optimize", "refine", "repeat-refinement")]
    original = root / "benchmarks/results/qfo_parameter_cpm_replay_v3/cpm_high/replay"
    stages = [dict(label=label, origin=origin, output=record(path)) for label, origin, path in (
        ("multipass", "reused", original / "orthogroups_multipass.txt"),
        ("multipass_refined", "reused", original / "orthogroups_multipass_refined.txt"),
        ("strict_profiles", "recovered", directory / "orthogroups_profiles.txt"),
        ("strict_profiles_refined", "recovered", directory / "orthogroups_profiles_refined.txt"))]
    parent = dict(status="cpm_checkpoint_recovered_pending_independent_admission", job_id="123", executor_commit="commit",
        source=source, optimizer_attempted=True, accuracy_evaluated=False, downstream_admitted=False, publication_ready=False,
        missing_original_statistics=["profile_counters", "successful_stage_timings", "full_run_time"], phases=phases, stages=stages)
    if problem in ("status", "job", "commit", "source", "attempt", "accuracy", "admitted", "invented_stats"):
        key, value = {"status": ("status", "failed"), "job": ("job_id", "999"), "commit": ("executor_commit", "changed"),
            "source": ("source", {}), "attempt": ("optimizer_attempted", False), "accuracy": ("accuracy_evaluated", True),
            "admitted": ("downstream_admitted", True), "invented_stats": ("missing_original_statistics", [])}[problem]
        parent[key] = value
    elif problem == "phase_order":
        parent["phases"].reverse()
    elif problem == "extra_optimizer":
        parent["phases"].append(copy.deepcopy(phases[0]))
    elif problem in ("phase_exit", "bool_exit", "command", "nan", "negative", "log"):
        key, value = {"phase_exit": ("returncode", 1), "bool_exit": ("returncode", False), "command": ("command", []),
            "nan": ("wall_s", float("nan")), "negative": ("wall_s", -1), "log": ("log", {})}[problem]
        parent["phases"][0][key] = value
    elif problem == "stage_origin":
        parent["stages"][0]["origin"] = "recovered"
    elif problem == "stage_output":
        parent["stages"][3]["output"] = stages[2]["output"]
    if problem:
        with pytest.raises(ValueError):
            module.parent_contract(parent, "123", "commit", source, root, directory)
    else:
        assert module.parent_contract(parent, "123", "commit", source, root, directory) == stages


def test_pending_recovery_creates_no_admission_output(tmp_path, monkeypatch):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
        "JobID|State|ExitCode|AllocCPUS|ReqMem|NodeList\n123|PENDING|0:0|0|64G|None assigned\n")
    with pytest.raises(ValueError, match="completed one-CPU"):
        module.admit(tmp_path, "123", "commit", tmp_path / "output")
    assert not (tmp_path / "output").exists()


def test_existing_output_is_preserved(tmp_path):
    with pytest.raises(FileExistsError):
        module.admit(tmp_path, "123", "commit", tmp_path)


def test_isolated_cli_help(tmp_path):
    done = subprocess.run([sys.executable, "-I", str(Path(module.__file__).resolve()), "--help"],
                          cwd=tmp_path, capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
