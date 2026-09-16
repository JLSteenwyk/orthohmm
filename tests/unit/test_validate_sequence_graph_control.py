from pathlib import Path

import pytest

from benchmark_tools import validate_sequence_graph_control as validator


def accounting(state="COMPLETED", code="0:0"):
    return ("JobID|JobIDRaw|State|ExitCode|Elapsed\n"
            f"21293_0|21296|{state}|{code}|00:01:00\n"
            "21293_1|21293|COMPLETED|0:0|00:01:00\n")


def test_array_identity_retains_raw_job_mapping():
    tasks = validator.completed_tasks(accounting())
    assert tasks["all_hits"]["JobIDRaw"] == "21296"
    assert tasks["top100"]["JobIDRaw"] == "21293"


@pytest.mark.parametrize("state,code", [("RUNNING", "0:0"), ("PENDING", "0:0"),
                                       ("FAILED", "1:0"), ("COMPLETED", "1:0")])
def test_rejects_partial_or_failed_execution(state, code):
    with pytest.raises(ValueError):
        validator.completed_tasks(accounting(state, code))


@pytest.mark.parametrize("duplicate", [False, True])
def test_rejects_incomplete_or_duplicate_tasks(duplicate):
    text = accounting()
    text = text + "21293_0|21296|COMPLETED|0:0|00:01:00\n" if duplicate else text.replace(
        "21293_1|21293|COMPLETED|0:0|00:01:00\n", "")
    with pytest.raises(ValueError):
        validator.completed_tasks(text)


@pytest.mark.parametrize("change", [None, "job", "variant", "command", "environment", "failed", "scored", "profile"])
def test_native_record_admission(change):
    command = ["python", "replay.py", "--cpu", "32"]
    preflight = {"variant": "all_hits", "job_id": "21296", "accuracy_evaluated": False,
                 "command": command, "environment_overrides": {
                     "PYTHONPATH": "/launcher", "PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1",
                     "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}}
    result = {"status": "graph_complete_pending_scoring", "exit_code": 0, "accuracy_evaluated": False}
    metrics = {"command": command, "cwd": "/launcher", "parameters": {"profile_expansion": False},
               "counts": {"genes": 251378, "species": 12},
               "stages": [{"label": "multipass"}, {"label": "multipass_refined"}]}
    if change == "job":
        preflight["job_id"] = "21293"
    elif change == "variant":
        preflight["variant"] = "top100"
    elif change == "command":
        metrics["command"] = command + ["--fasta-directory", "/input"]
    elif change == "environment":
        preflight["environment_overrides"]["PYTHONHASHSEED"] = "1"
    elif change == "failed":
        result["status"] = "verification_failed"
    elif change == "scored":
        result["accuracy_evaluated"] = True
    elif change == "profile":
        metrics["parameters"]["profile_expansion"] = True
    args = preflight, result, metrics, "all_hits", {"JobIDRaw": "21296"}, command, Path("/launcher")
    if change:
        with pytest.raises(ValueError):
            validator.check_records(*args)
    else:
        validator.check_records(*args)


def test_live_tasks_cannot_reach_runtime_or_outputs(tmp_path, monkeypatch):
    monkeypatch.setattr(validator.subprocess, "check_output", lambda *a, **k: accounting("RUNNING"))
    monkeypatch.setattr(validator, "verify", lambda *a: pytest.fail("Reached native admission before completion"))
    with pytest.raises(ValueError):
        validator.validate(tmp_path)
