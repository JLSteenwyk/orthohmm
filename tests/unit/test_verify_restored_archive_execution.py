import copy

import pytest

from benchmark_tools import verify_restored_archive_execution as module


def accounting(state="COMPLETED", cpus="32", node="bizon", memory="128G", exit_code="0:0"):
    header = "JobIDRaw|State|ExitCode|Elapsed|AllocCPUS|NodeList|ReqMem\n"
    return header + "".join(
        f"{job}|{state}|{exit_code}|01:00:00|{cpus}|{node}|{memory}\n"
        for job in ("22377", "22377.batch"))


def test_completed_scheduler():
    assert module.scheduler_gate(accounting())["State"] == "COMPLETED"


@pytest.mark.parametrize("changes", [dict(state="RUNNING"), dict(state="FAILED"),
    dict(exit_code="1:0"), dict(cpus="16"), dict(node="other"), dict(memory="128Gc")])
def test_wrong_scheduler(changes):
    with pytest.raises(ValueError):
        module.scheduler_gate(accounting(**changes))


def test_missing_batch():
    with pytest.raises(ValueError, match="Missing batch"):
        module.scheduler_gate("\n".join(accounting().splitlines()[:2]) + "\n")


def test_failed_step():
    with pytest.raises(ValueError, match="scheduler step"):
        module.scheduler_gate(accounting() + "22377.0|FAILED|1:0|00:01:00|32|bizon|128G\n")


def test_live_job_does_not_read_plan(tmp_path, monkeypatch):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: accounting(state="RUNNING"))
    with pytest.raises(ValueError, match="COMPLETED"):
        module.verify(tmp_path)
    assert list(tmp_path.iterdir()) == []


def binding(tmp_path):
    (tmp_path / "plan.json").write_text("{}")
    pin = module.record(tmp_path / "plan.json")
    plan = dict(command=["python", "workflow.py"])
    submission = dict(job_id="22377", stdout="22377\n", plan=pin)
    start = dict(plan=pin, job_id="22377")
    execution = dict(status="integrated_complete_pending_independent_admission",
        returncode=0, job_id="22377", plan=pin,
        command=["/usr/bin/time", "-v", "-o", str(tmp_path / "time.txt"), *plan["command"]])
    return plan, submission, start, execution


def test_execution_binding(tmp_path):
    module.execution_binding(tmp_path, *binding(tmp_path))


@pytest.mark.parametrize("index,key,value", [
    (1, "job_id", "22376"), (1, "stdout", "22376\n"),
    (2, "job_id", "22376"), (3, "returncode", 1),
    (3, "status", "running"), (3, "command", ["different"]),
    (3, "plan", {})])
def test_tampered_binding(tmp_path, index, key, value):
    values = copy.deepcopy(binding(tmp_path))
    values[index][key] = value
    with pytest.raises(ValueError, match="binding differs"):
        module.execution_binding(tmp_path, *values)
