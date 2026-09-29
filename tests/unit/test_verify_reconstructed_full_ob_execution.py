import copy

import pytest

from benchmark_tools import verify_reconstructed_full_ob_execution as module


def accounting(state="COMPLETED", cpus="32", node="bizon", memory="128G", exit_code="0:0"):
    return ("JobIDRaw|State|ExitCode|Elapsed|AllocCPUS|NodeList|ReqMem\n"
            f"22376|{state}|{exit_code}|01:00:00|{cpus}|{node}|{memory}\n")


def test_correct_scheduler():
    assert module.scheduler_gate(accounting(), 22376)["State"] == "COMPLETED"


@pytest.mark.parametrize("changes", [dict(state="RUNNING"), dict(state="FAILED"),
    dict(exit_code="1:0"), dict(cpus="16"), dict(node="other"), dict(memory="128Gc")])
def test_wrong_scheduler(changes):
    with pytest.raises(ValueError):
        module.scheduler_gate(accounting(**changes), 22376)


def test_failed_batch_step():
    text = accounting() + "22376.batch|FAILED|1:0|00:10:00|32|bizon|128G\n"
    with pytest.raises(ValueError, match="scheduler step"):
        module.scheduler_gate(text, 22376)


def test_wrong_job():
    with pytest.raises(ValueError, match="prespecified"):
        module.verify(None, 22337)


def test_live_job_does_not_read_plan(tmp_path, monkeypatch):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: accounting(state="RUNNING"))
    with pytest.raises(ValueError, match="COMPLETED"):
        module.verify(tmp_path)
    assert list(tmp_path.iterdir()) == []


def binding(tmp_path):
    (tmp_path / "plan.json").write_text("{}")
    pin = module.record(tmp_path / "plan.json")
    plan = dict(command=["python", "workflow.py"])
    submission = dict(job_id=22376, protocol_commit="83f0ebea", plan=pin)
    start = dict(plan=pin, job_id="22376")
    execution = dict(status="integrated_complete_pending_independent_admission",
        returncode=0, job_id="22376", plan=pin,
        command=["/usr/bin/time", "-v", "-o", str(tmp_path / "time.txt"), *plan["command"]])
    return plan, submission, start, execution


def test_execution_binding(tmp_path):
    module.execution_binding(tmp_path, *binding(tmp_path))


@pytest.mark.parametrize("index,key,value", [
    (1, "job_id", 22337), (1, "protocol_commit", "other"),
    (2, "job_id", "22337"), (3, "returncode", 1),
    (3, "status", "running"), (3, "command", ["different"]),
    (3, "plan", {})])
def test_tampered_binding(tmp_path, index, key, value):
    values = copy.deepcopy(binding(tmp_path))
    values[index][key] = value
    with pytest.raises(ValueError, match="binding differs"):
        module.execution_binding(tmp_path, *values)
