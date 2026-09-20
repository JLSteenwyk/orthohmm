from datetime import datetime
import json
import shlex
import subprocess
from types import SimpleNamespace
from zoneinfo import ZoneInfo

import pytest

from benchmark_tools import submit_root_context_overhead_session as submit
from benchmark_tools import audit_root_context_overhead_session as audit
from tests.unit.test_prepare_root_context_overhead import module as prepare, PARENT, PROTOCOL


def test_wait_bounds_match_frozen_plan_and_command():
    plan = prepare.build(PARENT, PROTOCOL)
    assert plan["waiting_session"] == dict(local_timeout_s=submit.LOCAL_TIMEOUT,
        remote_timeout_s=submit.REMOTE_TIMEOUT, remote_kill_grace_s=submit.KILL_GRACE)
    assert submit.LOCAL_TIMEOUT > submit.REMOTE_TIMEOUT+10 > plan["allocation"]["time_limit_s"]
    remote = shlex.split(submit.command("a"*64)[-1])
    assert remote[:7] == ["timeout", "--signal=TERM", "--kill-after=10s", "18120s", "sbatch", "--wait", "--parsable"]
    assert remote[-2] == str(submit.SUBMISSION_SCRIPT)
    assert not any(word in remote for word in ("systemctl", "loginctl", "squeue"))
    with pytest.raises(ValueError):
        submit.command(";invalid")


@pytest.mark.parametrize("mode", ["success", "busy", "queue_error", "queue_timeout", "timeout", "oserror", "remote_timeout", "job_failed"])
def test_single_submission_retains_all_observation_outcomes(tmp_path, monkeypatch, mode):
    calls = []

    def run(command, **kwargs):
        calls.append(command)
        if len(calls) == 1:
            assert command == submit.QUEUE_COMMAND and kwargs["timeout"] == 15
            if mode == "queue_timeout":
                raise subprocess.TimeoutExpired(command, 15)
            return SimpleNamespace(returncode=1 if mode == "queue_error" else 0,
                stdout="9 RUNNING other" if mode == "busy" else "", stderr="")
        assert kwargs["timeout"] == 18150
        launch = json.loads((tmp_path / "receipt/launch.json").read_text())
        assert launch["remote_timeout_s"] == 18120 and launch["local_timeout_s"] == 18150
        if mode == "timeout":
            raise subprocess.TimeoutExpired(command, 18150, output=b"123\n")
        if mode == "oserror":
            raise OSError("connection failed")
        return SimpleNamespace(returncode={"remote_timeout": 124, "job_failed": 1}.get(mode, 0), stdout="123\n", stderr="")

    monkeypatch.setattr(submit.subprocess, "run", run)
    result = submit.run(tmp_path / "receipt", "a"*64)
    assert result == json.loads((tmp_path / "receipt/result.json").read_text())
    assert result["scientific_timings_admitted"] is False and result["scheduler_terminal_verified"] is False
    assert len(calls) == (1 if mode in {"busy", "queue_error", "queue_timeout"} else 2)
    if mode == "timeout":
        assert result["status"] == "observation_timeout" and result["job_id"] == 123
    elif mode in {"success", "remote_timeout", "job_failed"}:
        assert result["status"] == "wait_returned"
        assert result["returncode"] == {"success": 0, "remote_timeout": 124, "job_failed": 1}[mode]
    elif mode == "oserror":
        assert result["status"] == "observation_error"
    with pytest.raises(FileExistsError):
        submit.run(tmp_path / "receipt", "a"*64)


@pytest.fixture
def evidence():
    start = int(datetime(2026, 9, 19, 20, 0, tzinfo=ZoneInfo("America/New_York")).timestamp())*10**9
    queue = dict(command=submit.QUEUE_COMMAND, returncode=0, stdout="", started_unix_ns=start-2*10**9,
                 finished_unix_ns=start-10**9)
    launch = dict(command=submit.command("a"*64), recipe_sha256="a"*64, local_timeout_s=18150,
        remote_timeout_s=18120, remote_kill_grace_s=10, scientific_timings_admitted=False,
        started_unix_ns=start-500000000)
    result = dict(source=dict(sha256="source", bytes=100), status="wait_returned", returncode=0, job_id=123,
        stdout="123\n", recipe_sha256="a"*64, scientific_timings_admitted=False,
        scheduler_terminal_verified=False, finished_unix_ns=start+18001*10**9)
    recipe = dict(records=[dict(kind="file", path=str(submit.RECIPE_ROOT / "benchmark_tools/submit_root_context_overhead_session.py"),
        sha256="source", bytes=100)])
    allocation = dict(JobId="123", ExitCode="0:0", StartTime="2026-09-19T20:00:00", EndTime="2026-09-20T01:00:00")
    return [queue, launch, result, allocation, recipe, "a"*64, 123, "America/New_York"]


@pytest.mark.parametrize("exit_code,returned", [("0:0", 0), ("1:0", 1), ("0:9", 1)])
def test_full_five_hour_receipt_and_terminal_outcomes(evidence, exit_code, returned):
    evidence[3]["ExitCode"] = exit_code
    evidence[2]["returncode"] = returned
    result = audit.validate(*evidence)
    assert result["status"] == "overhead_bounded_session_verified"
    assert result["scientific_timings_admitted"] is False


@pytest.mark.parametrize("fault", ["source", "duplicate", "queue", "queue_bool", "command", "short_bound",
    "job", "job_bool", "stdout", "exit_bool", "exit", "timeout", "early", "late", "overbound", "clock",
    "clock_bool", "admission", "offset"])
def test_inconsistent_receipt_rejected(evidence, fault):
    queue, launch, result, allocation, recipe, *_ = evidence
    if fault == "source":
        result["source"]["sha256"] = "other"
    elif fault == "duplicate":
        recipe["records"].append(recipe["records"][0])
    elif fault == "queue":
        queue["stdout"] = "123 RUNNING other"
    elif fault == "queue_bool":
        queue["returncode"] = False
    elif fault == "command":
        launch["command"] = []
    elif fault == "short_bound":
        launch["remote_timeout_s"] = 3720
    elif fault == "job":
        result["job_id"] = 124
    elif fault == "job_bool":
        evidence[-2] = True
    elif fault == "stdout":
        result["stdout"] = "warning\n123\n"
    elif fault == "exit_bool":
        result["returncode"] = False
    elif fault == "exit":
        result["returncode"] = 124
    elif fault == "timeout":
        result["status"] = "observation_timeout"
    elif fault == "early":
        result["finished_unix_ns"] -= 5*10**9
    elif fault == "late":
        launch["started_unix_ns"] += 5*10**9
    elif fault == "overbound":
        result["finished_unix_ns"] += 200*10**9
    elif fault == "clock":
        queue["finished_unix_ns"] = result["finished_unix_ns"]
    elif fault == "clock_bool":
        queue["started_unix_ns"] = True
    elif fault == "admission":
        result["scientific_timings_admitted"] = 0
    else:
        allocation["StartTime"] += "+00:00"
    with pytest.raises(ValueError):
        audit.validate(*evidence)
