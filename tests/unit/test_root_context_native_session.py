from datetime import datetime
import json
import shlex
import subprocess
from types import SimpleNamespace
from zoneinfo import ZoneInfo

import pytest

from benchmark_tools import submit_root_context_native_session as submit
from benchmark_tools import audit_root_context_native_session as audit


def test_one_hour_wait_command():
    command = submit.command("a"*64)
    remote = shlex.split(command[-1])
    assert remote[:7] == ["timeout", "--signal=TERM", "--kill-after=10s", "3720s", "sbatch", "--wait", "--parsable"]
    assert remote[-2] == str(submit.SUBMISSION_SCRIPT)
    assert submit.LOCAL_TIMEOUT > submit.REMOTE_TIMEOUT + submit.KILL_GRACE > 3600
    assert not any(word in remote for word in ("systemctl", "loginctl", "squeue"))
    with pytest.raises(ValueError):
        submit.command("a;unexpected")


@pytest.mark.parametrize("mode", ["success", "busy", "queue_error", "queue_timeout", "timeout", "oserror", "remote_timeout", "job_failed"])
def test_submit_once_and_retain_observation_failures(tmp_path, monkeypatch, mode):
    calls = []

    def run(command, **kwargs):
        calls.append(command)
        if len(calls) == 1:
            assert command == submit.QUEUE_COMMAND and kwargs["timeout"] == 15
            if mode == "queue_timeout":
                raise subprocess.TimeoutExpired(command, 15)
            return SimpleNamespace(returncode=1 if mode == "queue_error" else 0,
                stdout="99 RUNNING other" if mode == "busy" else "", stderr="")
        assert kwargs["timeout"] == 3750
        launch = json.loads((tmp_path / "receipt/launch.json").read_text())
        assert launch["remote_timeout_s"] == 3720 and launch["local_timeout_s"] == 3750
        if mode == "timeout":
            raise subprocess.TimeoutExpired(command, 3750, output=b"123\n")
        if mode == "oserror":
            raise OSError("connection failed")
        return SimpleNamespace(returncode={"remote_timeout": 124, "job_failed": 1}.get(mode, 0), stdout="123\n", stderr="")

    monkeypatch.setattr(submit.subprocess, "run", run)
    result = submit.run(tmp_path / "receipt", "a"*64)
    assert result == json.loads((tmp_path / "receipt/result.json").read_text())
    assert result["scientific_timings_admitted"] is False and result["scheduler_terminal_verified"] is False
    assert len(calls) == (1 if mode in {"busy", "queue_error", "queue_timeout"} else 2)
    if mode == "timeout":
        assert result["job_id"] == 123 and result["status"] == "observation_timeout"
    elif mode == "success":
        assert result["status"] == "wait_returned" and result["returncode"] == 0
    elif mode in {"remote_timeout", "job_failed"}:
        assert result["status"] == "wait_returned" and result["returncode"] != 0
    elif mode == "oserror":
        assert result["status"] == "observation_error"
    with pytest.raises(FileExistsError):
        submit.run(tmp_path / "receipt", "a"*64)


@pytest.fixture
def evidence():
    start = int(datetime(2026, 9, 19, 19, 0, tzinfo=ZoneInfo("America/New_York")).timestamp())*10**9
    queue = dict(command=submit.QUEUE_COMMAND, returncode=0, stdout="", started_unix_ns=start-2*10**9,
                 finished_unix_ns=start-10**9)
    launch = dict(command=submit.command("a"*64), recipe_sha256="a"*64, local_timeout_s=3750,
        remote_timeout_s=3720, remote_kill_grace_s=10, scientific_timings_admitted=False,
        started_unix_ns=start-500000000)
    result = dict(source=dict(sha256="source", bytes=100), status="wait_returned", returncode=0, job_id=123,
        stdout="123\n", recipe_sha256="a"*64, scientific_timings_admitted=False,
        scheduler_terminal_verified=False, finished_unix_ns=start+3601*10**9)
    recipe = dict(records=[dict(kind="file", path=str(submit.RECIPE_ROOT / "benchmark_tools/submit_root_context_native_session.py"),
        sha256="source", bytes=100)])
    allocation = dict(JobId="123", ExitCode="0:0", StartTime="2026-09-19T19:00:00", EndTime="2026-09-19T20:00:00")
    return [queue, launch, result, allocation, recipe, "a"*64, 123, "America/New_York"]


@pytest.mark.parametrize("exit_code,returned", [("0:0", 0), ("1:0", 1), ("0:9", 1)])
def test_receipt_accepts_terminal_success_failure_and_signal(evidence, exit_code, returned):
    evidence[3]["ExitCode"] = exit_code
    evidence[2]["returncode"] = returned
    result = audit.validate(*evidence)
    assert result["status"] == "native_bounded_session_verified"
    assert result["scientific_timings_admitted"] is False


@pytest.mark.parametrize("fault", ["source", "duplicate", "queue", "queue_bool", "command", "bound",
    "old_bound", "job", "job_bool", "stdout", "exit_bool", "exit", "timeout", "early", "late",
    "overbound", "clock", "clock_bool", "admission", "timezone_timestamp"])
def test_receipt_drift_rejected(evidence, fault):
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
    elif fault == "bound":
        launch["remote_timeout_s"] = 3600
    elif fault == "old_bound":
        launch["local_timeout_s"] = 1050
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
