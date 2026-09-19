import json
import shlex
import subprocess
from types import SimpleNamespace

import pytest

from benchmark_tools import submit_root_context_session as module


def test_bounded_wait_command_has_no_remote_polling_or_persistent_changes():
    value = module.command("a"*64)
    remote = shlex.split(value[-1])
    assert value[:3] == ["ssh", "-T", "-o"]
    assert remote[:7] == ["timeout", "--signal=TERM", "--kill-after=10s", "1020s", "sbatch", "--wait", "--parsable"]
    assert remote[-1] == "a"*64
    assert "recipe_v2" in remote[-2]
    assert not any(word in remote for word in ("loginctl", "systemctl", "squeue", "sleep"))
    with pytest.raises(ValueError):
        module.command("; unexpected")


@pytest.mark.parametrize("text,expected", [("123\n", 123), ("123;cluster\n", 123), ("0", None), ("warning\n123", None), ("", None)])
def test_unambiguous_job_identifier(text, expected):
    assert module.job_id(text) == expected


@pytest.mark.parametrize("queue_code,queue_text", [(0, ""), (0, "999 RUNNING other\n"), (1, "")])
def test_queue_guard_and_single_submission(tmp_path, monkeypatch, queue_code, queue_text):
    calls = []

    def run(command, **kwargs):
        calls.append((command, kwargs))
        if len(calls) == 1:
            return SimpleNamespace(returncode=queue_code, stdout=queue_text, stderr="")
        return SimpleNamespace(returncode=0, stdout="123\n", stderr="")

    monkeypatch.setattr(module.subprocess, "run", run)
    output = tmp_path / "launch"
    result = module.run(output, "a"*64)
    assert result["scheduler_terminal_verified"] is False
    assert json.loads((output / "result.json").read_text()) == result
    if queue_code == 0 and not queue_text:
        assert len(calls) == 2 and result["job_id"] == 123
        assert calls[1][1]["timeout"] == 1050
        assert (output / "launch.json").exists()
    else:
        assert len(calls) == 1 and result["status"] == "not_submitted_queue_not_verified_empty"
        assert not (output / "launch.json").exists()


@pytest.mark.parametrize("failure", ["timeout", "oserror", "remote_timeout", "job_failed"])
def test_observation_failure_never_retries_or_claims_terminal(tmp_path, monkeypatch, failure):
    calls = []

    def run(command, **kwargs):
        calls.append(command)
        if len(calls) == 1:
            return SimpleNamespace(returncode=0, stdout="", stderr="")
        if failure == "timeout":
            raise subprocess.TimeoutExpired(command, 1050, output=b"123\n")
        if failure == "oserror":
            raise OSError("connection failed")
        return SimpleNamespace(returncode=124 if failure == "remote_timeout" else 1, stdout="123\n", stderr="")

    monkeypatch.setattr(module.subprocess, "run", run)
    result = module.run(tmp_path / "launch", "a"*64)
    assert len(calls) == 2
    assert result["scheduler_terminal_verified"] is False
    assert result["scientific_timings_admitted"] is False
    if failure == "timeout":
        assert result["status"] == "observation_timeout" and result["job_id"] == 123
    elif failure == "oserror":
        assert result["status"] == "observation_error"
    else:
        assert result["status"] == "wait_returned" and result["returncode"] != 0


def test_existing_receipts_not_overwritten(tmp_path):
    with pytest.raises(FileExistsError):
        module.run(tmp_path, "a"*64)
