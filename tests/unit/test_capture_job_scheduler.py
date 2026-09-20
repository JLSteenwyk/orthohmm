import pytest

from benchmark_tools.capture_job_scheduler import terminal_record
from benchmark_tools.capture_array_scheduler import REQUIRED
from benchmark_tools import capture_job_scheduler as module
from types import SimpleNamespace
import json


def raw(**changes):
    fields = {key: "1" for key in REQUIRED - {"ArrayJobId", "ArrayTaskId"}}
    fields.update(JobId="21912", JobState="COMPLETED", ExitCode="0:0", **changes)
    return " ".join(f"{k}={v}" for k,v in fields.items())


def test_exact_terminal_job_only():
    text = raw()
    assert terminal_record(text, 21912) == text + "\n"
    assert terminal_record(text, 21913) is None
    assert terminal_record(text.replace("COMPLETED", "RUNNING"), 21912) is None


@pytest.mark.parametrize("fault", ["duplicate", "array", "missing", "exit"])
def test_malformed_record(fault):
    text = raw()
    if fault == "duplicate":
        text += " JobId=21912"
    elif fault == "array":
        text += " ArrayJobId=21912"
    elif fault == "missing":
        text = text.replace("Command=1", "Other=1")
    else:
        text = text.replace("ExitCode=0:0", "ExitCode=bad")
    with pytest.raises(ValueError):
        terminal_record(text, 21912)


def test_capture_retains_all_exact_terminal_records(tmp_path, monkeypatch):
    def run(argv, **kwargs):
        job = argv[3]
        return SimpleNamespace(returncode=0, stderr="", stdout=raw().replace("JobId=21912", f"JobId={job}"))
    monkeypatch.setattr(module.subprocess, "run", run)
    result = module.capture([21912, 21913, 21914], tmp_path / "capture")
    assert result["status"] == "complete"
    assert result["polls"] == 1
    assert result["missing"] == []
    assert not result["scientific_timings_admitted"]
    for job in result["jobs"]:
        assert f"JobId={job}" in (tmp_path / "capture" / f"scheduler_{job}.txt").read_text()
    assert json.loads((tmp_path / "capture/capture.json").read_text())["status"] == "complete"


@pytest.mark.parametrize("limit", [0, -1, True, "86580", float("inf"), float("nan")])
def test_invalid_observation_limit_leaves_no_directory(tmp_path, limit):
    output = tmp_path / "capture"
    with pytest.raises(ValueError, match="duration"):
        module.capture([21912], output, max_seconds=limit)
    assert not output.exists()


@pytest.mark.parametrize("arguments,expected", [([], 14400), (["--max-seconds", "86580"], 86580)])
def test_cli_forwards_explicit_long_run_bound(tmp_path, monkeypatch, arguments, expected):
    calls = []
    def capture(jobs, output, max_seconds):
        calls.append((jobs, output, max_seconds))
        return {"status": "incomplete"}
    monkeypatch.setattr(module, "capture", capture)
    assert module.main(["--jobs", "21912", "--output", str(tmp_path), *arguments]) == 1
    assert calls == [([21912], tmp_path, expected)]


def clock(monkeypatch):
    now = [0.0]
    monkeypatch.setattr(module.time, "monotonic", lambda: now[0])
    monkeypatch.setattr(module.time, "sleep", lambda seconds: now.__setitem__(0, now[0] + seconds))
    return now


def test_capture_can_observe_beyond_original_four_hour_cli_bound(tmp_path, monkeypatch):
    now = clock(monkeypatch)
    calls = []
    def run(argv, **kwargs):
        calls.append(argv)
        if len(calls) == 1:
            now[0] = 14401
            return SimpleNamespace(returncode=0, stderr="", stdout=raw().replace("COMPLETED", "RUNNING"))
        return SimpleNamespace(returncode=0, stderr="", stdout=raw())
    monkeypatch.setattr(module.subprocess, "run", run)
    result = module.capture([21912], tmp_path / "capture", max_seconds=86580)
    assert result["status"] == "complete" and len(calls) == 2
    assert result["elapsed_s"] > 14400 and result["max_seconds"] == 86580
    assert result["observation_limit_reached"] is False


def test_timeout_retains_same_job_as_missing_without_retry_or_other_commands(tmp_path, monkeypatch):
    now = clock(monkeypatch)
    calls = []
    def run(argv, **kwargs):
        calls.append((argv, kwargs["timeout"]))
        now[0] += kwargs["timeout"]
        raise module.subprocess.TimeoutExpired(argv, kwargs["timeout"])
    monkeypatch.setattr(module.subprocess, "run", run)
    result = module.capture([21912, 21913], tmp_path / "capture", max_seconds=2)
    assert calls == [(["scontrol", "show", "job", "21912", "--oneliner"], 2)]
    assert result["status"] == "incomplete" and result["missing"] == [21912, 21913]
    assert result["observation_limit_reached"] is True
    assert result["observation_errors"] == 1 and result["elapsed_s"] == 2
    assert json.loads((tmp_path / "capture/poll_000000_21912.json").read_text())["error_type"] == "TimeoutExpired"


def test_final_sleep_is_limited_and_nonterminal_is_not_promoted(tmp_path, monkeypatch):
    clock(monkeypatch)
    calls = []
    def run(argv, **kwargs):
        calls.append(argv)
        return SimpleNamespace(returncode=0, stderr="", stdout=raw().replace("COMPLETED", "RUNNING"))
    monkeypatch.setattr(module.subprocess, "run", run)
    result = module.capture([21912], tmp_path / "capture", max_seconds=0.25)
    assert len(calls) == 1 and result["elapsed_s"] == 0.25
    assert result["status"] == "incomplete" and result["retained"] == {}
