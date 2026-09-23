from types import SimpleNamespace
import json
import subprocess

import pytest

from benchmark_tools import observe_dgx_environment as module


@pytest.mark.parametrize("fault", [None, "timeout", "missing", "nonzero"])
def test_failures_are_retained_without_claiming_isolation(monkeypatch, fault):
    monkeypatch.setattr(module.os, "uname", lambda: SimpleNamespace(nodename="spark-7ff0"))
    calls = []
    def run(argv, **kwargs):
        assert kwargs["timeout"] == 10 and kwargs["env"]["LC_ALL"] == "C"
        calls.append(argv)
        if fault == "timeout": raise subprocess.TimeoutExpired(argv, 10)
        if fault == "missing": raise FileNotFoundError("fixture")
        return SimpleNamespace(returncode=7 if fault == "nonzero" else 0, stdout="[N/A]", stderr="")
    monkeypatch.setattr(module.subprocess, "run", run)
    result = module.collect(static=True)
    assert len(calls) == len(module.COMMANDS) + len(module.STATIC_COMMANDS)
    assert result["environmental_validity_established"] is False
    for row in result["commands"].values():
        if fault in {"timeout", "missing"}: assert "error_type" in row
        else:
            assert row["returncode"] == (7 if fault else 0)
            assert row["stdout"] == "[N/A]"
    assert "/proc/diskstats" in result["files"]


def test_other_host_rejected(monkeypatch):
    monkeypatch.setattr(module.os, "uname", lambda: SimpleNamespace(nodename="other"))
    with pytest.raises(ValueError): module.collect()


def test_session_recorder_cadence_stages_and_no_overwrite(tmp_path, monkeypatch):
    clock = [1.]
    monkeypatch.setattr(module.time, "monotonic", lambda: clock[0])
    monkeypatch.setattr(module, "collect", lambda static: {"static": static})
    recorder = module.Recorder(tmp_path / "environment")
    recorder.observe("before_submission", static=True)
    recorder.periodic()
    assert recorder.index == 1
    clock[0] = 31.
    recorder.periodic()
    recorder.observe("after_terminal_restoration", static=True)
    records = [json.loads(p.read_text()) for p in sorted(recorder.directory.glob("*.json"))]
    assert [r["session_stage"] for r in records] == ["before_submission", "job_wait", "after_terminal_restoration"]
    assert [r["index"] for r in records] == [0, 1, 2]
    with pytest.raises(FileExistsError):
        module.Recorder(recorder.directory).observe("before_submission")


def test_session_observation_failure_propagates(tmp_path, monkeypatch):
    def fail(static): raise OSError("fixture")
    monkeypatch.setattr(module, "collect", fail)
    recorder = module.Recorder(tmp_path / "environment")
    with pytest.raises(OSError): recorder.observe("before_submission")
    assert recorder.index == 0
