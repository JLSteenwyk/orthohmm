from types import SimpleNamespace
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
