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


def unit(path, dropins=""):
    return f"Id=test.service\nLoadState=loaded\nFragmentPath={path}\nDropInPaths={dropins}\nNeedDaemonReload=no\n"


def test_fingerprints_detect_content_and_symlink_changes_without_copying_contents(tmp_path):
    source = tmp_path / "test.service"
    source.write_text("private configuration")
    link = tmp_path / "alias.service"
    link.symlink_to(source)
    dropin = tmp_path / "override.conf"
    dropin.write_text("override")
    first = module.configuration_fingerprints(unit(link, str(dropin)))
    assert "private configuration" not in json.dumps(first)
    assert first["test.service"]["files"][0]["symlink_target"] == str(source)
    source.write_text("changed")
    second = module.configuration_fingerprints(unit(link, str(dropin)))
    assert first["test.service"]["files"][0]["sha256"] != second["test.service"]["files"][0]["sha256"]
    assert first["test.service"]["files"][1] == second["test.service"]["files"][1]


@pytest.mark.parametrize("raw", ["", "Id=test.service", "Id=x\nId=y", unit("/tmp/a") + "\n" + unit("/tmp/a")])
def test_malformed_configuration_rejected(raw):
    with pytest.raises(ValueError): module.configuration_fingerprints(raw)


@pytest.mark.parametrize("path", ["/nonexistent/unit", "/dev/null", "relative", r"/tmp/a\x20b"])
def test_unreadable_or_unsupported_paths_remain_errors(path):
    result = module.configuration_fingerprints(unit(path))
    assert "error_type" in result["test.service"]["files"][0]
