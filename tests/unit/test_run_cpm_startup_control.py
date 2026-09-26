import json
from pathlib import Path
import subprocess
from types import SimpleNamespace

import pytest

import benchmark_tools.run_cpm_startup_control as module


def setup(tmp_path, monkeypatch):
    python = tmp_path / "python"
    status = tmp_path / "status.json"
    def record(path):
        return dict(path=str(path), bytes=1, sha256="fixture")
    source = dict(scientific_child_command=[str(python), "-B", str(tmp_path / "tools/run.py")],
                  checked_records=[record(python), record("/usr/lib/x86_64-linux-gnu/libc.so.6")],
                  memchecker=dict(binary=record("/usr/bin/valgrind"), components=[]))
    status.write_text(json.dumps(source))
    monkeypatch.setattr(module, "STATUS_SHA", "fixture")
    monkeypatch.setattr(module, "digest", lambda path: "fixture")
    monkeypatch.setattr(module, "record", record)
    monkeypatch.setattr(module, "summarize", lambda *args: dict(error_records_by_kind={"InvalidRead": 1}))
    return status, tmp_path / "output", python


@pytest.mark.parametrize("returncode", [0, 97])
def test_one_minimal_child_and_no_admission(tmp_path, monkeypatch, returncode):
    status, output, python = setup(tmp_path, monkeypatch)
    calls = []
    def run(command, **kwargs):
        calls.append(command)
        assert command[-5:] == [str(python), "-B", "-S", "-c", "pass"]
        assert kwargs["timeout"] == 120
        assert kwargs["env"]["PYTHONMALLOC"] == "malloc"
        return SimpleNamespace(returncode=returncode)
    monkeypatch.setattr(module.subprocess, "run", run)
    report = module.run(status, output)
    assert report["returncode"] == returncode
    assert report["attempts"] == len(calls) == 1
    assert not report["accuracy_admitted"] and not report["publication_ready"]
    assert json.loads((output / "report.json").read_text()) == report
    with pytest.raises(FileExistsError):
        module.run(status, output)
    assert len(calls) == 1


def test_timeout_retains_failed_report_without_retry(tmp_path, monkeypatch):
    status, output, python = setup(tmp_path, monkeypatch)
    calls = []
    def timeout(command, **kwargs):
        calls.append(command)
        raise subprocess.TimeoutExpired(command, 120)
    monkeypatch.setattr(module.subprocess, "run", timeout)
    with pytest.raises(subprocess.TimeoutExpired):
        module.run(status, output)
    report = json.loads((output / "report.json").read_text())
    assert report["status"] == "startup_control_failed"
    assert len(calls) == 1


def test_changed_status_rejected_before_execution(tmp_path):
    status = tmp_path / "status.json"
    status.write_text("{}")
    output = tmp_path / "output"
    with pytest.raises(ValueError, match="Changed retained"):
        module.run(status, output)
    assert not output.exists()
