from types import SimpleNamespace

import pytest

from benchmark_tools import inspect_orthofinder_runtime as inspector


def test_missing_tools(monkeypatch):
    monkeypatch.setattr(inspector.shutil, "which", lambda name, path: None)
    assert all(row["status"] == "missing" for row in inspector.inspect_tools("/empty").values())


def test_child_path_used_and_version_exit_retained(monkeypatch, tmp_path):
    executable = tmp_path / "tool"
    executable.write_text("fixture")
    monkeypatch.setattr(inspector.shutil, "which", lambda name, path: str(executable))
    calls = []
    def run(argv, **kwargs):
        calls.append(kwargs)
        return SimpleNamespace(returncode=1, stdout="", stderr="version on stderr")
    monkeypatch.setattr(inspector.subprocess, "run", run)
    result = inspector.inspect_tools("/child/path")
    assert all(call["env"]["PATH"] == "/child/path" for call in calls)
    assert result["FastTree"]["exit_code"] == 1
    assert result["FastTree"]["stderr"] == "version on stderr"


def test_binary_change_rejected(monkeypatch, tmp_path):
    executable = tmp_path / "tool"
    executable.write_text("before")
    monkeypatch.setattr(inspector.shutil, "which", lambda name, path: str(executable))
    def run(*args, **kwargs):
        executable.write_text("after")
        return SimpleNamespace(returncode=0, stdout="version", stderr="")
    monkeypatch.setattr(inspector.subprocess, "run", run)
    with pytest.raises(ValueError, match="changed"):
        inspector.inspect_tools("/child/path")
