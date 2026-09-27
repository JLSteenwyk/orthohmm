import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import install_recovery_environment as module


def test_offline_hashed_commands():
    rows = module.commands(*map(Path, ("/base", "/installer", "/wheels", "/lock", "/out")))
    assert rows[0] == ["/base", "-I", "-m", "venv", "--without-pip", "/out/venv"]
    for option in ("--isolated", "--no-index", "--require-hashes", "--only-binary=:all:", "--no-cache-dir"):
        assert option in rows[1]
    assert rows[2] == ["/out/venv/bin/python", "-I", "-m", "pip", "check"]


def test_existing_destination_rejected(tmp_path):
    with pytest.raises(FileExistsError):
        module.install("/base", "/installer", "/wheels", "/lock", tmp_path)


def test_wrong_lock_rejected_before_creation(tmp_path):
    lock = tmp_path / "lock"
    lock.write_text("wrong")
    out = tmp_path / "out"
    with pytest.raises(ValueError):
        module.install("/base", "/installer", tmp_path, lock, out)
    assert not out.exists()


def test_missing_wheels_rejected(tmp_path, monkeypatch):
    lock = tmp_path / "lock"
    lock.write_text("fixture")
    monkeypatch.setattr(module, "LOCK_SHA", module.record(lock)["sha256"])
    with pytest.raises(ValueError):
        module.install("/base", "/installer", tmp_path, lock, tmp_path / "out")


@pytest.mark.parametrize("timeout", [False, True])
def test_failed_command_preserved_without_retry(tmp_path, monkeypatch, timeout):
    lock = tmp_path / "lock"
    lock.write_text("fixture")
    monkeypatch.setattr(module, "LOCK_SHA", module.record(lock)["sha256"])
    for i in range(11):
        (tmp_path / f"{i}.whl").write_bytes(b"fixture")
    calls = []

    def run(command, **kwargs):
        calls.append(command)
        if timeout:
            raise module.subprocess.TimeoutExpired(command, 600)
        return SimpleNamespace(returncode=17)

    monkeypatch.setattr(module.subprocess, "run", run)
    out = tmp_path / "out"
    with pytest.raises((RuntimeError, module.subprocess.TimeoutExpired)):
        module.install(lock, lock, tmp_path, lock, out)
    assert len(calls) == 1
    failure = json.loads((out / "failure.json").read_text())
    assert failure["retry"] is False
    assert not (out / "result.json").exists()
