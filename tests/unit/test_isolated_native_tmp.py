import os

import pytest

from benchmark_tools.isolated_native_tmp import fresh_tmp


def test_distinct_overrides_restored_after_failure(tmp_path, monkeypatch):
    monkeypatch.setenv("TMPDIR", "/shared")
    monkeypatch.setenv("TMP", "/other")
    monkeypatch.delenv("TEMP", raising=False)
    directory = tmp_path / "scratch"
    with pytest.raises(RuntimeError):
        with fresh_tmp(directory):
            assert all(os.environ[k] == str(directory) for k in ("TMPDIR", "TMP", "TEMP"))
            assert directory.is_dir() and not list(directory.iterdir())
            raise RuntimeError("failed")
    assert os.environ["TMPDIR"] == "/shared" and os.environ["TMP"] == "/other"
    assert "TEMP" not in os.environ
    with pytest.raises(ValueError):
        with fresh_tmp(directory):
            pass


def test_reject_tmpfs():
    with pytest.raises(ValueError):
        with fresh_tmp("/dev/shm/not_created_by_test"):
            pass
