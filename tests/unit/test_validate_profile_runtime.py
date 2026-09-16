import json
from types import SimpleNamespace

import pytest

from benchmark_tools import validate_profile_runtime as runtime


def mock_probe(monkeypatch, payload, code=0):
    def run(command, **kwargs):
        assert command[1] == "-c"
        assert kwargs["env"]["PYTHONPATH"] == str(kwargs["cwd"])
        assert kwargs["timeout"] == 60
        return SimpleNamespace(stdout=json.dumps(payload), stderr="", returncode=code)
    monkeypatch.setattr(runtime.subprocess, "run", run)


def test_success_records_runtime(monkeypatch, tmp_path):
    mock_probe(monkeypatch, {"status": "passed", "profile_length": 20})
    assert runtime.require_profile_runtime(tmp_path)["profile_length"] == 20


def test_missing_library_rejected(monkeypatch, tmp_path):
    mock_probe(monkeypatch, {"status": "failed", "error_type": "OSError"}, 1)
    with pytest.raises(ValueError, match="Profile runtime unavailable"):
        runtime.require_profile_runtime(tmp_path)


def test_nonzero_success_rejected(monkeypatch, tmp_path):
    mock_probe(monkeypatch, {"status": "passed"}, 1)
    with pytest.raises(ValueError, match="nonzero exit"):
        runtime.probe_profile_runtime(tmp_path)


def test_unknown_status_rejected(monkeypatch, tmp_path):
    mock_probe(monkeypatch, {"status": "unknown"})
    with pytest.raises(ValueError, match="Invalid profile runtime probe"):
        runtime.probe_profile_runtime(tmp_path)


def test_real_missing_checkout_rejected(tmp_path):
    # An isolated package prevents falling through to an installed OrthoHMM.
    package = tmp_path / "orthohmm"
    package.mkdir()
    (package / "__init__.py").write_text("")
    report = runtime.probe_profile_runtime(tmp_path)
    assert report["status"] == "failed"
    # Editable-install import hooks can expose the real package despite the
    # stub; the exact-source check must reject that fallback as well.
    assert report["error_type"] in {"ModuleNotFoundError", "ValueError"}
    if report["error_type"] == "ValueError":
        assert "wrong checkout" in report["error"]
    assert report["exit_code"] == 1
