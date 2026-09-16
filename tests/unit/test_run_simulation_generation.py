import json
import os
from pathlib import Path
import sys

import pytest

from benchmark_tools.benchmark_production import file_record
from benchmark_tools.run_simulation_generation import child_path, execute, preflight, verify_file


def test_hash_and_path_checks(tmp_path):
    path = tmp_path / "input"
    path.write_text("abc")
    record = file_record(path, tmp_path)
    verify_file(path, record)
    path.write_text("abd")
    with pytest.raises(ValueError, match="changed"):
        verify_file(path, record)
    with pytest.raises(ValueError, match="escapes"):
        child_path(tmp_path, "../elsewhere")


def test_manifest_mismatch_fails_before_external_commands(tmp_path):
    manifest = tmp_path / "manifest.json"
    manifest.write_text("{}")
    with pytest.raises(ValueError, match="Manifest checksum"):
        preflight(manifest, "0" * 64, tmp_path, "baseline_20261001")


def test_stage_failure_preserves_evidence_and_stops_following_stages(tmp_path):
    run = {"label": "fixture", "commands": [
        {"stage": "T", "argv": [sys.executable, "-c", "print('first')"]},
        {"stage": "G", "argv": [sys.executable, "-c", "raise SystemExit(7)"]},
        {"stage": "S", "argv": [sys.executable, "-c", "raise RuntimeError('must not run')"]}]}
    evidence = tmp_path / "evidence"
    with pytest.raises(RuntimeError, match="G failed: 7"):
        execute(run, os.environ.copy(), evidence, {"fixture_manifest": "hash"})
    status = json.loads((evidence / "status.json").read_text())
    assert status["status"] == "failed"
    assert [s["status"] for s in status["stages"]] == ["complete", "failed"]
    assert status["provenance"] == {"fixture_manifest": "hash"}
    assert "first" in (evidence / "T.log").read_text()
    assert (evidence / "G.time.log").exists()
    assert not (evidence / "S.log").exists()
    with pytest.raises(FileExistsError):
        execute(run, os.environ.copy(), evidence)


def test_success_records_commands_exit_status_and_no_method_inference(tmp_path):
    run = {"label": "fixture", "commands": [{"stage": "T", "argv": [sys.executable, "-c", "pass"]}]}
    result = execute(run, os.environ.copy(), tmp_path / "evidence")
    assert result["status"] == "complete"
    assert result["stages"][0]["exit_code"] == 0
    assert result["method_inference_executed"] is False
