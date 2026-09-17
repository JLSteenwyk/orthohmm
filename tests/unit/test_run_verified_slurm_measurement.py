import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "benchmark_tools"))
from run_verified_slurm_measurement import check_manifests, run_checked
from snapshot_runtime_trees import digest, inventory


def test_manifest_reverification_detects_changes(tmp_path):
    runtime = tmp_path / "runtime"
    runtime.mkdir()
    binary = runtime / "binary"
    binary.write_text("original")
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps(inventory([runtime])))
    specs = [(manifest, digest(manifest))]
    assert check_manifests(specs)[0]["status"] == "runtime_tree_identity_matches"
    binary.write_text("changed")
    with pytest.raises(ValueError, match="changed"):
        check_manifests(specs)


def test_wrong_manifest_hash_rejected(tmp_path):
    path = tmp_path / "manifest"
    path.write_text("{}")
    with pytest.raises(ValueError, match="digest differs"):
        check_manifests([(path, "0" * 64)])
    with pytest.raises(ValueError, match="at least one"):
        check_manifests([])


@pytest.mark.parametrize("status", ["command_exited_zero", "command_failed", "command_timed_out", "measurement_failed"])
def test_checks_bracket_every_terminal_status(tmp_path, status):
    events = []
    def checker(specs):
        events.append("check")
        return []
    def measure(path):
        events.append("measure")
        assert path.name == "measurement"
        return {"status": status}
    result = run_checked([], tmp_path / "run", measure, checker)
    assert events == ["check", "measure", "check"]
    assert result["status"] == status
    assert not result["scientific_results_admitted"]


def test_failed_preflight_never_launches(tmp_path):
    def fail(specs):
        raise ValueError("wrong runtime")
    result = run_checked([], tmp_path / "run", lambda _: pytest.fail("launched"), fail)
    assert result["status"] == "verified_wrapper_failed"
    assert "after" not in result


def test_postcheck_rejects_successful_command(tmp_path):
    calls = []
    def check(specs):
        calls.append(1)
        if len(calls) == 2:
            raise ValueError("changed")
        return []
    result = run_checked([], tmp_path / "run", lambda _: {"status": "command_exited_zero"}, check)
    assert result["status"] == "runtime_changed_or_unverifiable"


def test_measurement_exception_still_checks_after(tmp_path):
    calls = []
    def check(specs):
        calls.append(1)
        return []
    def measure(path):
        raise RuntimeError("collector failed")
    result = run_checked([], tmp_path / "run", measure, check)
    assert len(calls) == 2
    assert result["status"] == "verified_wrapper_failed"
    assert json.loads((tmp_path / "run/verification.json").read_text())["error"] == "collector failed"


def test_no_overwrite(tmp_path):
    with pytest.raises(FileExistsError):
        run_checked([], tmp_path, lambda _: None)
