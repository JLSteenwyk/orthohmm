import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from benchmark_tools import audit_scaling_root_context_measurement as module
from tests.unit.test_replay_scaling_root_context import archive, evidence, write
from tests.unit.test_verify_scaling_task_records import fixture, save


def composed(tmp_path, archive, index=0, code=0):
    metadata = tmp_path / "metadata"
    metadata.mkdir()
    args = fixture(metadata, index, code)
    directory = args[2]
    measured = archive[1]
    verified_path = directory / "verification.json"
    verified = json.loads(verified_path.read_text())
    measured["launched"] = verified["measurement"]["launched"]
    measured["native"]["exit_code"] = code
    measured["status"] = "command_exited_zero" if code == 0 else "command_failed"
    write((directory / "measurement", measured))
    prepared = json.loads((directory / "preparation.json").read_text())
    save(directory / "measurement/command.json", dict(command=prepared["measured_argv"],
        cpus=20, timeout_s=85800, interval_s=1.))
    verified.update(status=measured["status"], measurement=measured)
    save(verified_path, verified)
    args[-1] = 21816
    return args


@pytest.mark.parametrize("index", range(27))
def test_frozen_task_composes_with_raw_control_observations(tmp_path, archive, index):
    result = module.audit(*composed(tmp_path, archive, index))
    assert result["status"] == "scaling_task_measurement_audited"
    assert result["disposition"]["status"] == "native_exited_zero"
    assert result["replay"]["screening"] == archive[1]["screening"]
    for key in ("scientific_timings_admitted", "environmental_validity_established",
                "native_outputs_validated", "next_submission_authorized", "publication_ready"):
        assert result[key] is False


@pytest.mark.parametrize("code", [7, -9, 124])
def test_native_failure_is_audited_not_promoted_and_relocates(tmp_path, archive, code):
    args = composed(tmp_path, archive, 4, code)
    result = module.audit(*args)
    assert result["disposition"]["status"] == "native_exited_nonzero"
    assert result["replay"]["native_exit_code"] == code
    relocated = tmp_path / "relocated"
    shutil.copytree(args[2], relocated)
    args[2] = relocated
    other = module.audit(*args)
    assert other["disposition"] == result["disposition"]
    assert other["replay"]["screening"] == result["replay"]["screening"]


@pytest.mark.parametrize("fault", ["raw_command", "raw_point", "prepared_command", "wrapper", "job"])
def test_failure_in_either_layer_rejects_composed_audit(tmp_path, archive, fault):
    args = composed(tmp_path, archive)
    directory = args[2]
    if fault == "job":
        args[-1] += 1
    elif fault == "raw_point":
        save(directory / "measurement/point_000000.json", {})
    else:
        name = {"raw_command": "measurement/command.json", "prepared_command": "preparation.json",
                "wrapper": "verification.json"}[fault]
        path = directory / name
        data = json.loads(path.read_text())
        if fault == "raw_command": data["command"] = ["/usr/bin/true"]
        elif fault == "prepared_command": data["measured_argv"] = ["/usr/bin/true"]
        else: data["measurement"]["native_wall_s"] += 1
        save(path, data)
    with pytest.raises((ValueError, KeyError)):
        module.audit(*args)


def test_binding_evidence_changed_during_replay_rejected(tmp_path, archive, monkeypatch):
    args = composed(tmp_path, archive)
    original = module.replay
    def changed(*values):
        result = original(*values)
        save(args[2] / "preparation.json", {})
        return result
    monkeypatch.setattr(module, "replay", changed)
    with pytest.raises(ValueError):
        module.audit(*args)


def test_read_only_cli_writes_report_and_refuses_replacement(tmp_path, archive):
    plan, index, directory, recipe, sha, job = composed(tmp_path, archive, 4, 7)
    output = tmp_path / "audit.json"
    command = [sys.executable, "-m", module.__name__, "--plan", str(plan), "--index", str(index),
        "--directory", str(directory), "--recipe", str(recipe), "--recipe-sha", sha,
        "--job", str(job), "--output", str(output)]
    root = Path(__file__).resolve().parents[2]
    result = subprocess.run(command, cwd=root, capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr
    saved = output.read_bytes()
    assert json.loads(saved)["disposition"]["status"] == "native_exited_nonzero"
    repeated = subprocess.run(command, cwd=root, capture_output=True, text=True, timeout=30)
    assert repeated.returncode != 0
    assert "FileExistsError" in repeated.stderr
    assert output.read_bytes() == saved
