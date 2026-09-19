from copy import deepcopy
import json
from pathlib import Path
import sys

import pytest

from benchmark_tools import cpm_replay_context as module
from benchmark_tools.checked_replay_payload_worker import FILES, validate_payload
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.prepare_qfo_corrected_replay import command_for
from tests.unit.test_checked_replay_interceptor import payload_fixture


def baseline(root):
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    output = root / "benchmarks/results/qfo_corrected_checked_replay_v1"
    checkpoint = root / "checkpoint"
    fasta = root / "fasta"
    return {"input_fastas": [{"path": str(fasta / "a.fasta")}],
        "checkpoint_manifest": {"path": str(checkpoint / "manifest.json"), "sha256": "abc"},
        "native_command": command_for(Path(sys.executable), launcher, output, fasta, checkpoint, "abc"),
        "output_root": str(output), "cwd": str(launcher),
        "environment_overrides": {"PYTHONPATH": str(launcher), "PYTHONHASHSEED": "0",
            "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"},
        "expected_stages": ["multipass", "multipass_refined", "strict_profiles", "strict_profiles_refined"]}


@pytest.mark.parametrize("arm,resolution", [("control", .1), ("cpm_low", .08), ("cpm_high", .12)])
def test_only_prespecified_argument_changes(tmp_path, arm, resolution):
    plan = baseline(tmp_path)
    before = deepcopy(plan)
    context = module.derive(tmp_path, plan, arm)
    assert plan == before
    assert context["resolution"] == resolution
    old, new = plan["native_command"], context["native_command"]
    changed = {i for i, (a, b) in enumerate(zip(old, new)) if a != b}
    expected = {old.index(flag) + 1 for flag in ("--output-directory", "--json")}
    if arm != "control":
        expected.add(old.index("--cpm-resolution") + 1)
    assert changed == expected and len(old) == len(new)
    assert context["metadata"]["seed"] == 4 and context["metadata"]["include_isolates"] is True


@pytest.mark.parametrize("problem", ["arm", "resolution", "seed", "matrix", "profiles", "output", "cwd", "env", "stages", "fastas"])
def test_changed_baseline_rejected(tmp_path, problem):
    plan = baseline(tmp_path)
    arm = "cpm_low"
    if problem == "arm":
        arm = "norm_low"
    elif problem in ("resolution", "seed", "matrix", "profiles"):
        flag = {"resolution": "--cpm-resolution", "seed": "--leiden-seed", "matrix": "--matrix", "profiles": "--profile-iterations"}[problem]
        plan["native_command"][plan["native_command"].index(flag) + 1] = "other"
    elif problem == "output":
        plan["output_root"] = "/other"
    elif problem == "cwd":
        plan["cwd"] = "/other"
    elif problem == "env":
        plan["environment_overrides"]["OMP_NUM_THREADS"] = "32"
    elif problem == "stages":
        plan["expected_stages"].pop()
    else:
        plan["input_fastas"].append({"path": "/foreign/a.fasta"})
    with pytest.raises(ValueError):
        module.derive(tmp_path, plan, arm)


@pytest.mark.parametrize("arm", ["cpm_low", "cpm_high"])
def test_nondefault_payload_requires_explicit_context(tmp_path, arm):
    payload, _, manifest = payload_fixture(tmp_path)
    context = module.derive(tmp_path, baseline(tmp_path), arm)
    manifest["output_directory"] = context["metadata"]["output_directory"]
    (payload / "metadata.json").write_text(json.dumps(context["metadata"]))
    manifest["inputs"] = [record(payload / name) for name in FILES]
    with pytest.raises(ValueError, match="settings"):
        validate_payload(manifest, payload, manifest["inputs"][0])
    assert validate_payload(manifest, payload, manifest["inputs"][0], cpm_context=context) == context["metadata"]
    wrong = deepcopy(context)
    wrong["metadata"]["seed"] = 5
    with pytest.raises(ValueError, match="settings"):
        validate_payload(manifest, payload, manifest["inputs"][0], cpm_context=wrong)


def test_changed_context_plan_rejected_before_reading_neighborhood(tmp_path):
    with pytest.raises(ValueError, match="pinned"):
        module.evidence(tmp_path, baseline(tmp_path), {"path": "/other", "sha256": "wrong"}, "cpm_low")


def test_cpm_validator_rejects_missing_corrected_provenance_before_file_access(tmp_path):
    from benchmark_tools.validate_checked_replay_payload import validate
    with pytest.raises(ValueError, match="corrected HMM"):
        validate(tmp_path / "missing", {}, tmp_path, tmp_path, cpm_arm="cpm_low")


def test_frozen_context_binding_and_input_tamper(tmp_path, monkeypatch):
    plan = baseline(tmp_path)
    baseline_path = tmp_path / "benchmarks/work/qfo_corrected_replay_commands_20260918.json"
    baseline_path.parent.mkdir(parents=True)
    baseline_path.write_text(json.dumps(plan))
    baseline_record = record(baseline_path)
    monkeypatch.setattr(module, "REPLAY_SHA", baseline_record["sha256"])
    results = tmp_path / "benchmark_tools/results"
    results.mkdir(parents=True)
    protocol = results / "QFO_PARAMETER_NEIGHBORHOOD_PROTOCOL_20260919.md"
    protocol.write_text("frozen protocol")
    monkeypatch.setattr(module, "PROTOCOL_SHA", record(protocol)["sha256"])
    input_path = tmp_path / "baseline_input.json"
    input_path.write_text("{}")
    neighborhood = {"arms": [{"label": "control", "delta": {}},
        {"label": "cpm_low", "delta": {"cpm_resolution": .08}},
        {"label": "cpm_high", "delta": {"cpm_resolution": .12}}],
        "baseline_parameters": {"cpm_resolution": .1},
        "inputs": [{"path": "baseline_input.json", "sha256": record(input_path)["sha256"]}]}
    path = results / "qfo_parameter_neighborhood_plan_20260919.json"
    path.write_text(json.dumps(neighborhood))
    monkeypatch.setattr(module, "PLAN_SHA", record(path)["sha256"])
    context = module.evidence(tmp_path, plan, baseline_record, "cpm_low")
    assert context["resolution"] == .08 and len(context["checked_records"]) == 4
    input_path.write_text("changed")
    with pytest.raises(ValueError, match="input changed"):
        module.evidence(tmp_path, plan, baseline_record, "cpm_low")
