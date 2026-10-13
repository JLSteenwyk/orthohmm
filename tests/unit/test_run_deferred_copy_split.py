import json
from pathlib import Path

import numpy as np
import pytest

from benchmark_tools import run_deferred_copy_split as runner
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def fixture():
    species = np.array([0] * 14 + [1 + i % 77 for i in range(136)] + list(range(78)), dtype=np.int32)
    names = [f"gene{i:04d}" for i in range(len(species))]
    clusters = [list(range(150))] + [[i] for i in range(150, len(species))]
    edges = (np.array([0], dtype=np.int32), np.array([14], dtype=np.int32), np.array([2.]))
    return clusters, edges, species, names


def test_native_regeneration_and_label_blind_change_accounting(tmp_path):
    clusters, edges, species, names = fixture()
    control = runner.regenerate_refinement(clusters, edges, species, names, "default", tmp_path / "control")
    variant = runner.regenerate_refinement(clusters, edges, species, names, "phylogeny_deferred", tmp_path / "variant")
    assert clusters[0] not in control and clusters[0] in variant
    original = ([list(g) for g in clusters], [a.copy() for a in edges], species.copy(), list(names))
    repeated = runner.regenerate_refinement(clusters, edges, species, names, "default", tmp_path / "repeat")
    assert control == repeated
    assert record(tmp_path / "control")["sha256"] == record(tmp_path / "repeat")["sha256"]
    assert clusters == original[0] and all(np.array_equal(a, b) for a, b in zip(edges, original[1]))
    assert np.array_equal(species, original[2]) and names == original[3]
    result = runner.refinement_changes(clusters, control, variant, species, 500)
    assert result["retained_pre_refinement_families"] == 1
    assert result["genes_in_retained_families"] == 150
    assert result["retained_family_size_counts"] == {"150": 1}


@pytest.mark.parametrize("before,control,variant", [
    ([[0], [1]], [[0], [1]], [[0]]),
    ([[0], [1]], [[0], [1]], [[0], [0], [1]]),
    ([[0], [1]], [[0]], [[0], [1]]),
    ([[0]], [[0], [1]], [[0], [1]])])
def test_change_accounting_rejects_invalid_gene_universes(before, control, variant):
    with pytest.raises(ValueError):
        runner.refinement_changes(before, control, variant, np.array([0, 1]), 500)


def test_change_accounting_rejects_out_of_bound_restoration(tmp_path):
    clusters, edges, species, names = fixture()
    control = runner.regenerate_refinement(clusters, edges, species, names, "default", tmp_path / "control")
    with pytest.raises(ValueError, match="safeguard"):
        runner.refinement_changes(clusters, control, clusters, species, 149)


def test_failed_native_control_never_runs_variant_and_preserves_attempt(tmp_path, monkeypatch):
    source = Path(runner.__file__).resolve().parent.parent
    output = tmp_path / "output"
    output.mkdir()
    plan = output / "execution_plan.json"
    runtime = {"files": []}
    plan.write_text(json.dumps({"source_root": str(source), "repo": str(source), "output": str(output),
        "inputs": [], "sources": [], "runtime": runtime, "source_revision": "test", "scope": "test",
        "tool_path": "/usr/bin", "cases": [{"label": label, "output": str(output / label), "command": [label]}
                                             for label in runner.LABELS]}))
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    monkeypatch.setattr(runner, "runtime_records", lambda: runtime)
    monkeypatch.setattr(runner, "host_activity", lambda: {"scope": "test"})
    def construct(plan, report, destination):
        report["arms"] = [{"label": label, "status": "candidates_verified_unscored"} for label in runner.LABELS]
    monkeypatch.setattr(runner, "construct", construct)
    calls = []
    def failed(command, *args):
        calls.append(command)
        return {"returncode": 1}
    monkeypatch.setattr(runner, "timed_stage", failed)
    result = runner.run(plan)
    assert calls == [["control"]]
    assert result["status"] == "finished_with_explicit_gaps"
    assert result["arms"][0]["status"] == "failed"
    assert result["arms"][1]["status"] == "not_run_control_gate_failed"
    assert result["defaults_promoted"] is False
    with pytest.raises(ValueError, match="existing attempt"):
        runner.run(plan)


def test_candidate_preparation_failure_is_retained_without_inference(tmp_path, monkeypatch):
    source = Path(runner.__file__).resolve().parent.parent
    plan = tmp_path / "plan.json"
    runtime = {"files": []}
    plan.write_text(json.dumps({"source_root": str(source), "repo": str(source), "output": str(tmp_path),
        "inputs": [], "sources": [], "runtime": runtime, "source_revision": "test", "scope": "test",
        "tool_path": "/usr/bin", "cases": [{"label": label} for label in runner.LABELS]}))
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    monkeypatch.setattr(runner, "runtime_records", lambda: runtime)
    monkeypatch.setattr(runner, "host_activity", lambda: {})
    def failed(*args):
        raise ValueError("Default control mismatch")
    monkeypatch.setattr(runner, "construct", failed)
    result = runner.run(plan)
    assert result["status"] == "failed_preparation_or_input_gate"
    assert result["error"] == "Default control mismatch"
    assert result["arms"] == []


@pytest.mark.parametrize("cpu", [0, 17, True])
def test_preparation_rejects_unsafe_resources_before_reading_inputs(tmp_path, cpu):
    with pytest.raises(ValueError, match="1-16 CPUs"):
        runner.prepare(tmp_path, tmp_path, tmp_path / "output", tmp_path / "graph", tmp_path / "comparison", cpu)
