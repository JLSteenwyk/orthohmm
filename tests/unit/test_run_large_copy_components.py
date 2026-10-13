from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pytest

from benchmark_tools import run_large_copy_components as runner
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def comparison():
    methods = []
    for label in runner.BRANCH_LABELS:
        scores = dict.fromkeys(runner.ENDPOINTS, .5)
        cells = [{"dataset": panel, "status": "independent_content_readback_verified",
                  "scores": {k: v for k, v in scores.items() if k.startswith(panel + ": ")}}
                 for panel in runner.EXPECTED_GENES]
        methods.append({"key": "copy_branch:" + label, "scores": scores, "cells": cells})
    return {"methods": methods}


@pytest.mark.parametrize("dataset", list(runner.EXPECTED_GENES))
def test_all_complete_branch_cells_are_required_before_next_round(dataset):
    result = runner.admitted_branch_cells(comparison(), dataset)
    assert tuple(result) == runner.BRANCH_LABELS
    assert all(c["dataset"] == dataset for c in result.values())


@pytest.mark.parametrize("mutation", ["missing_method", "duplicate_method", "missing_endpoint", "missing_cell",
                                    "duplicate_cell", "live_cell", "failed_cell", "invalid_metric", "vector_mismatch"])
def test_partial_unverified_or_inconsistent_branch_vector_cannot_enable_inference(mutation):
    value = comparison()
    method = value["methods"][-1]
    if mutation == "missing_method":
        value["methods"].pop()
    elif mutation == "duplicate_method":
        value["methods"].append(deepcopy(method))
    elif mutation == "missing_endpoint":
        method["scores"].pop(next(iter(method["scores"])))
    elif mutation == "missing_cell":
        method["cells"].pop()
    elif mutation == "duplicate_cell":
        method["cells"].append(deepcopy(method["cells"][0]))
    elif mutation in ("live_cell", "failed_cell"):
        method["cells"][0]["status"] = "native_running" if mutation == "live_cell" else "failed"
    elif mutation == "invalid_metric":
        key = next(iter(method["cells"][0]["scores"]))
        method["cells"][0]["scores"][key] = float("nan")
    else:
        method["scores"][next(iter(method["scores"]))] = .7
    with pytest.raises(ValueError):
        runner.admitted_branch_cells(value, "QfO")


def trace():
    baseline = {"members": list(range(600)), "size": 600, "large_copy_eligible": True,
                "medium_copy_eligible": False, "deferred": False, "maximum_species_copies": 14,
                "dataset_species_count": 78}
    variant = {**baseline, "component_min_edge_weight": 1.5, "retain_large_copy_components_max_size": 500,
               "component_sizes": [501, 2, *[1] * 97], "retained_component_sizes": [2],
               "oversized_component_sizes": [501], "isolated_genes": 97}
    return [baseline], [variant]


def test_executed_trace_counts_include_oversized_fallback_and_isolates():
    assert runner.executed_component_counts(*trace(), "QfO") == {
        "families": 1, "genes": 600, "components": 99, "bounded_components": 1,
        "bounded_component_genes": 2, "oversized_components": 1, "oversized_component_genes": 501,
        "isolated_genes": 97}


@pytest.mark.parametrize("field,value", [("deferred", True), ("component_min_edge_weight", 1.),
    ("retain_large_copy_components_max_size", 501), ("component_sizes", [500, 100]),
    ("component_sizes", [True, 599]), ("retained_component_sizes", [501, 2]),
    ("oversized_component_sizes", []), ("isolated_genes", 0), ("maximum_species_copies", 20)])
def test_changed_executed_policy_is_rejected(field, value):
    control, variant = trace()
    variant[0][field] = value
    with pytest.raises(ValueError):
        runner.executed_component_counts(control, variant, "QfO")


def test_narrow_trace_must_be_identical_and_large_activation_is_rejected():
    with pytest.raises(ValueError, match="narrow panel"):
        runner.executed_component_counts(*trace(), "OrthoBench")
    event = {"members": list(range(70)), "size": 70, "large_copy_eligible": False,
             "medium_copy_eligible": True, "deferred": False}
    assert runner.executed_component_counts([event], [dict(event)], "OrthoBench") == {}
    with pytest.raises(ValueError, match="Non-large"):
        runner.executed_component_counts([event], [{**event, "component_sizes": [70]}], "OrthoBench")


def candidate(seed, target, *args):
    target.mkdir(parents=True)
    path, constraints = target / "groups.txt", target / "constraints.json"
    path.write_bytes(seed.read_bytes())
    constraints.write_text(json.dumps({"groups": seed.read_text()}))
    return {"seed_partition": record(seed), "candidate_partition": record(path),
            "membership_constraints": record(constraints), "expansion": {"parameters": {"unchanged": True}}}


@pytest.mark.parametrize("dataset", ["QfO", "OrthoBench"])
def test_construction_uses_actual_numeric_refinement_and_exact_control(tmp_path, monkeypatch, dataset):
    species_count = 78 if dataset == "QfO" else 12
    species = np.array([0] * 14 + [1 + i % (species_count - 1) for i in range(136)]
                       + list(range(species_count)), dtype=np.int32)
    names = ["g" + str(i) for i in range(len(species))]
    clusters = [list(range(150)), *[[i] for i in range(150, len(species))]]
    empty = np.array([], dtype=np.int32)
    edges = (np.array([0]), np.array([1]), np.array([2.]))
    monkeypatch.setattr(runner, "load_graph", lambda *a: (names, species, (empty, empty, []),
        clusters, edges, None, {}, {"genes": len(names), "species": species_count}))
    monkeypatch.setattr(runner, "prepare_partition", candidate)
    seed = tmp_path / "old_seed.txt"
    runner.regenerate_refinement(clusters, edges, species, names, "default", seed)
    original = {"seed": record(seed), "candidate": candidate(seed, tmp_path / "old_candidate")}
    output = tmp_path / "new"
    output.mkdir()
    expected = {"families": 1, "genes": 150, "components": 149, "bounded_components": 1,
                "bounded_component_genes": 2, "oversized_components": 0, "oversized_component_genes": 0,
                "isolated_genes": 148}
    rows, diagnostics = runner.construct({"dataset": dataset, "structural_control": original,
        "expected_component_counts": expected}, output)
    assert rows["control"]["seed"]["sha256"] == original["seed"]["sha256"]
    assert runner.native_inputs_unchanged(rows) is (dataset != "QfO")
    assert diagnostics["component_counts"] == (expected if dataset == "QfO" else {})
    assert all(Path(r["copy_split_trace"]["path"]).is_file() for r in rows.values())


def worker_fixture(tmp_path, monkeypatch, dataset="OrthoBench", changed=False):
    source = Path(runner.__file__).resolve().parent.parent
    pin = tmp_path / "pin.txt"
    pin.write_text("fixture\n")
    artifact = record(pin)
    scores = {k.split(": ", 1)[1]: .5 for k in runner.ENDPOINTS if k.startswith(dataset + ": ")}
    native = {"genes": runner.EXPECTED_GENES[dataset], "species": 78 if dataset == "QfO" else 12,
              "groups": 1, "pairs": 1, "products": {}, "checked": [], "summary": {"species_tree_families": 2}}
    output = tmp_path / "components500"
    (output / "output/orthohmm_phylogeny").mkdir(parents=True)
    for name in runner.PRODUCTS:
        p = output / "output/orthohmm_phylogeny" / name
        p.write_text("fixture\n")
    row = {"seed": artifact, "copy_split_trace": artifact,
           "candidate": {"seed_partition": artifact, "candidate_partition": artifact, "membership_constraints": artifact}}
    rows = {"control": deepcopy(row), "components500": deepcopy(row)}
    if changed:
        rows["components500"]["candidate"]["candidate_partition"]["sha256"] = "changed_fixture"
    plan = {"source_root": str(source), "repo": str(source), "output": str(tmp_path), "dataset": dataset,
            "source_revision": "fixture", "scope": "unit fixture", "tool_path": "/usr/bin", "case": {
                "output": str(output), "command": [], "scoring_id": "short_fixture"},
            "parent_native_row": {"label": "control"}, "parent_native_case": {"output": str(tmp_path / "old"), "scoring_id": "old_fixture"},
            "matched_control": {"scores": runner.endpoint_vector(dataset, scores), "products": native["products"],
                                "native_summary": native["summary"]}}
    path = tmp_path / "plan.json"
    path.write_text(json.dumps(plan))
    calls = []
    monkeypatch.setattr(runner, "verify", lambda *a: None)
    monkeypatch.setattr(runner, "host_activity", lambda: "shared-host unit fixture")
    monkeypatch.setattr(runner, "construct", lambda *a: (rows, {}))
    monkeypatch.setattr(runner, "replay_native", lambda *a: calls.append("native_readback") or native)
    monkeypatch.setattr(runner, "readback_scores", lambda *a: calls.append("score_readback") or (scores, {}, []))
    return path, calls


@pytest.mark.parametrize("dataset", list(runner.EXPECTED_GENES))
def test_exact_native_input_reuse_rechecks_native_and_scores_without_rerunning(tmp_path, monkeypatch, dataset):
    path, calls = worker_fixture(tmp_path, monkeypatch, dataset)
    monkeypatch.setattr(runner, "timed_stage", lambda *a: pytest.fail("No duplicate native inference"))
    monkeypatch.setattr(runner, "score_qfo", lambda *a: pytest.fail("No duplicate scoring"))
    result = runner.run(path)
    assert result["status"] == "independent_content_readback_verified"
    assert calls == ["native_readback", "score_readback"]
    assert "native_execution" not in result
    assert "reuse_no_native_rerun" in result["execution_kind"]
    with pytest.raises(ValueError, match="untouched"):
        runner.run(path)


def test_changed_qfo_component_inputs_run_own_tree_and_read_all_scores(tmp_path, monkeypatch):
    path, calls = worker_fixture(tmp_path, monkeypatch, "QfO", changed=True)
    monkeypatch.setattr(runner, "timed_stage", lambda *a: calls.append("native_execution") or {"returncode": 0})
    monkeypatch.setattr(runner, "validate_partition", lambda *a: {})
    monkeypatch.setattr(runner, "score_qfo", lambda *a: calls.append("new_scoring") or {})
    result = runner.run(path)
    assert result["status"] == "independent_content_readback_verified"
    assert calls == ["native_execution", "native_readback", "new_scoring", "score_readback"]
    assert len(result["scores"]) == 6 and result["execution_kind"] == "new_cached_graph_component_inference_own_tree"


def test_failed_input_gate_prevents_candidate_native_and_score_work(tmp_path, monkeypatch):
    path, calls = worker_fixture(tmp_path, monkeypatch, "QfO")
    def fail(*a):
        raise ValueError("changed pinned input")
    monkeypatch.setattr(runner, "verify", fail)
    monkeypatch.setattr(runner, "construct", lambda *a: pytest.fail("No candidates before verification"))
    result = runner.run(path)
    assert result["status"] == "failed_with_explicit_gap" and result["scores"] == {} and calls == []


def test_failed_native_attempt_is_preserved_without_scoring(tmp_path, monkeypatch):
    path, calls = worker_fixture(tmp_path, monkeypatch, "QfO", changed=True)
    monkeypatch.setattr(runner, "timed_stage", lambda *a: {"returncode": 1})
    monkeypatch.setattr(runner, "score_qfo", lambda *a: pytest.fail("No scoring after native failure"))
    result = runner.run(path)
    assert result["status"] == "failed_with_explicit_gap" and result["native_execution"]["returncode"] == 1
    assert calls == [] and result["scores"] == {}


def original_command():
    return ["python", "old/replay_phylogeny.py", "--fasta-directory", "fasta", "--candidate-clusters", "candidate",
        "--membership-constraints", "constraints", "--output-directory", "old/output", "--json", "old/metrics.json",
        "--cpu", "16", "--species-tree-mode", "infer", "--root-rule", "species_overlap",
        "--pair-rule", "positive_paralogy", "--species-tree-rooting", "min_variance",
        "--satellite-membership-rule", "high_confidence_pair", "--checkpoint-source", "old/cache",
        "--aligner", "mafft", "--tree-builder", "fasttree"]


def protocol_plan(tmp_path):
    source = Path(runner.__file__).resolve().parent.parent
    folder = tmp_path / "components500"
    candidate = folder / "candidate/orthohmm_working_res/orthohmm_edges_clustered.txt"
    return {"source_root": str(source), "output": str(tmp_path), "cpu": 16,
        "parent_native_case": {"command": original_command()}, "checkpoint_source": "old/cache",
        "case": {"label": "components500", "profile": runner.PROFILE, "output": str(folder),
            "command": runner.native_command(original_command(), source, folder, candidate,
                candidate.with_name("phylogeny_candidate_merges.json"), Path("old/cache"))},
        "inputs": [], "sources": [], "runtime": {"files": []}}


@pytest.mark.parametrize("mutation", ["none", "cpu", "profile", "label", "output", "tree_rule", "extra_flag", "runtime", "cuda"])
def test_verification_enforces_exact_native_protocol_and_runtime(tmp_path, monkeypatch, mutation):
    plan = protocol_plan(tmp_path)
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    monkeypatch.setattr(runner, "runtime", lambda: {"files": []})
    if mutation == "cpu":
        plan["cpu"] = 32
    elif mutation in ("profile", "label", "output"):
        plan["case"][mutation] = "wrong"
    elif mutation == "tree_rule":
        command = plan["case"]["command"]
        command[command.index("--root-rule") + 1] = "supported_children"
    elif mutation == "extra_flag":
        plan["case"]["command"].append("--unconstrained-membership")
    elif mutation == "runtime":
        monkeypatch.setattr(runner, "runtime", lambda: {"files": [], "changed": True})
    elif mutation == "cuda":
        monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "0")
    if mutation == "none":
        runner.verify(plan)
    else:
        with pytest.raises(ValueError):
            runner.verify(plan)


def test_preparation_stops_on_live_vectors_before_loading_parent_or_creating_output(tmp_path, monkeypatch):
    source = Path(runner.__file__).resolve().parent.parent
    value = comparison()
    value["methods"][-1]["scores"].pop(next(iter(value["methods"][-1]["scores"])))
    path = tmp_path / "comparison.json"
    path.write_text(json.dumps(value))
    monkeypatch.setattr(runner.subprocess, "run", lambda *a, **k: None)
    output = tmp_path / "new"
    with pytest.raises(ValueError, match="complete independently admitted"):
        runner.prepare(source, source, output, "QfO", tmp_path / "missing_parent", path)
    assert not output.exists()


def test_corrupt_matched_native_reuse_is_not_admitted(tmp_path, monkeypatch):
    path, calls = worker_fixture(tmp_path, monkeypatch)
    plan = json.loads(path.read_text())
    plan["matched_control"]["native_summary"]["species_tree_families"] = 3
    path.write_text(json.dumps(plan))
    result = runner.run(path)
    assert result["status"] == "failed_with_explicit_gap" and result["scores"] == {}
    assert "products differ" in result["error"]
