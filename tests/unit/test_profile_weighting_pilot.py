from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from benchmark_tools import run_profile_weighting_pilot as pilot


def test_strict_f1_gate_and_precision_recall_tolerances():
    control = {"f1": .7, "precision": .8, "recall": .6}
    assert not pilot.advancement(control, dict(control))["eligible"]
    assert pilot.advancement(control, {**control, "f1": .701})["eligible"]
    assert not pilot.advancement(control, {**control, "f1": .701, "precision": .79})["eligible"]
    assert not pilot.advancement(control, {**control, "f1": .701, "recall": .59})["eligible"]
    assert pilot.advancement(control, {**control, "f1": .701, "precision": .796})["eligible"]
    assert pilot.advancement(control, {**control, "f1": .701, "precision": .795})["eligible"]
    assert pilot.advancement(control, {**control, "f1": .701, "recall": .595})["eligible"]


@pytest.mark.parametrize("value", [True, None, np.nan, np.inf, -1, 1.1, "0.7"])
def test_gate_rejects_invalid_values(value):
    control = {"f1": .7, "precision": .8, "recall": .6}
    with pytest.raises(ValueError, match="finite normalized"):
        pilot.advancement(control, {**control, "f1": value})


def graph_fixture(mode):
    return {"parameters": {
        "profile_sequence_weighting": mode, "profile_iterations": 1, "profile_expansion": True,
        "matrix": "BLOSUM62", "cpm_resolution": .1, "leiden_seed": 4,
        "jackknife_profile_thresholds": False, "jackknife_single_copy_profiles": False,
        "profile_min_species": 1, "singleton_tie_policy": "first_hit",
    }, "counts": {"profiles_built": 2, "strict_profile_edges": 3}, "profile_iterations": [{
        "profiles_built": 2, "profile_effective_counts": ([] if mode == "none" else [
            {"cluster_id": 0, "effective_sequence_count": 1.},
            {"cluster_id": 1, "effective_sequence_count": 3.},
        ]),
    }]}


@pytest.mark.parametrize("mode", pilot.ARMS)
def test_activation_checks_exact_settings_and_actual_profiles(mode):
    result = pilot.validate_activation(graph_fixture(mode), mode)
    assert result["profiles_built"] == 2
    assert result["weighted_profiles"] == (0 if mode == "none" else 2)


@pytest.mark.parametrize("key,value", [
    ("profile_sequence_weighting", "none"), ("profile_iterations", 2), ("profile_expansion", False),
    ("matrix", "BLOSUM45"), ("cpm_resolution", .08), ("leiden_seed", 5),
    ("jackknife_profile_thresholds", True), ("jackknife_single_copy_profiles", True),
    ("profile_min_species", 3), ("singleton_tie_policy", "defer_ambiguous"),
])
def test_changed_or_inactive_settings_rejected(key, value):
    graph = graph_fixture("henikoff")
    graph["parameters"][key] = value
    with pytest.raises(ValueError, match="differs from frozen"):
        pilot.validate_activation(graph, "henikoff")


@pytest.mark.parametrize("kind", ["empty", "duplicate", "nan", "zero", "too_large", "none"])
def test_invalid_effective_counts_rejected(kind):
    graph = graph_fixture("henikoff")
    rows = graph["profile_iterations"][0]["profile_effective_counts"]
    if kind == "empty":
        rows.clear()
    elif kind == "duplicate":
        rows[1]["cluster_id"] = 0
    elif kind == "none":
        with pytest.raises(ValueError, match="Unweighted control"):
            pilot.validate_activation({**graph, "parameters": graph_fixture("none")["parameters"]}, "none")
        return
    else:
        rows[0]["effective_sequence_count"] = {"nan": np.nan, "zero": 0, "too_large": 201}[kind]
    with pytest.raises(ValueError, match="effective counts"):
        pilot.validate_activation(graph, "henikoff")


@pytest.mark.parametrize("mode", pilot.ARMS)
def test_only_graph_command_changes_across_arms(mode):
    source, folder = Path("/frozen"), Path("/out")
    plan = {"cpu": 16, "case": {"checkpoint": "/hits", "checkpoint_manifest": {"sha256": "pin"},
        "fasta_directory": "/proteomes", "aligner": "mafft", "tree_builder": "FastTree",
        "baseline_output": "/old/native"}}
    graph, native = pilot.arm_commands(source, plan, folder, mode)
    control_graph, control_native = pilot.arm_commands(source, plan, folder, "none")
    assert native == control_native
    assert graph == (control_graph if mode == "none" else control_graph + ["--profile-sequence-weighting", mode])
    assert "--species-tree-mode" in native and native[native.index("--species-tree-mode") + 1] == "infer"
    assert graph[graph.index("--cpu") + 1] == "16"


def test_no_weighting_policy_can_be_silently_substituted():
    with pytest.raises(ValueError, match="Unexpected weighting arm"):
        pilot.arm_commands(Path("/source"), {}, Path("/out"), "other")


def test_selected_baseline_requires_unique_verified_canonical_row(tmp_path):
    path = tmp_path / pilot.RETAINED
    path.parent.mkdir(parents=True)
    row = {"dataset": "OrthoBench", "label": "canonical", "status": "failed"}
    path.write_text(json.dumps({"arms": [row]}))
    with pytest.raises(ValueError, match="exactly one admitted"):
        pilot.selected_baseline(tmp_path)
    row["status"] = "independent_content_readback_verified"
    path.write_text(json.dumps({"arms": [row, deepcopy(row)]}))
    with pytest.raises(ValueError, match="exactly one admitted"):
        pilot.selected_baseline(tmp_path)


def fake_execution(tmp_path, monkeypatch, failed_stage=None, control_mismatch=False):
    source = Path(pilot.__file__).resolve().parents[1]
    plan = {"source_root": str(source), "repo": str(source), "output": str(tmp_path),
        "cpu": 16, "arms": list(pilot.ARMS), "sources": [], "inputs": [],
        "runtime": {"files": []}, "source_revision": "frozen", "tool_path": "/tools", "scope": "test",
        "case": {"fastas": [{"path": "/input/a.fa"}], "checkpoint": "/cache",
            "checkpoint_manifest": {"sha256": "pin"}, "candidate": {}, "constraints": {}, "golden": {},
            "expected": {"f1": .7, "precision": .8, "recall": .6}}}
    plan_path = tmp_path / "execution_plan.json"
    plan_path.write_text(json.dumps(plan))
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    monkeypatch.setattr(pilot, "UNIVERSE", (3, 1))
    monkeypatch.setattr(pilot, "runtime", lambda: {"files": []})
    monkeypatch.setattr(pilot, "check", lambda pin: None)
    monkeypatch.setattr(pilot, "record", lambda path: {"path": str(path)})
    monkeypatch.setattr(pilot, "host_activity", lambda: {})
    monkeypatch.setattr(pilot.SeqIO, "parse", lambda *args: [SimpleNamespace(id=name, seq="ACDE") for name in ("a", "b", "c")])
    monkeypatch.setattr(pilot, "audit", lambda *args: {})
    monkeypatch.setattr(pilot, "load_accuracy_checkpoint", lambda *args, **kwargs: (
        ["a", "b", "c"], np.zeros(3, dtype=np.int32), [], [], []))
    monkeypatch.setattr(pilot, "validate_ownership", lambda *args: {"genes": 3, "species": 1})
    monkeypatch.setattr(pilot, "arm_commands", lambda source, plan, folder, mode: ([mode], [mode]))
    executions = []

    def timed(command, folder, name, source, env):
        mode = command[0]
        executions.append((mode, name))
        if name == "graph":
            graph = graph_fixture(mode)
            graph["stages"] = []
            (folder / "graph_metrics.json").write_text(json.dumps(graph))
        return {"returncode": int((mode, name) == failed_stage)}

    monkeypatch.setattr(pilot, "timed_stage", timed)
    monkeypatch.setattr(pilot, "prepare_partition", lambda seed, folder, *args: {
        "candidate_partition": {"path": str(folder / "groups")}, "membership_constraints": {"path": str(folder / "constraints")}})

    def compare(*args):
        if control_mismatch:
            raise ValueError("Unchanged candidate control mismatch")
        return {"match": True}

    monkeypatch.setattr(pilot, "candidate_control_matches", compare)
    monkeypatch.setattr(pilot, "compare_golden", compare)
    monkeypatch.setattr(pilot, "validate_partition", lambda *args: {"all_candidate_genes_preserved": True})
    monkeypatch.setattr(pilot, "replay_native", lambda *args: {"genes": 3, "species": 1, "summary": {}, "checked": []})
    monkeypatch.setattr(pilot, "score_orthobench", lambda *args: {"scores": {"OrthoBench F1": .7}})
    monkeypatch.setattr(pilot, "readback_scores", lambda *args: (
        {"OrthoBench F1": .7}, {"precision": .8, "recall": .6}, []))
    return plan_path, executions


def test_all_arms_run_sequentially_and_neutral_scores_do_not_advance(tmp_path, monkeypatch):
    path, executed = fake_execution(tmp_path, monkeypatch)
    result = pilot.run(path)
    assert result["status"] == "finished_all_arms_verified"
    assert executed == [(mode, stage) for mode in pilot.ARMS for stage in ("graph", "native")]
    assert all(not row["advancement"]["eligible"] for row in result["arms"][1:])


def test_failed_control_prevents_all_weighted_runs(tmp_path, monkeypatch):
    path, executed = fake_execution(tmp_path, monkeypatch, control_mismatch=True)
    result = pilot.run(path)
    assert result["status"] == "finished_with_explicit_gaps"
    assert executed == [("none", "graph")]
    assert [row["status"] for row in result["arms"]] == ["failed"] + ["not_run_failed_control_gate"] * 3


def test_one_failed_weighted_arm_does_not_suppress_other_prespecified_arms(tmp_path, monkeypatch):
    path, executed = fake_execution(tmp_path, monkeypatch, failed_stage=("henikoff", "graph"))
    result = pilot.run(path)
    assert result["arms"][1]["status"] == "failed"
    assert all(row["status"] == "independent_content_readback_verified" for row in result["arms"][2:])
    assert ("species_henikoff_cluster", "native") in executed


def test_existing_execution_receipt_is_never_overwritten(tmp_path, monkeypatch):
    path, executed = fake_execution(tmp_path, monkeypatch)
    receipt = tmp_path / "results.json"
    receipt.write_text("preserve")
    with pytest.raises(ValueError, match="existing attempt"):
        pilot.run(path)
    assert receipt.read_text() == "preserve" and not executed


def test_preparation_requires_frozen_source_and_fresh_namespace(tmp_path):
    with pytest.raises(ValueError, match="frozen scientific checkout"):
        pilot.prepare(tmp_path, tmp_path, tmp_path / "new")
    with pytest.raises(ValueError, match="fresh output"):
        pilot.prepare(Path("/repo"), Path("/source"), tmp_path)


def alignment_graph_fixture(mode):
    graph = graph_fixture("none")
    graph["parameters"]["profile_alignment_policy"] = mode
    iteration = graph["profile_iterations"][0]
    iteration.update(alignment_profiles=[], alignment_evidence=[])
    if mode == "off":
        return graph
    iteration["alignment_profiles"] = [{"cluster_id": 0, "member_ids": ["anchor0", "anchor1", "anchor2"]}]
    for index, matched in enumerate((8, 3)):
        keep = matched >= 8
        iteration["alignment_evidence"].append({
            "gene_id": 3 + index, "gene_name": f"query{index}", "profile_cluster_id": 0,
            "anchor_gene_id": index, "anchor_name": f"anchor{index}", "anchor_score": .8,
            "model_length": 10, "target_length": 10, "score": matched * 5,
            "model_start": 0, "model_end": 10, "target_start": 0, "target_end": 10,
            "matched_known_residues": matched, "passes_guard": keep,
            "category": "both_sides_covered" if keep else "partial_shared_region",
            "aligned_model_coverage": matched / 10, "aligned_target_coverage": matched / 10,
            "bounding_model_coverage": 1., "bounding_target_coverage": 1.,
            "coordinate_system": "zero_based_half_open",
        })
    graph["counts"]["strict_profile_edges"] = 1
    return graph


def test_alignment_experiment_uses_its_own_frozen_protocol_and_two_arms():
    arms, protocol = pilot.experiment_settings("alignment")
    assert arms == ("off", "strict_span_guard")
    assert protocol == "hmm_aligned_profile_guard_20261018.json"
    with pytest.raises(ValueError, match="Unknown profile-evidence experiment"):
        pilot.experiment_settings("unknown")


@pytest.mark.parametrize("mode", ["off", "strict_span_guard"])
def test_alignment_command_changes_only_the_guard_option(mode):
    plan = {"experiment": "alignment", "cpu": 16, "case": {
        "checkpoint": "/hits", "checkpoint_manifest": {"sha256": "pin"},
        "fasta_directory": "/proteomes", "aligner": "mafft", "tree_builder": "FastTree", "baseline_output": "/old/native"}}
    graph, native = pilot.arm_commands(Path("/frozen"), plan, Path("/out"), mode)
    control_graph, control_native = pilot.arm_commands(Path("/frozen"), plan, Path("/out"), "off")
    assert native == control_native
    assert graph == (control_graph if mode == "off" else control_graph + ["--profile-alignment-policy", mode])
    assert "--profile-sequence-weighting" not in graph


@pytest.mark.parametrize("mode", ["off", "strict_span_guard"])
def test_alignment_activation_is_recounted_from_coordinates_and_unique_edges(mode):
    result = pilot.validate_activation(alignment_graph_fixture(mode), mode, "alignment")
    assert result["weighted_profiles"] == 0
    if mode != "off":
        assert result["alignment_traced"] == 2
        assert result["alignment_rejected"] == result["alignment_retained"] == 1


@pytest.mark.parametrize("field,value", [
    ("passes_guard", False), ("matched_known_residues", 3), ("aligned_target_coverage", .2),
    ("bounding_model_coverage", .5), ("category", "asymmetric_coverage"),
    ("coordinate_system", "one_based"), ("model_end", 11), ("score", 0),
    ("anchor_name", "wrong_member"), ("anchor_gene_id", 3), ("gene_id", -1),
])
def test_tampered_alignment_evidence_cannot_pass_activation(field, value):
    graph = alignment_graph_fixture("strict_span_guard")
    graph["profile_iterations"][0]["alignment_evidence"][0][field] = value
    with pytest.raises(ValueError):
        pilot.validate_activation(graph, "strict_span_guard", "alignment")


@pytest.mark.parametrize("kind", ["empty", "missing_profiles", "duplicate_gene", "edge_count", "weighted"])
def test_inactive_incomplete_or_combined_alignment_evidence_rejected(kind):
    graph = alignment_graph_fixture("strict_span_guard")
    iteration = graph["profile_iterations"][0]
    if kind == "empty":
        iteration["alignment_evidence"] = []
    elif kind == "missing_profiles":
        iteration["alignment_profiles"] = []
    elif kind == "duplicate_gene":
        iteration["alignment_evidence"][1]["gene_id"] = iteration["alignment_evidence"][0]["gene_id"]
    elif kind == "edge_count":
        graph["counts"]["strict_profile_edges"] += 1
    else:
        graph["parameters"]["profile_sequence_weighting"] = "henikoff"
    with pytest.raises(ValueError):
        pilot.validate_activation(graph, "strict_span_guard", "alignment")


def test_alignment_off_control_cannot_silently_run_an_active_guard():
    graph = alignment_graph_fixture("strict_span_guard")
    graph["parameters"]["profile_alignment_policy"] = "off"
    with pytest.raises(ValueError, match="Off control"):
        pilot.validate_activation(graph, "off", "alignment")


def test_common_runner_executes_two_arm_alignment_panel_with_original_control_gate(tmp_path, monkeypatch):
    path, executed = fake_execution(tmp_path, monkeypatch)
    plan = json.loads(path.read_text())
    plan.update(experiment="alignment", arms=["off", "strict_span_guard"])
    path.write_text(json.dumps(plan))
    original = pilot.timed_stage

    def timed(command, folder, name, source, env):
        result = original(command, folder, name, source, env)
        if name == "graph":
            graph = alignment_graph_fixture(command[0])
            graph["stages"] = []
            (folder / "graph_metrics.json").write_text(json.dumps(graph))
        return result

    monkeypatch.setattr(pilot, "timed_stage", timed)
    result = pilot.run(path)
    assert result["experiment"] == "alignment" and result["status"] == "finished_all_arms_verified"
    assert executed == [(mode, stage) for mode in plan["arms"] for stage in ("graph", "native")]
    assert result["arms"][1]["activation"]["alignment_rejected"] == 1
    assert not result["arms"][1]["advancement"]["eligible"]
