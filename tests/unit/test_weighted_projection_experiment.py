from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools import run_profile_weighting_pilot as pilot
from benchmark_tools import weighted_projection_experiment as experiment
from tests.unit.test_profile_weighting_pilot import fake_execution, graph_fixture


def graph(label):
    arm = experiment.configuration(label)
    result = graph_fixture(arm["profile_sequence_weighting"])
    result["parameters"]["profile_projection_policy"] = arm["profile_projection_policy"]
    if arm["role"] == "candidate":
        iteration = result["profile_iterations"][0]
        iteration["profile_construction_policies"] = [
            {**r, "sequence_weighting": arm["profile_sequence_weighting"],
             "projection_policy": arm["profile_projection_policy"]} for r in iteration["profile_effective_counts"]]
    result["stages"] = []
    return result


@pytest.mark.parametrize("label", experiment.LABELS)
def test_registered_structured_policies_and_activation(label):
    arm = experiment.configuration(label)
    result = pilot.validate_activation(graph(label), label, "weighted_projection")
    assert result["arm"] == label and result["role"] == arm["role"]
    assert result["mode"] == arm["profile_sequence_weighting"]
    assert result["projection_policy"] == arm["profile_projection_policy"]
    assert result["joint_profiles"] == (2 if arm["role"] == "candidate" else 0)


@pytest.mark.parametrize("kind", ["empty", "duplicate", "mode", "policy", "nan", "zero", "effective", "requested_policy"])
def test_requested_flags_do_not_replace_actual_joint_construction_evidence(kind):
    result = graph("repair_henikoff")
    records = result["profile_iterations"][0]["profile_construction_policies"]
    if kind == "empty": records.clear()
    elif kind == "duplicate": records[1]["cluster_id"] = records[0]["cluster_id"]
    elif kind == "mode": records[0]["sequence_weighting"] = "none"
    elif kind == "policy": records[0]["projection_policy"] = "legacy"
    elif kind == "requested_policy": result["parameters"]["profile_projection_policy"] = "legacy"
    else: records[0]["effective_sequence_count"] = {"nan": float("nan"), "zero": 0., "effective": 1.1}[kind]
    with pytest.raises(ValueError):
        pilot.validate_activation(result, "repair_henikoff", "weighted_projection")


@pytest.mark.parametrize("label", experiment.LABELS[:4])
def test_constituents_cannot_silently_enable_joint_construction(label):
    result = graph(label)
    result["profile_iterations"][0]["profile_construction_policies"] = [{"unregistered": True}]
    with pytest.raises(ValueError, match="Constituent control"):
        pilot.validate_activation(result, label, "weighted_projection")


@pytest.mark.parametrize("label", experiment.LABELS)
def test_commands_use_structured_flags_and_matched_cache(label):
    plan = {"experiment": "weighted_projection", "cpu": 16, "case": {
        "checkpoint": "/hits", "checkpoint_manifest": {"sha256": "pin"}, "fasta_directory": "/proteomes",
        "aligner": "mafft", "tree_builder": "FastTree", "baseline_output": "/canonical"},
        "constituent_controls": {c: {"baseline_output": "/cache/" + c} for c in experiment.LABELS[:4]}}
    command, native = pilot.arm_commands(Path("/source"), plan, Path("/out"), label)
    arm = experiment.configuration(label)
    for flag in ("profile-sequence-weighting", "profile-projection-policy"):
        assert command[command.index("--" + flag) + 1] == arm[flag.replace("-", "_")]
    assert native[native.index("--checkpoint-source") + 1] == (
        "/cache/" + label if label in experiment.LABELS[:4] else "/canonical")
    assert native[native.index("--species-tree-mode") + 1] == "infer"
    assert "--selective-profile-mode" not in command and "--profile-alignment-policy" not in command


@pytest.mark.parametrize("label", experiment.LABELS[4:])
@pytest.mark.parametrize("kind", ["passes", "ties_stronger", "loses_stronger", "precision", "recall"])
def test_candidate_must_beat_each_constituent_without_relaxing_safeguards(label, kind):
    controls = {c: {"f1": .7 + .001 * i, "precision": .8, "recall": .6}
                for i, c in enumerate(experiment.LABELS[:4])}
    candidate = {"f1": .71, "precision": .8, "recall": .6}
    if kind == "ties_stronger": candidate["f1"] = .703
    elif kind == "loses_stronger": candidate["f1"] = .702
    elif kind == "precision": candidate["precision"] -= .006
    elif kind == "recall": candidate["recall"] -= .006
    result = experiment.matched_gate(label, candidate, controls, pilot.advancement)
    assert result["eligible"] == (kind == "passes")
    assert list(result["matched_controls"]) == experiment.configuration(label)["matched_controls"]


@pytest.mark.parametrize("label", experiment.LABELS[:4])
def test_controls_never_acquire_candidate_roles(label):
    with pytest.raises(ValueError, match="cannot become"):
        experiment.matched_gate(label, {}, {}, pilot.advancement)


def factorial_execution(tmp_path, monkeypatch, failed_stage=None, failed_graph=None, wrong_score=None):
    path, executions = fake_execution(tmp_path, monkeypatch, failed_stage=failed_stage)
    plan = json.loads(path.read_text())
    plan.update(experiment="weighted_projection", arms=list(experiment.LABELS),
                arm_parameters=experiment.specification()["arms"])
    coverage = {"genes": 3, "species": 1, "groups": 1, "pairs": 1}
    plan["case"]["baseline_output"] = "/canonical"
    plan["constituent_controls"] = {label: {
        "expected": {"f1": .7, "precision": .8, "recall": .6}, "coverage": coverage,
        "graph_stages": [], "candidate": {}, "constraints": {}, "golden": {}}
        for label in experiment.LABELS[:4]}
    # Production write_json sorts mappings. Membership, not JSON object order, is meaningful.
    path.write_text(json.dumps(plan, sort_keys=True))
    monkeypatch.setattr(pilot, "native_pair_turnover", lambda *a: {"added_pairs": 0, "lost_pairs": 0})
    monkeypatch.setattr(pilot, "replay_native", lambda *a: {**coverage, "summary": {}, "checked": []})
    original = pilot.timed_stage

    def timed(command, folder, stage, source, env):
        result = original(command, folder, stage, source, env)
        if stage == "graph":
            actual = graph(command[0])
            if command[0] == failed_graph:
                actual["stages"] = [{"label": "changed", "output": {"bytes": 1, "sha256": "different"}}]
            (folder / "graph_metrics.json").write_text(json.dumps(actual))
        return result

    monkeypatch.setattr(pilot, "timed_stage", timed)
    monkeypatch.setattr(pilot, "readback_scores", lambda repo, dataset, row, *a: (
        {"OrthoBench F1": .71 if row["label"] == wrong_score else .7}, {"precision": .8, "recall": .6}, []))
    return path, executions


def test_complete_six_arm_panel_runs_in_order_without_retuning(tmp_path, monkeypatch):
    path, executions = factorial_execution(tmp_path, monkeypatch)
    result = pilot.run(path)
    assert result["status"] == "finished_all_arms_verified"
    assert executions == [(label, stage) for label in experiment.LABELS for stage in ("graph", "native")]
    for row in result["arms"][:4]:
        assert "advancement" not in row and row["graph_control_matches"]
    assert all(not row["advancement"]["eligible"] for row in result["arms"][4:])


@pytest.mark.parametrize("control", experiment.LABELS[:4])
@pytest.mark.parametrize("stage", ["graph", "native", "graph_products", "score"])
def test_any_failed_control_prevents_both_candidates(tmp_path, monkeypatch, control, stage):
    path, executed = factorial_execution(tmp_path, monkeypatch,
        failed_stage=(control, stage) if stage in ("graph", "native") else None,
        failed_graph=control if stage == "graph_products" else None,
        wrong_score=control if stage == "score" else None)
    result = pilot.run(path)
    assert result["status"] == "finished_with_explicit_gaps"
    assert next(row for row in result["arms"] if row["label"] == control)["status"] == "failed"
    assert all(row["status"] == "not_run_failed_constituent_control_gate" for row in result["arms"][4:])
    assert not any(label in experiment.LABELS[4:] for label, _ in executed)


def test_one_failed_candidate_does_not_hide_the_other_registered_combination(tmp_path, monkeypatch):
    path, executed = factorial_execution(tmp_path, monkeypatch, failed_stage=(experiment.LABELS[4], "native"))
    result = pilot.run(path)
    assert result["status"] == "finished_with_explicit_gaps"
    assert result["arms"][4]["status"] == "failed"
    assert result["arms"][5]["status"] == "independent_content_readback_verified"
    assert (experiment.LABELS[5], "native") in executed


def test_changed_structured_configuration_fails_before_any_work(tmp_path, monkeypatch):
    path, executed = factorial_execution(tmp_path, monkeypatch)
    plan = json.loads(path.read_text())
    plan["arm_parameters"][-1]["profile_sequence_weighting"] = "none"
    path.write_text(json.dumps(plan))
    with pytest.raises(ValueError, match="structured joint arms"):
        pilot.run(path)
    assert not executed and not (tmp_path / "results.json").exists()


def test_restarting_a_completed_panel_preserves_its_first_attempt(tmp_path, monkeypatch):
    path, executed = factorial_execution(tmp_path, monkeypatch)
    pilot.run(path)
    before = (tmp_path / "results.json").read_bytes()
    with pytest.raises(ValueError, match="existing attempt"):
        pilot.run(path)
    assert (tmp_path / "results.json").read_bytes() == before


@pytest.mark.parametrize("kind", [None, "admission", "status", "missing", "point", "graph", "tree"])
def test_original_controls_require_admitted_cells_not_only_matching_score_points(tmp_path, monkeypatch, kind):
    from benchmark_tools.metric_frontier import ENDPOINTS
    from benchmark_tools.prepare_ob_candidate_neighborhood import record
    source = Path(experiment.__file__).resolve().parents[1]
    spec = experiment.specification()
    methods = []
    for arm in spec["arms"][:4]:
        label = arm["label"]
        native = tmp_path / label / "native/orthohmm_phylogeny"
        native.mkdir(parents=True)
        graph_path = native.parents[1] / "graph_metrics.json"
        graph_path.write_text(json.dumps(graph(label)))
        point = spec["control_points"][label]
        tree = {"path": str(native / "species_tree.rooted.nwk"), "bytes": 1, "sha256": "tree"}
        cell = {"dataset": "OrthoBench", "status": "independent_content_readback_verified",
            "scores": {ENDPOINTS[0]: point["f1"]}, "details": {k: point[k] for k in ("precision", "recall")},
            "checked": [record(graph_path), tree], "coverage": {},
            "candidate": {"candidate_partition": {}, "membership_constraints": {}},
            "products": {"orthohmm_root_hogs.tsv": {"path": str(native / "orthohmm_root_hogs.tsv")}}}
        methods.append({"key": arm["reference_key"], "cells": [cell]})
    if kind == "missing": methods[1]["cells"] = []
    elif kind == "point": methods[1]["cells"][0]["details"]["recall"] += .001
    elif kind == "graph": methods[1]["cells"][0]["checked"].pop(0)
    elif kind == "tree": methods[1]["cells"][0]["checked"].pop()
    compiler = tmp_path / "compiler.json"
    compiler.write_text(json.dumps({"status": "wrong" if kind == "status" else
        "admitted_individual_hmm_methods_and_projection_transfers_collated", "accuracy_result_admitted": True,
        "defaults_promoted": False, "goal_completed": False, "methods": methods}))
    archive_path = tmp_path / spec["individual_admission"]["path"]
    archive_path.parent.mkdir(parents=True)
    archive_path.write_text(json.dumps({"sources": {"compiler": {"reader": record(compiler)}}}))
    def check(pin):
        if kind == "admission" and pin["path"] == str(archive_path):
            raise ValueError("Changed immutable archive")
    monkeypatch.setattr(experiment, "check", check)
    original = (source / "benchmark_tools" / experiment.PROTOCOL).read_bytes()
    monkeypatch.setattr(experiment.subprocess, "check_output", lambda *a, **kw: original)
    if kind:
        with pytest.raises(ValueError): experiment.control_references(tmp_path, source)
    else:
        controls, pins = experiment.control_references(tmp_path, source)
        assert list(controls) == list(experiment.LABELS[:4]) and pins
        assert all(c["expected"] == spec["control_points"][label] for label, c in controls.items())
