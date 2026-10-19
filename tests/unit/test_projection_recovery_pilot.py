from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools import run_profile_weighting_pilot as pilot
from tests.unit.test_profile_weighting_pilot import fake_execution
from tests.unit.test_selective_profiles import real_stage


ARMS = ("legacy_off", "projection_off", "projection_builtin_reselect",
        "projection_hmmer_score_only", "projection_hmmer_reselect")


@pytest.mark.parametrize("mode", ARMS)
def test_recovery_has_exact_projection_and_backend_flags_with_unchanged_native(mode):
    source, folder = Path("/frozen"), Path("/out")
    plan = {"cpu": 16, "experiment": "projection_recovery", "case": {
        "checkpoint": "/hits", "checkpoint_manifest": {"sha256": "pin"}, "fasta_directory": "/proteomes",
        "aligner": "mafft", "tree_builder": "FastTree", "baseline_output": "/old/native"}}
    graph, native = pilot.arm_commands(source, plan, folder, mode)
    off_graph, off_native = pilot.arm_commands(source, plan, folder, "legacy_off")
    policy, backend = pilot.recovery_configuration(mode)
    expected = off_graph.copy()
    if policy != "legacy":
        expected += ["--profile-projection-policy", "residue_preserving"]
    if backend != "off":
        expected += ["--selective-profile-mode", backend]
    assert graph == expected and native == off_native
    assert pilot.experiment_settings("projection_recovery")[0] == ARMS


def fake_recovery(tmp_path, monkeypatch, *, failed_stage=None, changed_digest=False, changed_graph=False):
    path, executions = fake_execution(tmp_path, monkeypatch, failed_stage=failed_stage)
    plan = json.loads(path.read_text())
    plan.update(experiment="projection_recovery", arms=list(ARMS), hmmer_tools={})
    plan["case"]["baseline_output"] = "/baseline"
    path.write_text(json.dumps(plan))

    def activation(graph, mode, experiment, **kwargs):
        assert experiment == "projection_recovery"
        return {"candidate_digest": "changed" if changed_digest and mode == ARMS[3] else "same",
                "completely_measured_queries": 2, "applicability": "completely_measured"}

    monkeypatch.setattr(pilot, "validate_activation", activation)
    monkeypatch.setattr(pilot, "native_pair_turnover", lambda *a: {"added_pairs": 0, "lost_pairs": 0})
    if changed_graph:
        original = pilot.timed_stage

        def timed(command, folder, stage, source, environment):
            result = original(command, folder, stage, source, environment)
            if command[0] == ARMS[3] and stage == "graph":
                graph = json.loads((folder / "graph_metrics.json").read_text())
                graph["stages"] = [{"label": "changed", "output": {"bytes": 1, "sha256": "changed"}}]
                (folder / "graph_metrics.json").write_text(json.dumps(graph))
            return result

        monkeypatch.setattr(pilot, "timed_stage", timed)
    return path, executions


def test_all_five_arms_run_with_distinct_repair_and_decision_contrasts(tmp_path, monkeypatch):
    path, executions = fake_recovery(tmp_path, monkeypatch)
    result = pilot.run(path)
    assert result["status"] == "finished_all_arms_verified"
    assert executions == [(mode, stage) for mode in ARMS for stage in ("graph", "native")]
    legacy, projection, builtin, measurement, hmmer = result["arms"]
    assert "advancement" not in legacy
    assert "advancement" in projection
    assert measurement["measurement_control_matches"] and "advancement" not in measurement
    for row in (builtin, hmmer):
        assert "advancement" in row and "matched_projection_advancement" in row


@pytest.mark.parametrize("wrong_measurement_score", [False, True])
def test_measurement_compares_to_repaired_not_legacy_control(tmp_path, monkeypatch, wrong_measurement_score):
    path, executions = fake_recovery(tmp_path, monkeypatch)
    original_timed = pilot.timed_stage

    def timed(command, folder, stage, source, environment):
        result = original_timed(command, folder, stage, source, environment)
        if stage == "graph":
            graph = json.loads((folder / "graph_metrics.json").read_text())
            graph["stages"] = [{"label": "profiles", "output": {"path": str(folder / "groups"),
                "bytes": 1, "sha256": "legacy" if command[0] == ARMS[0] else "repaired"}}]
            (folder / "graph_metrics.json").write_text(json.dumps(graph))
        return result

    def score(mode):
        return .7 if mode == ARMS[0] or (wrong_measurement_score and mode == ARMS[3]) else .71

    monkeypatch.setattr(pilot, "timed_stage", timed)
    monkeypatch.setattr(pilot, "read_partition", lambda *a: {})
    monkeypatch.setattr(pilot, "score_orthobench", lambda repo, prediction: {
        "scores": {"OrthoBench F1": score(prediction.parents[2].name)}})
    monkeypatch.setattr(pilot, "readback_scores", lambda repo, dataset, row, *a: (
        {"OrthoBench F1": score(row["label"])}, {"precision": .8, "recall": .6}, []))
    references = []

    def candidate_compare(candidate, constraints, original_candidate, original_constraints):
        references.append((candidate, original_candidate))
        return {"match": True}

    monkeypatch.setattr(pilot, "candidate_control_matches", candidate_compare)
    result = pilot.run(path)
    assert "projection_off" in references[1][1]["path"]
    if wrong_measurement_score:
        assert result["arms"][3]["status"] == "failed"
        assert result["arms"][4]["status"] == "not_run_failed_measurement_control_gate"
    else:
        assert result["status"] == "finished_all_arms_verified"
        assert result["arms"][3]["measurement_control_matches"]
        assert result["arms"][2]["advancement"]["eligible"]
        assert not result["arms"][2]["matched_projection_advancement"]["eligible"]


@pytest.mark.parametrize("stage", ["graph", "native"])
def test_failed_projection_control_blocks_all_repaired_decision_arms(tmp_path, monkeypatch, stage):
    path, executions = fake_recovery(tmp_path, monkeypatch, failed_stage=(ARMS[1], stage))
    result = pilot.run(path)
    assert result["status"] == "finished_with_explicit_gaps"
    assert result["arms"][0]["status"] == "independent_content_readback_verified"
    assert result["arms"][1]["status"] == "failed"
    assert all(r["status"] == "not_run_failed_projection_control_gate" for r in result["arms"][2:])
    assert not any(mode in ARMS[2:] for mode, _ in executions)


@pytest.mark.parametrize("kind", ["graph", "native", "digest", "graph_products"])
def test_failed_measurement_control_blocks_hmmer_without_hiding_failure(tmp_path, monkeypatch, kind):
    failed = (ARMS[3], kind) if kind in ("graph", "native") else None
    path, executions = fake_recovery(tmp_path, monkeypatch, failed_stage=failed,
        changed_digest=kind == "digest", changed_graph=kind == "graph_products")
    result = pilot.run(path)
    assert result["status"] == "finished_with_explicit_gaps"
    assert all(r["status"] == "independent_content_readback_verified" for r in result["arms"][:3])
    assert result["arms"][3]["status"] == "failed"
    assert result["arms"][4]["status"] == "not_run_failed_measurement_control_gate"
    assert not any(mode == ARMS[4] for mode, _ in executions)


def test_recovery_activation_maps_roles_and_rejects_hidden_projection_changes(monkeypatch):
    from tests.unit.test_profile_weighting_pilot import graph_fixture
    # Retain the ordinary off validation while checking policy activation.
    graph = graph_fixture("none")
    graph["parameters"]["profile_projection_policy"] = "residue_preserving"
    actual = pilot.validate_activation(graph, "projection_off", "projection_recovery")
    assert actual["projection_policy"] == "residue_preserving" and actual["selective_mode"] == "off"
    for experiment, mode in [("weighting", "none"), ("selective_hmmer", "off"),
                             ("projection_recovery", "legacy_off")]:
        with pytest.raises(ValueError, match="differs from frozen settings"):
            pilot.validate_activation(deepcopy(graph), mode, experiment)


@pytest.mark.parametrize("mode", ["builtin_reselect", "hmmer_score_only", "hmmer_reselect"])
def test_actual_backend_receipts_bind_the_profile_projection_policy(real_stage, tmp_path, mode):
    from benchmark_tools.validate_selective_profile_evidence import validate
    from orthohmm.search.matrices import get_matrix, get_background_freqs
    from orthohmm.search.msa_profile import build_msa_profile
    from orthohmm.search.profile_expansion import _decode_sequence, compute_profile_self_thresholds
    from orthohmm.search.selective_profiles import selective_profile_edges
    args, database = real_stage
    index = {gene: i for i, gene in enumerate(args[0])}
    for family, profile in args[2].items():
        args[2][family] = build_msa_profile(
            [_decode_sequence(database.get_sequence(index[g])) for g in profile.member_ids],
            profile.member_ids, get_matrix("BLOSUM62"), get_background_freqs("BLOSUM62"),
            projection_policy="residue_preserving")
    args[4] = compute_profile_self_thresholds(args[2], args[0], database, cpu=1)
    graphdir = tmp_path / "graph"
    graphdir.mkdir()
    (graphdir / "orthogroups_multipass.txt").write_text(
        "".join(" ".join(args[0][g] for g in group) + "\n" for group in args[1]))
    edges, evidence = selective_profile_edges(*args, database, mode=mode,
        output_directory=graphdir / "selective_profiles", cpu=1)
    graph = {"parameters": {"selective_profile_mode": mode, "profile_projection_policy": "residue_preserving"},
             "counts": {"strict_profile_edges": len(edges)}, "selective_profile_evidence": [evidence]}
    actual = validate(graph, mode, tmp_path, args[0])
    assert actual["completely_measured_queries"] == 1 and actual["selected_profiles"] == 2
    graph["parameters"]["profile_projection_policy"] = "legacy"
    with pytest.raises(ValueError, match="scope differs"):
        validate(graph, mode, tmp_path, args[0])
