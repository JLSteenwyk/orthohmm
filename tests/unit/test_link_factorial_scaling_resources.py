"""Resource linkage must preserve output differences and every repeat."""

from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools.link_factorial_scaling_resources import (
    COMMON, SCOPES, STAGES, compare, partition, render, summaries,
    validate_attempt, validate_metrics,
)
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def test_partition_relabels_without_changing_membership(tmp_path):
    old = tmp_path / "old.txt"
    new = tmp_path / "new.txt"
    old.write_text("a b\nc\n")
    new.write_text("OG9: c\nOG2: b a\n")
    result = compare(partition(old, "space_separated_groups"), partition(new, "named_groups"), {"a"})
    assert result["partition_equal"] is True
    assert result["genes_in_changed_groups"] == 0


@pytest.mark.parametrize("text", ["a a\nb\n", "a b\nb\n", ""])
def test_partition_rejects_duplicates_and_empty(tmp_path, text):
    path = tmp_path / "bad.txt"
    path.write_text(text)
    with pytest.raises(ValueError):
        partition(path, "space_separated_groups")


def test_count_equality_does_not_prove_partition_equality():
    old = ({"a", "b", "c"}, {frozenset({"a", "b"}), frozenset({"c"})})
    new = ({"a", "b", "c"}, {frozenset({"a"}), frozenset({"b", "c"})})
    result = compare(old, new, {"a"})
    assert result["groups"] == 2
    assert not result["partition_equal"]
    assert not result["reference_touching_groups_equal"]
    assert result["reference_genes_in_changed_groups"] == ["a"]


def test_nonreference_changes_do_not_affect_reference_touching_groups():
    old = ({"a", "b", "c", "d"}, {frozenset({"a"}), frozenset({"b", "c"}), frozenset({"d"})})
    new = ({"a", "b", "c", "d"}, {frozenset({"a"}), frozenset({"b"}), frozenset({"c", "d"})})
    result = compare(old, new, {"a"})
    assert not result["partition_equal"]
    assert result["reference_touching_groups_equal"]
    assert result["genes_in_changed_groups"] == 3


def test_universe_mismatch_rejected():
    with pytest.raises(ValueError, match="universe"):
        compare(({"a"}, {frozenset({"a"})}), ({"b"}, {frozenset({"b"})}), set())


def native_fixture():
    argv = ["python", "-m", "orthohmm", "input", "--refinement_profile", "default", "--stop", "infer"]
    return ({"status": "complete", "metadata": dict(COMMON),
             "command": ["python", "/frozen/orthohmm/__main__.py", *argv[3:]], "cwd": "/frozen",
             "counts": {"genes": 251378, "species": 12},
             "stages": {name: {"wall_s": 1.0} for name in STAGES}},
            {"native_argv": argv, "cwd": "/frozen"})


def test_complete_native_settings():
    metrics, run = native_fixture()
    validate_metrics(metrics, run, None)


@pytest.mark.parametrize("change", ["matrix", "stage", "command", "refinement", "size", "nonfinite"])
def test_wrong_native_settings_rejected(change):
    metrics, run = native_fixture()
    if change == "matrix":
        metrics["metadata"]["substitution_matrix"] = "BLOSUM45"
    elif change == "stage":
        del metrics["stages"]["search"]
    elif change == "command":
        metrics["command"].append("--resume")
    elif change == "refinement":
        run["native_argv"][5] = "other"
        metrics["command"][4] = "other"
    elif change == "size":
        metrics["counts"]["genes"] -= 1
    else:
        metrics["stages"]["search"]["wall_s"] = float("nan")
    with pytest.raises(ValueError):
        validate_metrics(metrics, run, None)


@pytest.mark.parametrize("change", [None, "candidate_parameter", "pair_rule", "checkpoint", "species_checkpoint"])
def test_reconciliation_settings_and_checkpoint_checks(change):
    metrics, run = native_fixture()
    expansion = {"parameters": {"min_margin": 1.5}, "profile": "satellite_v2",
                 "membership_policy": "high_confidence_pair", "candidate_families": 54445,
                 "seed_families": 62885, "merges": 8440, "iterations": 2}
    metrics["metadata"].update({"phylogeny": "reconcile", "species_tree_mode": "infer",
        "species_tree_rooting": "min_variance", "phylogeny_candidates": "satellite_v2",
        "phylogeny_root_rule": "species_overlap", "phylogeny_pair_rule": "positive_paralogy",
        "phylogeny_candidate_profile": deepcopy(expansion)})
    metrics["stages"].update({name: {"wall_s": 1.0} for name in ("phylogeny", "phylogeny_candidates")})
    metrics["counts"].update({"phylogeny_checkpoint_hits": 0, "phylogeny_remapped_checkpoint_hits": 0,
                              "phylogeny_species_tree_checkpoint_hit": False})
    if change == "candidate_parameter":
        metrics["metadata"]["phylogeny_candidate_profile"]["parameters"]["min_margin"] = 1.0
    elif change == "pair_rule":
        metrics["metadata"]["phylogeny_pair_rule"] = "other"
    elif change == "checkpoint":
        metrics["counts"]["phylogeny_checkpoint_hits"] = 1
    elif change == "species_checkpoint":
        metrics["counts"]["phylogeny_species_tree_checkpoint_hit"] = True
    if change:
        with pytest.raises(ValueError):
            validate_metrics(metrics, run, expansion)
    else:
        validate_metrics(metrics, run, expansion)


def attempt_fixture():
    row = {"index": 7, "job_id": 22403, "method": "hs", "proteomes": 12, "repeat": 0,
           "comparative_timing_eligible": True, "status": "reviewed_shared_observation",
           "preflight_foreign_average_cores": 58.8, "whole_run_maximum_foreign_average_cores": 60.3,
           "resources": {"wall_seconds": 3000.0, "cpu_seconds": 86000.0, "peak_memory_bytes": 11000000000}}
    summary = {**deepcopy(row), "scheduler_state": "COMPLETED", "scheduler_exit_code": "0:0",
               "primary_resources_replayed": True, "shared_host_resources_reviewed": True,
               "execution_scope": "shared_host_matched_resources", "resource_scopes": dict(SCOPES)}
    return row, summary


def test_contention_is_accepted():
    row, summary = attempt_fixture()
    validate_attempt(row, summary, SCOPES)


@pytest.mark.parametrize("change", ["failed", "value", "scope", "unreviewed", "negative", "boolean"])
def test_invalid_resource_linkage_rejected(change):
    row, summary = attempt_fixture()
    if change == "failed":
        summary["scheduler_state"] = "FAILED"
    elif change == "value":
        summary["resources"]["wall_seconds"] += 1
    elif change == "scope":
        summary["resource_scopes"]["peak_memory_bytes"] = "sampled_rss"
    elif change == "unreviewed":
        row["comparative_timing_eligible"] = False
    else:
        value = -1 if change == "negative" else True
        row["resources"]["wall_seconds"] = summary["resources"]["wall_seconds"] = value
    with pytest.raises(ValueError):
        validate_attempt(row, summary, SCOPES)


def points_fixture():
    return [{"cell": cell, "repeat": repeat,
             "resources": {key: [100, 20, 30][repeat] for key in SCOPES},
             "partition": {"partition_equal": repeat == 0, "reference_touching_groups_equal": True}}
            for cell in ("p1_c0_r0", "p1_c1_r1") for repeat in range(3)]


def test_all_repeats_used_not_matching_or_fastest_only():
    result = summaries(points_fixture())
    assert all(row["resources"]["wall_seconds"] == {"median": 30, "minimum": 20, "maximum": 100}
               for row in result)
    assert all(row["exact_partition_repeats"] == 1 for row in result)


def test_incomplete_or_duplicate_repeats_rejected():
    with pytest.raises(ValueError, match="repeats"):
        summaries(points_fixture()[:-1])
    points = points_fixture()
    points[-1]["repeat"] = 1
    with pytest.raises(ValueError, match="repeats"):
        summaries(points)


def test_actual_report_readback_when_available():
    root = Path(__file__).resolve().parents[2]
    path = root / "benchmark_tools/results/factorial_native_resource_linkage_20261004/linkage.json"
    if not path.exists():
        pytest.skip("Actual collection not yet available")
    report = json.loads(path.read_text())
    assert report["summaries"] == summaries(report["points"])
    assert len(report["unmatched_cells"]) == 14
    assert [p["partition"]["genes_in_changed_groups"] for p in report["points"]] == [0, 0, 0, 149, 172, 0]
    assert all(p["partition"]["reference_touching_groups_equal"] for p in report["points"])
    candidates = [p["candidate_partition"] for p in report["points"] if p["candidate_partition"]]
    assert [p["genes_in_changed_groups"] for p in candidates] == [70, 160, 253]
    assert [len(p["original_only_groups"]) for p in candidates] == [10, 2, 14]
    assert all(not p["reference_genes_in_changed_groups"] for p in candidates)
    assert [s["exact_partition_repeats"] for s in report["summaries"]] == [3, 1]
    assert report["native_inference_repeated"] is False
    assert report["accuracy_recomputed"] is False
    assert report["original_factorial_stage_costs_modified"] is False
    assert report["publication_ready"] is False
    assert path.with_suffix(".md").read_text() == render(report)
    for reference in [*report["inputs"], report["source"]]:
        assert record(reference["path"]) == reference
    panel = json.loads((root / "benchmark_tools/results/threadripper_shared_panel_snapshot_20261004_v27/panel.json").read_text())
    by_index = {row["index"]: row for row in panel["runs"]}
    for point in report["points"]:
        for key in ("resources", "preflight_foreign_average_cores", "whole_run_maximum_foreign_average_cores"):
            assert point[key] == by_index[point["index"]][key]
