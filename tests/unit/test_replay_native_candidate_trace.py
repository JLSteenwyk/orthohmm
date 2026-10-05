"""Accepted-union replay is whole-partition evidence, not score replay."""

import json
from pathlib import Path

import pytest

from benchmark_tools.link_factorial_scaling_resources import compare
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.replay_native_candidate_trace import collect, render, replay
from tests.unit.test_diagnose_candidate_trace_variation import row


def seeds(*groups):
    return {frozenset(g) for g in groups}


def test_two_rounds_reconstruct_components():
    initial = seeds(["s"], ["t"], ["a", "b"], ["untouched"])
    trace = [row(), row("t", target=["s", "a", "b"], iteration=1)]
    universe, result = replay(initial, trace)
    assert universe == {"s", "t", "a", "b", "untouched"}
    assert result == seeds(["s", "t", "a", "b"], ["untouched"])


def test_multiple_attachments_use_round_start_not_partial_union():
    assert replay(seeds(["s"], ["t"], ["a", "b"]),
                  [row(), row("t", source_cluster=2)])[1] == seeds(["s", "t", "a", "b"])


def test_conflicting_round_labels_rejected():
    with pytest.raises(ValueError, match="Conflicting"):
        replay(seeds(["s"], ["t"], ["a", "b"]), [row(), row("t")])


@pytest.mark.parametrize("initial", [set(), seeds([]), seeds(["s", "a"], ["a", "b"])])
def test_invalid_partition_rejected(initial):
    with pytest.raises(ValueError, match="seed partition"):
        replay(initial, [row()])


def test_partial_or_unknown_endpoint_rejected():
    with pytest.raises(ValueError, match="round-start"):
        replay(seeds(["s", "extra"], ["a", "b"]), [row()])


def test_current_round_union_not_valid_as_round_start_endpoint():
    with pytest.raises(ValueError, match="round-start"):
        replay(seeds(["s"], ["t"], ["a", "b"]),
               [row(), row("t", target=["s", "a", "b"], source_cluster=2)])


def test_cycle_rejected_even_with_valid_endpoints():
    trace = [row("s", target=["a"], target_size=1, target_cluster=0),
             row("a", target=["t"], target_size=1, source_cluster=0, target_cluster=2),
             row("s", target=["t"], target_size=1, target_cluster=2)]
    with pytest.raises(ValueError, match="Redundant"):
        replay(seeds(["s"], ["t"], ["a"]), trace)


@pytest.fixture
def retained(tmp_path):
    def write(name, value):
        path = tmp_path / name
        path.write_text(json.dumps(value) if isinstance(value, (dict, list)) else value)
        return record(path)

    old_trace = [row(), row("x", source_cluster=2)]
    new_trace = [row(), row("y", source_cluster=3)]
    seed = write("seed.txt", "s\nx\ny\na b\n")
    old = write("original.txt", "s x a b\ny\n")
    new = write("prediction.txt", "s y a b\nx\n")
    old_trace_ref = write("original_trace.json", old_trace)
    new_trace_ref = write("phylogeny_candidate_merges.json", new_trace)
    prep = {"candidate_arms": {"p0_c1": {"seed_partition": seed, "candidate_partition": old,
            "membership_constraints": old_trace_ref,
            "expansion": {"parameters": {"max_satellites_per_anchor": 2}}}}}
    old_groups, new_groups = seeds(["s", "x", "a", "b"], ["y"]), seeds(["s", "y", "a", "b"], ["x"])
    universe = {"s", "x", "y", "a", "b"}
    delta = compare((universe, old_groups), (universe, new_groups), set())
    delta["native_only_groups"] = delta.pop("scaling_only_groups")
    score = {"schema": "native_factorial_orthobench_score_v1", "status": "terminal_native_orthobench_scored",
             "cell": "p0_c1_r0", "prediction_format": "space_separated_groups", "job_id": 1,
             "native_outputs_validated": True, "original_prediction": old, "prediction": new,
             "evidence": [new_trace_ref], "canonical_partition_comparison": delta}
    return write, prep, score


def test_collect_checks_whole_partitions_and_preserves_limitations(retained):
    write, prep, score = retained
    report = collect(write("preparation.json", prep), write("score.json", score))
    assert report["full_original_partition_reconstructed"]
    assert report["full_native_partition_reconstructed"]
    assert report["trace_comparison"]["common_semantic_merges"] == 1
    assert len(report["native_only_accepted_merges"]) == 1
    assert not report["accuracy_rescored"]
    assert not report["frozen_method_modified"]
    assert "not replayed" in report["limitations"][0]
    assert "Both complete partitions" in render(report)
    json.dumps(report, allow_nan=False)


@pytest.mark.parametrize("updates", [{"cell": "p0_c1_r1"}, {"native_outputs_validated": False},
                                     {"prediction_format": "root_hogs"}, {"schema": "unreviewed"}])
def test_invalid_scored_cell_rejected(retained, updates):
    write, prep, score = retained
    score.update(updates)
    with pytest.raises(ValueError, match="retained scored"):
        collect(write("preparation.json", prep), write("score.json", score))


@pytest.mark.parametrize("duplicate", [False, True])
def test_unbound_or_ambiguous_native_trace_rejected(retained, duplicate):
    write, prep, score = retained
    score["evidence"] = score["evidence"] * 2 if duplicate else []
    with pytest.raises(ValueError, match="uniquely bound"):
        collect(write("preparation.json", prep), write("score.json", score))


def test_trace_not_explaining_saved_partition_rejected(retained):
    write, prep, score = retained
    score["prediction"] = write("prediction.txt", "s x y a b\n")
    with pytest.raises(ValueError, match="reconstruct complete"):
        collect(write("preparation.json", prep), write("score.json", score))


def test_disagreement_with_scored_changed_genes_rejected(retained):
    write, prep, score = retained
    score["canonical_partition_comparison"]["genes_in_changed_groups"] = 0
    with pytest.raises(ValueError, match="disagree"):
        collect(write("preparation.json", prep), write("score.json", score))


def test_changed_trace_checksum_rejected(retained):
    write, prep, score = retained
    Path(score["evidence"][0]["path"]).write_text("[]")
    with pytest.raises(ValueError):
        collect(write("preparation.json", prep), write("score.json", score))


def test_actual_retained_report_reproduces_at_json_boundary():
    root = Path(__file__).resolve().parents[2]
    directory = root / "benchmark_tools/results/native_candidate_trace_replay_22429"
    path = directory / "diagnostic.json"
    if not path.exists():
        pytest.skip("Actual local trace diagnostic not installed")
    retained = json.loads(path.read_text())
    if any(not Path(ref["path"]).is_file() for ref in retained["inputs"]):
        pytest.skip("Retained local raw trace assets not installed")
    current = collect(record(root / "benchmark_tools/results/orthobench_factorial_prepared_20260916.json"),
                      record(root / "benchmark_tools/results/native_factorial_orthobench_score_22429.json"))
    # Semantic trace keys use tuples in memory; JSON arrays load as lists.
    current = json.loads(json.dumps(current, allow_nan=False))
    current["diagnostic_python"] = retained["diagnostic_python"]
    assert current == retained
    assert (directory / "diagnostic.md").read_text() == render(retained)
    assert retained["trace_comparison"]["common_semantic_merges"] == 8482
    assert retained["partition_comparison"]["genes_in_changed_groups"] == 58
    assert retained["partition_comparison"]["reference_genes_in_changed_groups"] == []
    assert len(retained["original_only_accepted_merges"]) == 18
    independent = json.loads((directory / "independent_readback.json").read_text())
    assert independent["report_sha256"] == record(path)["sha256"]
    assert independent["unchanged_native_helper_pins_checked"] == 920
    assert all(p["whole_partition_equal"] for p in independent["partitions"])
