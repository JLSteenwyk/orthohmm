import copy

import pytest

from benchmark_tools import trace_ob_families as trace


def test_overlapping_reference_assignments_are_preserved():
    refs = {"family1": {"a", "b"}, "family2": {"a", "c"}}
    assert trace.reference_inventory(refs) == {"families": 2, "memberships": 4, "unique_genes": 3,
                                              "shared_genes": {"a": ["family1", "family2"]}}
    groups = {"all": {"a", "b", "c"}}
    index = {gene: "all" for gene in "abc"}
    for genes in refs.values():
        result = trace.family_group_summary(genes, groups, index, set("abc"))
        assert result["within_family_pairs"] == 1
        assert result["cross_reference_pairs_incident"] == 2


def test_partition_rejects_duplicate_and_incomplete_inputs(tmp_path):
    path = tmp_path / "groups.txt"
    for text in ("a b\na\n", "a\n"):
        path.write_text(text)
        with pytest.raises(ValueError):
            trace.partition(path, "plain", {"a", "b"})
    path.write_text("a b\nc\n")
    groups, index = trace.partition(path, "plain", {"a", "b", "c"})
    assert groups == {"0": {"a", "b"}, "1": {"c"}}
    assert index == {"a": "0", "b": "0", "c": "1"}


def test_family_summary_distinguishes_other_reference_and_unlabelled():
    groups = {"g": {"a", "b", "other", "u"}, "h": {"c", "v"}}
    index = {gene: key for key, members in groups.items() for gene in members}
    row = trace.family_group_summary({"a", "b", "c"}, groups, index, {"a", "b", "c", "other"})
    assert row["within_family_pairs"] == 1
    assert row["cross_reference_pairs_incident"] == 2
    assert row["unlabelled_pairs_incident"] == 3
    assert row["largest_family_component"] == 2
    assert row["groups_touching_family"] == 2


def test_pair_trace_directionality_and_nonmonotone_transitions():
    indices = {stage: {"a": "1", "b": "1" if i in (0, 2, 4) else "2", "c": "3"}
               for i, stage in enumerate(trace.STAGES)}
    rows = trace.pair_trace({"a", "b", "c"}, {("b", "a"): .5}, indices, {"a": "x", "b": "x", "c": "y"})
    assert len(rows) == 3
    assert rows[0]["forward_normalized_hit"] is None
    assert rows[0]["reverse_normalized_hit"] == .5
    assert rows[0]["same_species"] is True
    counts = trace.transitions(rows)
    assert counts["multipass_to_multipass_refined"] == {"lost": 1, "gained": 0, "retained": 0, "absent_both": 2}
    assert counts["multipass_refined_to_strict_profiles"]["gained"] == 1
    assert all(sum(row.values()) == 3 for row in counts.values())


@pytest.mark.parametrize("score", [0., -1., float("nan"), float("inf")])
def test_invalid_hit_score_rejected(score):
    indices = {stage: {"a": "1", "b": "1"} for stage in trace.STAGES}
    with pytest.raises(ValueError, match="hit score"):
        trace.pair_trace({"a", "b"}, {("a", "b"): score}, indices, {"a": "x", "b": "y"})


def merge(source, target, iteration=0):
    return {"source_genes": source, "target_genes": target, "source_size": len(source), "target_size": len(target),
            "iteration": iteration, "support": .5, "margin": 1.5}


def test_merge_reconstruction_accepts_iteration_start_side_snapshots():
    events = [merge(["a"], ["b"]), merge(["c"], ["b"]), merge(["d"], ["a", "b", "c"], 1)]
    retained = trace.validate_merge_reconstruction(events, {g: {g} for g in "abcd"}, {"all": set("abcd")}, {"a"})
    assert [row["event_index"] for row in retained] == [0, 2]


@pytest.mark.parametrize("change", ["unknown", "duplicate", "overlap", "size", "order", "cycle", "final_partition"])
def test_invalid_merge_evidence_rejected(change):
    events = [merge(["a"], ["b"]), merge(["c"], ["a", "b"], 1)]
    candidates = {"all": set("abc")}
    if change == "unknown":
        events[0]["source_genes"] = ["z"]
    elif change == "duplicate":
        events[0]["source_genes"] = ["a", "a"]
    elif change == "overlap":
        events[0]["target_genes"] = ["a"]
    elif change == "size":
        events[0]["source_size"] = 2
    elif change == "order":
        events[0]["iteration"] = 2
    elif change == "cycle":
        events.insert(1, copy.deepcopy(events[0]))
    elif change == "final_partition":
        candidates = {"ab": set("ab"), "c": {"c"}}
    with pytest.raises(ValueError):
        trace.validate_merge_reconstruction(events, {g: {g} for g in "abc"}, candidates, {"a"})
