"""Independent native-row projections without importing the trace worker."""

from collections import Counter

import pytest

from benchmark_tools import readback_simulation_upstream as readback


def context():
    return {"truth": {("a", "b")}, "ancestors": {"a": "1", "b": "1"},
        "hits": {("b", "a")}, "edges": {("a", "b")}, "components": {"a": 0, "b": 0},
        "candidates": {"a": 0, "b": 1}, "native": set()}


def row():
    return {"gene_a": "a", "gene_b": "b", "ancestor": "1",
        **dict(zip(readback.FLAGS, ("0", "1", "1", "1", "0", "0")))}


def test_reverse_hit_orientation_and_group_separation_are_preserved():
    assert readback.row_flags(row(), context()) == (0, 1, 1, 1, 0, 0)


@pytest.mark.parametrize("field", readback.FLAGS)
def test_each_stage_flag_is_independently_verified(field):
    value = row()
    value[field] = "1" if value[field] == "0" else "0"
    with pytest.raises(ValueError, match="stage flags"):
        readback.row_flags(value, context())


@pytest.mark.parametrize("value", ["true", "False", "", "2"])
def test_nonbinary_tsv_flag_rejects(value):
    item = row()
    item["hit_forward"] = value
    with pytest.raises(ValueError, match="Nonbinary"):
        readback.row_flags(item, context())


def test_ancestral_labels_are_verified_against_truth():
    item = row()
    item["ancestor"] = "2"
    with pytest.raises(ValueError, match="Ancestral"):
        readback.row_flags(item, context())


def test_nonreference_pair_rejects():
    item = row()
    item["gene_b"] = "c"
    with pytest.raises(ValueError, match="nonreference"):
        readback.row_flags(item, context())


def test_igraph_uses_all_vertices_and_undirected_paths(tmp_path):
    path = tmp_path / "edges"
    path.write_text("a\tb\t1\nc\tb\t1\n")
    edges, labels, count = readback.component_labels(path, ["a", "b", "c", "d"])
    assert count == 2 and labels["a"] == labels["c"] != labels["d"]
    assert edges == {("a", "b"), ("b", "c")}


def test_complete_tsv_projection_preserves_overlap_and_partition_boundaries():
    counter = Counter({(0, 0, 0, 0, 0, 0): 5, (1, 1, 1, 1, 0, 0): 2,
        (0, 1, 0, 1, 1, 0): 3, (1, 1, 1, 1, 1, 1): 7})
    totals = readback.project(counter)
    assert totals["true_pairs"] == 17
    assert totals["native_fn"] == 10 and totals["native_tp"] == 7
    assert totals["across_candidates"] == 7 and totals["native_fn_within_candidates"] == 3
    assert totals["different_graph_components_hit_orientations_0"] == 5
    assert totals["connected_but_separated_hit_orientations_2"] == 2
    assert totals["direct_graph_edge_but_separated"] == 2


def test_omitted_zero_counter_values_are_not_missing_positive_counts():
    assert readback.equal_counts({"a": 0, "b": 2}, {"b": 2})
    assert not readback.equal_counts({"a": 1, "b": 2}, {"b": 2})
