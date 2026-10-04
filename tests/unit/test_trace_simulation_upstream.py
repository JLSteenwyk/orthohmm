"""Stage-boundary tests; graph connectivity is not equated with orthology."""

import pytest

from benchmark_tools import trace_simulation_upstream as trace


def test_components_include_isolates_and_undirected_paths(tmp_path):
    path = tmp_path / "edges"
    path.write_text("a\tb\t1\nc\tb\t2\n")
    edges, labels, count = trace.graph(path, ["a", "b", "c", "d"])
    assert edges == {("a", "b"), ("b", "c")}
    assert count == 2 and labels["a"] == labels["c"] != labels["d"]


@pytest.mark.parametrize("data", ["a\tx\t1\n", "a\ta\t1\n", "a\tb\t0\n",
    "a\tb\t-1\n", "a\tb\tnan\n", "a\tb\tinf\n", "a\tb\t1\nb\ta\t1\n", "a\tb\n"])
def test_invalid_native_graph_rejects(tmp_path, data):
    path = tmp_path / "edges"
    path.write_text(data)
    with pytest.raises(ValueError):
        trace.graph(path, ["a", "b"])


def test_empty_graph_has_one_component_per_gene(tmp_path):
    path = tmp_path / "edges"
    path.write_text("")
    edges, labels, count = trace.graph(path, ["a", "b", "c"])
    assert not edges and count == 3 and len(set(labels.values())) == 3


def test_direct_hit_absence_does_not_mean_no_graph_path_or_prediction():
    owners, ancestors = {"a": "A", "b": "B", "c": "C"}, {"a": "1", "b": "1", "c": "2"}
    rows = list(trace.trace_pairs([("a", "b")], owners, ancestors, set(), {("a", "c"), ("b", "c")},
        {g: 0 for g in owners}, {g: "f" for g in owners}, {("a", "b")}))
    assert rows[0]["graph_connected"] and rows[0]["native_predicted"]
    assert not rows[0]["hit_forward"] and not rows[0]["graph_direct"]
    assert trace.summarize(rows)["totals"]["native_fn"] == 0


def test_direct_hit_and_graph_edge_can_still_be_separated_by_grouping():
    owners, ancestors = {"a": "A", "b": "B"}, {"a": "1", "b": "1"}
    rows = list(trace.trace_pairs([("b", "a")], owners, ancestors, {("a", "b")}, {("a", "b")},
        {"a": 0, "b": 0}, {"a": "f1", "b": "f2"}, set()))
    result = trace.summarize(rows)["totals"]
    assert result["connected_but_separated"] == result["direct_graph_edge_but_separated"] == 1
    assert result["connected_but_separated_hit_orientations_1"] == 1


def test_candidate_expansion_can_bridge_retained_graph_components():
    owners, ancestors = {"a": "A", "b": "B"}, {"a": "1", "b": "1"}
    rows = list(trace.trace_pairs([("a", "b")], owners, ancestors, set(), set(),
        {"a": 0, "b": 1}, {"a": "f", "b": "f"}, set()))
    result = trace.summarize(rows)["totals"]
    assert result["native_fn_within_candidates"] == 1 and result["across_candidates"] == 0


@pytest.mark.parametrize("pairs", [[("a", "a")], [("a", "b"), ("b", "a")], [("a", "x")]])
def test_malformed_or_duplicate_truth_rejects(pairs):
    owners, ancestors = {"a": "A", "b": "B"}, {"a": "1", "b": "1"}
    with pytest.raises(ValueError):
        list(trace.trace_pairs(pairs, owners, ancestors, set(), set(), {"a": 0, "b": 0}, {"a": "f", "b": "f"}, set()))


@pytest.mark.parametrize("field", ["graph_direct", "native_predicted"])
def test_impossible_stage_combinations_reject(field):
    row = {**{k: False for k in trace.FLAGS}, "ancestor": "1"}
    row[field] = True
    with pytest.raises(ValueError):
        trace.summarize([row])


def test_seed_sidecar_does_not_invent_gene_to_seed_partition(tmp_path):
    path = tmp_path / "seeds.tsv"
    path.write_text("candidate_family\tseed_families\nf\tSeed0000000,Seed0000001\n")
    result = trace.seed_sidecar(path, {"f": {"a", "b", "c"}},
        {"seed_families": 2, "candidate_families": 1, "merges": 1},
        [{"source_genes": ["a"], "target_genes": ["b", "c"]}])
    assert result["merges"] == 1 and result["gene_to_seed_inside_merged_candidates_retained"] is False


@pytest.mark.parametrize("data", ["f\tSeed0000000,Seed0000000\n", "f\tSeed-000001\n",
    "f\tSeed0000000\ng\tSeed0000000\n", "x\tSeed0000000\n"])
def test_invalid_seed_sidecar_rejects(tmp_path, data):
    path = tmp_path / "seeds.tsv"
    path.write_text("candidate_family\tseed_families\n" + data)
    with pytest.raises(ValueError):
        trace.seed_sidecar(path, {"f": {"a"}, "g": {"b"}},
            {"seed_families": 2, "candidate_families": 2, "merges": 0}, [])


def test_merge_crossing_final_candidates_rejects(tmp_path):
    path = tmp_path / "seeds.tsv"
    path.write_text("candidate_family\tseed_families\nf\tSeed0000000,Seed0000001\ng\tSeed0000002\n")
    with pytest.raises(ValueError, match="boundary"):
        trace.seed_sidecar(path, {"f": {"a", "b"}, "g": {"c"}},
            {"seed_families": 3, "candidate_families": 2, "merges": 1},
            [{"source_genes": ["a"], "target_genes": ["c"]}])


def test_renderer_displays_zero_counts_omitted_by_counter_sum():
    report = {"summary": [{"condition": "baseline", "totals": {"true_pairs": 10}}], "limitations": []}
    text = trace.render(report)
    assert "| baseline | 10 | 0 | 0 | 0 | 0 | 0 |" in text
    assert "0 / 0 / 0" in text
