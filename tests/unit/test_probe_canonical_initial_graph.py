import numpy as np
import pytest

from benchmark_tools.probe_canonical_initial_graph import selection_trace, trace_comparison, same_arrays, verify_trace_edges
from orthohmm.accuracy import build_rbnh_edges, build_singleton_assignment_edges


def fixture():
    names = ["singleton", "a", "b", "c", "d"]
    species = np.arange(5, dtype=np.int32)
    clusters = [[0], [1, 2], [3, 4]]
    queries = np.array([1, 2, 3, 4, 0, 0, 0, 0], dtype=np.int32)
    targets = np.array([2, 1, 4, 3, 1, 3, 2, 4], dtype=np.int32)
    scores = np.array([1, 1, 1, 1, .8, .8, .7, .6], dtype=np.float64)
    order = np.array([0, 1, 2, 3, 5, 4, 6, 7])
    return names, species, clusters, queries, targets, scores, order


def test_matched_rbnh_graph_can_have_order_sensitive_singleton_attachment():
    names, species, clusters, q, t, s, order = fixture()
    graphs = [build_rbnh_edges(names, species, a, b, c) for a, b, c in ((q, t, s), (q[order], t[order], s[order]))]
    assert graphs[0].gene_names == graphs[1].gene_names
    for field in ("sources", "targets", "weights"):
        assert getattr(graphs[0], field).tobytes() == getattr(graphs[1], field).tobytes()
    edges = [build_singleton_assignment_edges(names, clusters, a, b, c)
             for a, b, c in ((q, t, s), (q[order], t[order], s[order]))]
    traces = [selection_trace(len(names), clusters, a, b, c)
              for a, b, c in ((q, t, s), (q[order], t[order], s[order]))]
    assert trace_comparison(*traces).tolist() == [0]
    assert traces[0][1][0] == 1 and traces[1][1][0] == 2
    assert traces[0][2][0] == traces[1][2][0] == 2
    assert set(edges[0].targets.tolist()) == {1, 2}
    assert set(edges[1].targets.tolist()) == {3, 4}
    for trace, output in zip(traces, edges):
        verify_trace_edges(len(names), clusters, trace, output)
    with pytest.raises(ValueError, match="actual singleton edges"):
        verify_trace_edges(len(names), clusters, traces[0], edges[1])


def test_duplicate_winning_genes_in_one_cluster_are_not_multicluster_ties():
    _, _, clusters, q, t, s, _ = fixture()
    s[-2] = .8
    s[5] = .7
    best, chosen, count = selection_trace(5, clusters, q, t, s)
    assert best[0] == .8 and chosen[0] == 1 and count[0] == 1


def test_nonpositive_nonfinite_and_self_cluster_hits_do_not_select_targets():
    q = np.array([0, 0, 0, 0, 1], dtype=np.int32)
    t = np.array([1, 2, 1, 0, 0], dtype=np.int32)
    scores = np.array([0, -1, np.nan, 1, 10], dtype=np.float64)
    trace = selection_trace(3, [[0], [1, 2]], q, t, scores)
    assert trace[1].tolist() == [-1, -1, -1]
    assert trace[2].tolist() == [0, 0, 0]


def test_empty_input_and_all_singletons_do_not_select_clusters():
    q = np.array([], dtype=np.int32)
    trace = selection_trace(2, [[0], [1]], q, q, np.array([], dtype=np.float64))
    assert trace[1].tolist() == [-1, -1]
    assert trace_comparison(trace, trace).size == 0


@pytest.mark.parametrize("field", [0, 2])
def test_score_or_tie_distribution_difference_rejects_matched_claim(field):
    left = (np.array([1.]), np.array([1]), np.array([2]))
    right = list(left)
    right[field] = np.array([3])
    with pytest.raises(ValueError, match="distributions differ"):
        trace_comparison(left, right)


def test_changed_selection_without_a_tie_rejected():
    left = (np.array([1.]), np.array([1]), np.array([1]))
    right = (np.array([1.]), np.array([2]), np.array([1]))
    with pytest.raises(ValueError, match="no multi-cluster"):
        trace_comparison(left, right)


def test_exact_graph_comparison_requires_all_array_bytes_and_hashes():
    left = {"count": 2, "arrays": {k: {"bytes": 10, "sha256": "same"} for k in ("sources", "targets", "weights")}}
    right = {"count": 2, "arrays": {k: dict(v) for k, v in left["arrays"].items()}}
    assert same_arrays(left, right)
    right["arrays"]["weights"]["sha256"] = "different"
    assert not same_arrays(left, right)
