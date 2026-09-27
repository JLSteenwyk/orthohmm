from types import SimpleNamespace

import numpy as np
import pytest

from benchmark_tools.probe_installed_ob_graph import align_scores, edge_comparison


def test_alignment_preserves_both_orders():
    old = (np.array([0, 1]), np.array([1, 2]), np.array([1., 2.]))
    new = (np.array([1, 0]), np.array([2, 1]), np.array([3., 4.]))
    a, b = align_scores(3, *old, *new)
    np.testing.assert_array_equal(a, [2., 1.])
    np.testing.assert_array_equal(b, [4., 3.])


@pytest.mark.parametrize("q,t,s", [([0, 0], [1, 1], [1., 1.]),
    ([0], [2], [1.]), ([0], [0], [1.]), ([-1], [1], [1.]),
    ([0], [1], [float("nan")]), ([0], [1], [-1.])])
def test_invalid_hit_sets(q, t, s):
    old = (np.array([0]), np.array([1]), np.array([1.]))
    with pytest.raises(ValueError):
        align_scores(3, *old, np.array(q), np.array(t), np.array(s))


def graph(q, t, s):
    return SimpleNamespace(sources=np.array(q), targets=np.array(t), weights=np.array(s))


def test_edge_comparison_with_order_and_weights():
    a = graph([0, 1], [1, 2], [1., 2.])
    b = graph([1, 0], [2, 2], [2.5, 4.])
    r = edge_comparison(a, b, 3)
    assert (r["left_only"], r["right_only"], r["shared_edges"]) == (1, 1, 1)
    assert r["shared_weights_changed"] == 1
    assert r["max_shared_weight_difference"] == .5


@pytest.mark.parametrize("a", [graph([1], [0], [1.]), graph([0, 0], [1, 1], [1., 1.])])
def test_invalid_edges(a):
    with pytest.raises(ValueError):
        edge_comparison(a, graph([0], [1], [1.]), 3)
