import sys

import numpy as np
import pytest

from benchmark_tools.trace_ob_initial_edges import capture_edges, classify
from orthohmm.accuracy import build_rbnh_edges


@pytest.mark.parametrize("f,r,a,b,edge,expected", [
    (None, None, 1, 2, False, "no_direct_hit"),
    (2, None, 2, 3, True, "initial_edge"),
    (None, 2, 3, 2, True, "initial_edge"),
    (1, 1, 2, 3, False, "below_endpoint_threshold"),
    (1, None, float("inf"), float("inf"), False, "no_finite_endpoint_threshold"),
])
def test_decisions(f, r, a, b, edge, expected):
    assert classify(f, r, a, b, edge) == expected


def test_edge_disagreement_rejected():
    with pytest.raises(ValueError, match="differs"):
        classify(2, None, 1, 1, False)


@pytest.mark.parametrize("value", [float("nan"), -1, 0])
def test_invalid_threshold(value):
    with pytest.raises(ValueError):
        classify(1, None, value, 1, True)


@pytest.mark.parametrize("seed", range(8))
def test_native_threshold_decision_matches_all_pairs(seed):
    rng = np.random.default_rng(seed)
    n = 12
    pairs = [(a, b) for a in range(n) for b in range(n) if rng.random() < .4]
    q = np.array([a for a, _ in pairs], dtype=np.int32)
    t = np.array([b for _, b in pairs], dtype=np.int32)
    scores = rng.integers(1, 5, len(pairs)).astype(float)
    previous = sys.gettrace()
    edges, thresholds = capture_edges(build_rbnh_edges, list(map(str, range(n))),
                                      np.arange(n, dtype=np.int32) % 3, q, t, scores)
    assert sys.gettrace() is previous
    hits = dict(zip(pairs, scores))
    present = {tuple(sorted((int(a), int(b)))) for a, b in zip(edges.sources, edges.targets)}
    for a in range(n):
        for b in range(a + 1, n):
            classify(hits.get((a, b)), hits.get((b, a)), thresholds[a], thresholds[b], (a, b) in present)


def test_capture_restored_after_native_failure():
    previous = sys.gettrace()
    with pytest.raises(ValueError):
        capture_edges(build_rbnh_edges, ["a"], np.array([0]), np.array([2]), np.array([0]), np.array([1.]))
    assert sys.gettrace() is previous
