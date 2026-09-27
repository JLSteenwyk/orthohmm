import numpy as np
import pytest

from benchmark_tools.audit_qfo_candidate_order import order_statistics


@pytest.mark.parametrize("chunk", [1, 2, 3, 20])
def test_sorted_unique(chunk):
    r = order_statistics(np.array([0, 0, 1]), np.array([0, 2, 0]), np.ones(3), 3, chunk)
    assert r["canonical_input_unchanged"] and r["global_uniqueness_proven"]
    assert r["self_hits"] == 1


@pytest.mark.parametrize("chunk", [1, 2, 3])
def test_unsorted_and_duplicate_boundaries(chunk):
    r = order_statistics(np.array([1, 0, 0]), np.array([0, 2, 2]), np.ones(3), 3, chunk)
    assert r["adjacent_descents"] == 1
    assert r["adjacent_duplicates"] == 1
    assert not r["canonical_input_unchanged"]


@pytest.mark.parametrize("scores", [[1, float('nan')], [1, 0], [1, -1]])
def test_invalid_scores(scores):
    with pytest.raises(ValueError):
        order_statistics(np.array([0, 1]), np.array([1, 0]), np.array(scores), 2)


def test_invalid_index_and_shape():
    with pytest.raises(ValueError):
        order_statistics(np.array([0, 2]), np.array([1, 0]), np.ones(2), 2)
    with pytest.raises(ValueError):
        order_statistics(np.array([0]), np.array([1, 0]), np.ones(2), 2)


def test_target_descent_and_nonadjacent_duplicate():
    r = order_statistics(np.array([0, 0, 0]), np.array([1, 0, 1]), np.ones(3), 2, 1)
    assert r["adjacent_descents"] == 1 and r["adjacent_duplicates"] == 0
    assert not r["global_uniqueness_proven"]


def test_empty_arrays():
    r = order_statistics(np.array([], dtype=int), np.array([], dtype=int), np.array([]), 1)
    assert r["hits"] == 0 and r["canonical_input_unchanged"]


def test_matches_independent_python_pair_comparison():
    rng = np.random.default_rng(1234)
    q, t = rng.integers(0, 20, size=(2, 123))
    pairs = list(zip(q.tolist(), t.tolist()))
    expected = sum(b < a for a, b in zip(pairs, pairs[1:]))
    duplicates = sum(b == a for a, b in zip(pairs, pairs[1:]))
    for chunk in (1, 7, 123, 200):
        r = order_statistics(q, t, np.ones(123), 20, chunk)
        assert r["adjacent_descents"] == expected
        assert r["adjacent_duplicates"] == duplicates
