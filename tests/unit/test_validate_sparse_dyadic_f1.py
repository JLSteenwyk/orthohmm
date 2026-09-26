import numpy as np
import pytest

from benchmark_tools.validate_dyadic_f1_uncertainty import contrast
from benchmark_tools.validate_sparse_dyadic_f1 import sparse_contrast, sample


@pytest.mark.parametrize("case", ["unequal_regular", "unequal_node", "rare_perfect"])
def test_sparse_equals_full_grid(case):
    n = 16
    diagonal, edges, left, right, _ = sample(np.random.default_rng(37), n, case,
                                           np.array([1, 3, 6]), np.array([.8, .15, .05]))
    all_left, all_right = np.triu_indices(n, 1)
    dense = np.zeros((len(all_left), 2, 3))
    lookup = {(i, j): k for k, (i, j) in enumerate(zip(all_left, all_right))}
    for row, i, j in zip(edges, left, right):
        dense[lookup[i, j]] = row
    expected = contrast(diagonal, dense, all_left, all_right)
    observed = sparse_contrast(diagonal, edges, left, right)
    assert observed == pytest.approx(expected[:2], abs=1e-15)


def test_empty_edges_and_identical_methods():
    diagonal = np.ones((3, 2, 3))
    result = sparse_contrast(diagonal, np.zeros((0, 2, 3)), np.array([], dtype=int), np.array([], dtype=int))
    assert result == (0., 0.)


@pytest.mark.parametrize("left,right", [([0, 0], [1, 1]), ([1], [0]), ([-1], [2]), ([0], [3])])
def test_invalid_sparse_indices(left, right):
    with pytest.raises(ValueError):
        sparse_contrast(np.ones((3, 2, 3)), np.ones((len(left), 2, 3)), np.array(left), np.array(right))


def test_negative_counts_rejected():
    with pytest.raises(ValueError):
        sparse_contrast(-np.ones((3, 2, 3)), np.zeros((0, 2, 3)), np.array([], dtype=int), np.array([], dtype=int))
