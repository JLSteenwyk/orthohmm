import numpy as np
import pytest

from benchmark_tools import validate_dyadic_f1_uncertainty as module


def test_gradient_finite_difference():
    counts = np.array([43., 17., 9.])
    _, gradient = module.f1_gradient(counts)
    for i in range(3):
        step = np.eye(3)[i] * 1e-4
        numeric = (module.f1_gradient(counts + step)[0] - module.f1_gradient(counts - step)[0]) / 2e-4
        assert numeric == pytest.approx(gradient[i], abs=1e-10)


def test_overlap_sum_matches_brute_force_with_diagonals():
    left, right = np.triu_indices(5, 1)
    rng = np.random.default_rng(123)
    single, edge = rng.normal(size=5), rng.normal(size=len(left))
    sets = [{i} for i in range(5)] + [{i, j} for i, j in zip(left, right)]
    values = np.concatenate([single, edge])
    brute = sum(values[i] * values[j] for i in range(len(values))
                for j in range(len(values)) if sets[i] & sets[j])
    assert module.shared_endpoint_variance(single, edge, left, right) == pytest.approx(brute)


def test_identical_methods_cancel():
    left, right = np.triu_indices(8, 1)
    diagonal, edges, _ = module.sample(np.random.default_rng(1), 8, "node_sparse", left, right)
    diagonal[:, 1] = diagonal[:, 0]
    edges[:, 1] = edges[:, 0]
    assert module.contrast(diagonal, edges, left, right) == (0., 0., 0.)


def test_method_swap_preserves_variance():
    left, right = np.triu_indices(8, 1)
    diagonal, edges, _ = module.sample(np.random.default_rng(4), 8, "node_sparse", left, right)
    a = module.contrast(diagonal, edges, left, right)
    b = module.contrast(diagonal[:, ::-1], edges[:, ::-1], left, right)
    assert a[0] == pytest.approx(-b[0])
    assert a[1:] == pytest.approx(b[1:])


def test_missing_zero_dyads_rejected():
    left, right = np.triu_indices(8, 1)
    diagonal, edges, _ = module.sample(np.random.default_rng(4), 8, "node_sparse", left, right)
    with pytest.raises(ValueError, match="complete"):
        module.contrast(diagonal, edges[:-1], left[:-1], right[:-1])


@pytest.mark.parametrize("counts", [[0, 0, 0], [-1, 1, 1], [np.nan, 1, 1]])
def test_invalid_counts(counts):
    with pytest.raises(ValueError):
        module.f1_gradient(counts)


def test_small_run_deterministic_and_not_admitted(monkeypatch):
    monkeypatch.setattr(module, "SIZES", (8,))
    a, b = module.run(5), module.run(5)
    assert a == b
    assert not a["benchmark_intervals_admitted"]
    assert len(a["rows"]) == 4
