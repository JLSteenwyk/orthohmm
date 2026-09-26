import pytest

from benchmark_tools.audit_dyadic_f1_screen import oracle


def test_independent_diagonal_variance_by_enumeration():
    n = 64
    pairs = n * (n - 1) / 2
    fp_a, fp_b = 8 * pairs / n, 5 * pairs / n
    a = 2 * (20 * n + fp_a) / (35 * n + fp_a)**2
    b = 2 * (20 * n + fp_b) / (36 * n + fp_b)**2
    outcomes = [(a - b, .75), (-b, .05), (0., .20)]
    mean = sum(value * prob for value, prob in outcomes)
    variance = 20 * n * sum(prob * (value - mean)**2 for value, prob in outcomes)
    result = oracle(n, "independent_sparse")
    assert result["conditional_diagonal_variance"] == pytest.approx(variance)
    assert result["latent_effect_variance"] == 0


def test_clade_multiplier_only_changes_latent_variance():
    independent = oracle(64, "node_sparse")
    clade = oracle(64, "clade_sparse")
    assert clade["target"] == independent["target"]
    assert clade["latent_effect_variance"] == pytest.approx(8 * independent["latent_effect_variance"])
    assert clade["conditional_diagonal_variance"] == independent["conditional_diagonal_variance"]
    assert clade["conditional_poisson_variance"] == independent["conditional_poisson_variance"]


def test_unsupported_design_rejected():
    with pytest.raises(ValueError):
        oracle(12, "node_sparse")
