import math

import numpy as np
import pytest

from benchmark_tools.probe_frozen_null_scores import (
    ENDPOINTS, LENGTHS, REGIMES, THRESHOLDS, approximate_e, critical_score,
    pairs, probabilities, run, tail,
)


@pytest.mark.parametrize("regime", REGIMES)
def test_sampling_is_reproducible_and_independent_of_band(regime):
    background = np.arange(1, 21)
    first = pairs(0, regime, 50, background, count=20)
    second = pairs(0, regime, 50, background, count=20)
    for a, b in zip(first, second):
        np.testing.assert_array_equal(a, b)
        assert a.shape == (20, 50) and a.dtype == np.uint8
    assert not np.array_equal(first[0], first[1])
    assert not np.array_equal(first[0], pairs(1, regime, 50, background, count=20)[0])
    assert math.isclose(probabilities(background, regime).sum(), 1)
    if regime == "half_glutamine":
        assert probabilities(background, regime)[13] > 0.5


@pytest.mark.parametrize("length", LENGTHS)
@pytest.mark.parametrize("threshold", THRESHOLDS)
def test_strict_integer_boundary_preserves_frozen_e_gate(length, threshold):
    score = critical_score(length, threshold)
    assert approximate_e(score, length) < threshold
    assert approximate_e(score - 1, length) >= threshold


def test_null_tail_distinguishes_poisson_reference_and_observed_frequency():
    length, threshold = 150, 0.1
    boundary = critical_score(length, threshold)
    result = tail(np.array([0, boundary - 1, boundary, boundary + 1]), length, threshold)
    assert result["hits"] == 2 and result["trials"] == 4 and result["fraction"] == 0.5
    assert result["boundary_approximate_e"] < threshold
    assert result["poisson_model_tail_reference"] < threshold
    assert math.isclose(result["poisson_model_tail_reference"], -math.expm1(-result["boundary_approximate_e"]))
    assert result["bonferroni_clopper_pearson"][0] < result["nominal_clopper_pearson"][0]
    assert result["bonferroni_clopper_pearson"][1] > result["nominal_clopper_pearson"][1]


def test_zero_hits_does_not_establish_zero_probability():
    result = tail(np.zeros(10000, dtype=np.int32), 150, 1e-4)
    assert result["hits"] == 0 and result["fraction"] == 0
    assert result["nominal_clopper_pearson"][0] == 0
    assert math.isclose(result["nominal_clopper_pearson"][1], 1 - 0.025 ** (1 / 10000), abs_tol=1e-12)
    assert result["bonferroni_clopper_pearson"][1] > 1e-4
    assert ENDPOINTS == 90


@pytest.mark.parametrize("scores", [[], [1.5], [-1], [[1, 2]]])
def test_invalid_scores_fail_closed(scores):
    with pytest.raises(ValueError):
        tail(scores, 150, 0.1)


def test_output_is_not_overwritten(tmp_path):
    with pytest.raises(FileExistsError):
        run(tmp_path, tmp_path / "missing-protocol", tmp_path)
