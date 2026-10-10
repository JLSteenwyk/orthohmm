from fractions import Fraction
from itertools import combinations
import math

import pytest

from benchmark_tools import native_fas_two_sample_interval as current
from benchmark_tools.native_fas_sampling_interval import difference_interval, expected_native_mean


@pytest.mark.parametrize("pre,missing", [
    ([0., .2, 1.], [None, .1, .8, None]),
    ([0., 1., 1., 1.], [None] * 4),
    ([.1], [.9]),
    ([.1, .9], [None, 0., 1.]),
    ([0., 0., 1., 1.], [0., 1., 1., None]),
])
@pytest.mark.parametrize("error", [.2, .05, .05 / 12])
def test_unknown_both_means_coverage_over_actual_two_stratum_subsets(pre, missing, error):
    # All native subsets are enumerated, not independent resamples of gene pairs.
    G = sum(x is not None for x in missing)
    muP = sum(pre) / len(pre)
    muG = sum(x for x in missing if x is not None) / G if G else 0.
    for k in range(1, len(pre) + 1):
        for c in range(len(missing) + 1):
            target = expected_native_mean(len(missing), G, c, k, muP, muG)
            total = covered = 0
            for old in combinations(pre, k):
                for selected in combinations(missing, c):
                    new = [x for x in selected if x is not None]
                    interval = current.method_interval(len(missing), c, len(pre), k,
                        sum(old) / k, len(new), sum(new) / len(new) if new else None, error)
                    lo, hi = interval["expected_native_mean_bounds"]
                    covered += lo - 1e-14 <= target <= hi + 1e-14
                    total += 1
            assert Fraction(covered, total) >= Fraction(str(1 - 3 * error))


def test_corners_cover_every_accepted_success_count_and_mean_endpoint():
    for r in range(5):
        result = current.method_interval(6, 4, 20, 9, .2, r, .7 if r else None, .1)
        lo, hi = result["success_count_bounds"]
        values = [expected_native_mean(6, g, 4, 9, a, b)
                  for g in range(lo, hi + 1)
                  for a in result["precomputed_mean_bounds"]
                  for b in result["return_mean_bounds"]]
        assert result["expected_native_mean_bounds"] == pytest.approx([min(values), max(values)], abs=2e-15)


def test_precomputed_census_and_zero_missing_draws():
    assert current.method_interval(0, 0, 2, 2, .4, 0, None, .01)["expected_native_mean_bounds"] == [.4, .4]
    partial = current.method_interval(0, 0, 100, 20, .4, 0, None, .01)
    assert partial["expected_native_mean_bounds"] == partial["precomputed_mean_bounds"]
    assert partial["expected_native_mean_bounds"][0] < .4 < partial["expected_native_mean_bounds"][1]


def test_added_uncertainty_does_not_condition_on_unknown_precomputed_population_mean():
    from benchmark_tools.native_fas_sampling_interval import method_interval as known
    old = known(100, 30, 10, .4, 28, .6, .05 / 12)
    new = current.method_interval(100, 30, 200, 10, .4, 28, .6, .05 / 12)
    assert new["expected_native_mean_bounds"][0] < old["expected_native_mean_bounds"][0]
    assert new["expected_native_mean_bounds"][1] > old["expected_native_mean_bounds"][1]


def test_joint_allocation_and_difference_projection_allow_arbitrary_dependence():
    error = Fraction(1, 240)
    assert 4 * 3 * error == Fraction(1, 20)
    rows = [current.method_interval(6, 4, 20, 9, value, 4, .7, float(error))
            for value in (.1, .3, .6, .9)]
    for a, b in combinations(rows, 2):
        lo, hi = difference_interval(a, b)
        assert difference_interval(b, a) == [-hi, -lo]


@pytest.mark.parametrize("size,mean,error", [
    (True, .1, .05), (-1, None, .05), (1.5, .4, .05), (0, .5, .05),
    (1, None, .05), (1, math.nan, .05), (2, 1.1, .05), (2, .4, 0), (2, .4, math.inf),
])
def test_invalid_sample_inputs(size, mean, error):
    with pytest.raises(ValueError):
        current.bounded_mean_interval(size, mean, error)


@pytest.mark.parametrize("P,k", [(0, 0), (20, 0), (1, 2), (2.0, 1), (2, True)])
def test_invalid_precomputed_counts(P, k):
    with pytest.raises(ValueError):
        current.method_interval(4, 2, P, k, .4, 1, .5, .05)
