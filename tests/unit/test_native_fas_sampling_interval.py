from fractions import Fraction
from itertools import combinations
from math import comb

import pytest

from benchmark_tools import native_fas_sampling_interval as current


def mass(M, G, c, r):
    if not max(0, c - M + G) <= r <= min(c, G):
        return Fraction(0)
    return Fraction(comb(G, r) * comb(M - G, c - r), comb(M, c))


@pytest.mark.parametrize("error", [Fraction(1,2),Fraction(1,8),Fraction(1,16),Fraction(1,20)])
def test_integer_tail_inversion_and_count_coverage_against_exact_mass(error):
    for M in range(9):
        for c in range(M + 1):
            bounds = {}
            for r in range(c + 1):
                accepted = [G for G in range(M + 1)
                            if sum((mass(M, G, c, n) for n in range(r, c + 1)), Fraction()) >= error / 2
                            and sum((mass(M, G, c, n) for n in range(r + 1)), Fraction()) >= error / 2]
                assert accepted
                bounds[r] = current.count_interval(M, c, r, float(error))
                assert bounds[r] == (min(accepted), max(accepted))
            for G in range(M + 1):
                covered = sum((mass(M, G, c, r) for r in range(c + 1)
                               if bounds[r][0] <= G <= bounds[r][1]), Fraction())
                assert covered >= 1 - error


def test_inclusive_boundary_tails_and_full_population_counts():
    assert current.count_interval(32,1,1,1/16) == (1,32)
    assert current.count_interval(32,1,0,1/16) == (0,31)
    for r in range(9):
        assert current.count_interval(8,8,r,.05) == (r,r)


def test_random_denominator_target_matches_actual_subset_ratios_not_count_plugin():
    pre = [Fraction(1, 10), Fraction(3, 10), Fraction(8, 10)]
    missing = [None, Fraction(1, 5), Fraction(9, 10), None]
    ratios = []
    for old in combinations(pre, 2):
        for selected in combinations(missing, 2):
            new = [score for score in selected if score is not None]
            ratios.append((sum(old) + sum(new)) / (2 + len(new)))
    exact = sum(ratios) / len(ratios)
    assert exact == Fraction(107, 240) and exact != Fraction(9, 20)
    actual = current.expected_native_mean(4, 2, 2, 2, float(sum(pre)/3), .55)
    assert actual == pytest.approx(float(exact), abs=1e-15)


def test_weights_and_rectangle_corners_against_exact_full_support():
    for M in range(1, 9):
        for c in range(M + 1):
            weights = [current.mixture_weight(M, G, c, 3) for G in range(M + 1)]
            assert all(a >= b for a, b in zip(weights, weights[1:]))
            for G, weight in enumerate(weights):
                exact = sum((mass(M, G, c, r) * Fraction(3, 3 + r)
                             for r in range(c + 1)), Fraction())
                assert weight == pytest.approx(float(exact), abs=2e-15)
            for r in range(c + 1):
                interval = current.method_interval(M, c, 3, .4, r, .7 if r else None, .05)
                gl, gh = interval["success_count_bounds"]
                ml, mh = interval["return_mean_bounds"]
                values = [current.expected_native_mean(M, G, c, 3, .4, mean)
                          for G in range(gl, gh + 1) for mean in (ml, mh)]
                assert interval["expected_native_mean_bounds"] == pytest.approx([min(values), max(values)], abs=2e-15)


@pytest.mark.parametrize("missing", [[None]*4, [.2], [None,.2], [0.,1.,None,None], [.2,.8,.9,None,None]])
def test_conditional_target_coverage_over_actual_subsets_with_score_dependent_omissions(missing):
    M = len(missing)
    returns = [Fraction(str(x)) for x in missing if x is not None]
    G = len(returns)
    mean = float(sum(returns)/G) if G else 0.
    for c in range(M + 1):
        target = current.expected_native_mean(M, G, c, 2, .4, mean)
        covered = 0
        subsets = list(combinations(missing, c))
        for subset in subsets:
            selected = [x for x in subset if x is not None]
            interval = current.method_interval(M, c, 2, .4, len(selected),
                                               sum(selected)/len(selected) if selected else None, .025)
            lo, hi = interval["expected_native_mean_bounds"]
            covered += lo - 1e-14 <= target <= hi + 1e-14
        assert Fraction(covered, len(subsets)) >= Fraction(19, 20)


def test_zero_returns_are_not_zero_width_plugin_and_G_zero_G_one_boundaries():
    interval = current.method_interval(20, 5, 2, .4, 0, None, .025)
    assert interval["return_mean_bounds"] == [0, 1]
    assert interval["expected_native_mean_bounds"][0] < .4 < interval["expected_native_mean_bounds"][1]
    assert current.expected_native_mean(20, 0, 5, 2, .4, 1.) == .4
    assert current.method_interval(0, 0, 2, .4, 0, None, .025)["expected_native_mean_bounds"] == [.4, .4]
    assert current.expected_native_mean(1, 1, 1, 2, .4, .7) == pytest.approx(.5)


def test_difference_projection_and_swapping():
    left = current.method_interval(6, 3, 2, .2, 2, .4, .05/16)
    right = current.method_interval(5, 4, 7, .9, 4, .8, .05/16)
    lo, hi = current.difference_interval(left, right)
    assert current.difference_interval(right, left) == [-hi, -lo]


@pytest.mark.parametrize("M,c,r", [(-1,0,0),(4,5,1),(4,2,3),(4,2,-1),(4,2,True),(4.,2,1),(10000,9001,0)])
def test_invalid_count_support_refused(M,c,r):
    with pytest.raises(ValueError):
        current.count_interval(M,c,r,.05)


@pytest.mark.parametrize("error", [0,1,-.1,float("nan"),float("inf")])
def test_invalid_error_refused(error):
    with pytest.raises(ValueError):
        current.count_interval(4,2,1,error)


def test_invalid_score_and_missingness_inputs_refused():
    for r, mean in ((0,.2),(1,None),(1,float("nan")),(1,1.1)):
        with pytest.raises(ValueError):
            current.method_interval(4,2,2,.4,r,mean,.025)
    with pytest.raises(ValueError):
        current.mixture_weight(4,2,2,0)
