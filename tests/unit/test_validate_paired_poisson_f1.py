"""Independent tail identities and nonlinear paired-rectangle checks."""

from itertools import product
import math

import numpy as np
import pytest
from scipy.stats import poisson

from benchmark_tools import validate_paired_poisson_f1 as candidate


@pytest.mark.parametrize("count", [0, 1, 2, 9, 50])
def test_poisson_tail_inversion_not_only_chi_square_agreement(count):
    alpha = .05 / 3
    lo, hi = candidate.garwood(count, alpha)
    assert poisson.cdf(count, hi) == pytest.approx(alpha / 2, abs=1e-13)
    if count:
        assert poisson.sf(count - 1, lo) == pytest.approx(alpha / 2, abs=1e-13)
    else:
        assert lo == 0 and hi == pytest.approx(-math.log(alpha / 2))


@pytest.mark.parametrize("truth", [1, 512, 23934])
def test_actual_ratio_difference_and_swapped_methods(truth):
    for c, a, b in ((1., 3., 0.), (4., .5, 8.), (0., 0., 0.)):
        expected = 2 * truth / (2 * truth + c + a) - 2 * truth / (2 * truth + c + b)
        assert candidate.contrast(truth, c, a, b) == pytest.approx(expected, abs=2e-16)
        assert candidate.contrast(truth, c, a, b) == -candidate.contrast(truth, c, b, a)
    lo, hi = candidate.interval(truth, (2, 0, 4))
    swap_lo, swap_hi = candidate.interval(truth, (2, 4, 0))
    assert swap_lo == pytest.approx(-hi) and swap_hi == pytest.approx(-lo)


@pytest.mark.parametrize("left,right", [((0, 2), (3, 5)), ((3, 5), (0, 2)), ((1, 4), (2, 3)), ((1, 1), (1, 1))])
def test_projected_extrema_match_all_eight_corners_and_cover_interiors(left, right):
    bounds = [(0, 7), left, right]
    lo, hi = candidate.project_box(512, bounds)
    values = [float(candidate.contrast(512, *point)) for point in product(*bounds)]
    assert lo == min(values) and hi == max(values)
    for point in product(*(np.linspace(a, b, 5) for a, b in bounds)):
        value = float(candidate.contrast(512, *point))
        assert lo - 1e-16 <= value <= hi + 1e-16


def test_zero_observed_disagreement_does_not_collapse_the_interval():
    lo, hi = candidate.interval(23934, (1, 0, 0))
    assert lo < 0 < hi
    assert lo <= candidate.contrast(23934, 1., 1., 0.) <= hi
    assert math.exp(-1) > .36  # Zero extra events retain positive probability.


@pytest.mark.parametrize("counts", [(-1, 0, 0), (1.5, 0, 0), (True, 0, 0), (2**53, 0, 0), (0, 0)])
def test_invalid_category_counts_are_refused(counts):
    with pytest.raises(ValueError):
        candidate.interval(512, counts)


@pytest.mark.parametrize("truth", [0, -1, 1.5, True, 2**53])
def test_invalid_truth_mass_is_refused(truth):
    with pytest.raises(ValueError):
        candidate.interval(truth, (0, 0, 0))


@pytest.mark.parametrize("alpha", [0, 1, -1, math.nan, math.inf, True])
def test_invalid_error_allocation_is_refused(alpha):
    with pytest.raises(ValueError):
        candidate.interval(512, (0, 0, 0), alpha)


def test_mean_and_box_validation():
    for means in ((-1, 0, 0), (math.inf, 0, 0), (math.nan, 0, 0), (0, 0)):
        with pytest.raises(ValueError):
            candidate.enumerate_cell(512, means)
    with pytest.raises(ValueError, match="Reversed"):
        candidate.project_box(512, [(1, 0), (0, 1), (0, 1)])


def test_vectorized_coverage_matches_independent_scalar_probability_sum():
    truth, means = 5, (.2, .3, .4)
    result = candidate.enumerate_cell(truth, means, tail_budget=1e-8)
    target = 2 * truth / (2 * truth + means[0] + means[1]) - 2 * truth / (2 * truth + means[0] + means[2])
    mass = covered = width_mass = 0.
    for counts in product(*(range(cutoff + 1) for cutoff in result["cutoffs"])):
        probability = math.prod(math.exp(-mean) * mean**count / math.factorial(count)
                                for mean, count in zip(means, counts))
        lo, hi = candidate.interval(truth, counts)
        mass += probability
        width_mass += probability * float(hi - lo)
        if lo <= target <= hi:
            covered += probability
    assert result["enumerated_mass"] == pytest.approx(mass, abs=1e-13)
    assert result["covered_enumerated_mass"] == pytest.approx(covered, abs=1e-13)
    assert result["mean_width_enumerated_mass"] == pytest.approx(width_mass, abs=1e-13)
    assert result["covered_enumerated_mass"] >= result["rectangle_covered_mass"] - 1e-14
    assert result["numerical_rounding_certified"] is False


def test_coverage_of_deterministic_zero_law_is_one():
    result = candidate.enumerate_cell(512, (0, 0, 0))
    assert result["enumerated_triples"] == 1
    assert result["covered_enumerated_mass"] == result["enumerated_mass"] == 1
    assert result["omitted_tail_mass"] == 0
    assert result["mean_width_enumerated_mass"] > 0


@pytest.mark.parametrize("truth", [512, 23934])
def test_common_shock_control_rejects_blanket_coverage_claim(truth):
    result = candidate.shock_counterexample(truth)
    assert sum(r["probability"] for r in result["outcomes"]) == 1
    assert result["mean_left_only"] == sum(r["probability"] * r["left_only_count"] for r in result["outcomes"])
    assert result["coverage"] == 0
    assert result["poisson_marginal_assumption_satisfied"] is False
    assert result["native_dependence_model"] is False


def test_cell_bound_refuses_before_large_allocation(monkeypatch):
    monkeypatch.setattr(candidate, "MAX_TRIPLES", 1)
    with pytest.raises(ValueError, match="Too many"):
        candidate.enumerate_cell(512, (1, 1, 1))


def test_changed_protocol_refuses_before_enumeration(tmp_path, monkeypatch):
    path = tmp_path / candidate.PROTOCOL
    path.parent.mkdir(parents=True)
    path.write_text("unit protocol")
    def forbidden(*args):
        pytest.fail("Should not enumerate after binding failure")
    monkeypatch.setattr(candidate, "enumerate_cell", forbidden)
    with pytest.raises(ValueError, match="Changed prospective"):
        candidate.run(tmp_path, "0" * 64)


def test_scope_flags_cannot_promote_conditional_success_to_native(tmp_path, monkeypatch):
    path = tmp_path / candidate.PROTOCOL
    path.parent.mkdir(parents=True)
    path.write_text("unit protocol")
    monkeypatch.setattr(candidate, "TRUTH_MASSES", (5,))
    monkeypatch.setattr(candidate, "MEANS", ((0., 0., 0.),))
    result = candidate.run(tmp_path, candidate.pin(path)["sha256"])
    assert result["all_conditional_cells_pass"] is True
    assert result["native_intervals_admitted"] is False
    assert result["overall_uncertainty_method_admitted"] is False
    assert result["independent_biological_confirmation"] is False
    assert "shared_clade" in result["excluded_laws"]


def test_failed_conditional_cell_stays_failed_without_native_admission(tmp_path, monkeypatch):
    path = tmp_path / candidate.PROTOCOL
    path.parent.mkdir(parents=True)
    path.write_text("unit protocol")
    monkeypatch.setattr(candidate, "TRUTH_MASSES", (5,))
    monkeypatch.setattr(candidate, "MEANS", ((0., 0., 0.),))
    monkeypatch.setattr(candidate, "enumerate_cell", lambda *args: dict(numerical_check_pass=False))
    result = candidate.run(tmp_path, candidate.pin(path)["sha256"])
    assert not result["all_conditional_cells_pass"] and not result["native_intervals_admitted"]


def test_existing_output_is_never_replaced(tmp_path, monkeypatch):
    output = tmp_path / "result.json"
    output.write_text("preserve")
    monkeypatch.setattr("sys.argv", ["test", "--root", str(tmp_path), "--protocol-sha256", "0" * 64, "--output", str(output)])
    with pytest.raises(FileExistsError):
        candidate.main()
    assert output.read_text() == "preserve"
