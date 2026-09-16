from copy import deepcopy

import numpy as np
import pytest

from benchmark_tools.bootstrap_orthobench import METRICS
from benchmark_tools.bootstrap_orthobench_factorial import CELLS, conditional_contrasts, factorial_bootstrap, render_report


def cells():
    return {cell: {"status": "complete", "refog_records": [
        {"refog": "small", "genes": 2, "true_positive": 1, "false_positive": 0, "false_negative": 0},
        {"refog": "large", "genes": 4, "true_positive": 1, "false_positive": 2, "false_negative": 5}]}
        for cell in CELLS}


def test_all_twelve_edges_flip_exactly_one_factor():
    contrasts = conditional_contrasts()
    assert len(contrasts) == len({(c["on"], c["off"]) for c in contrasts}) == 12
    for row in contrasts:
        assert sum(a != b for a, b in zip(row["on"], row["off"])) == 1
        assert len(row["fixed"]) == 2


def test_identical_cells_zero_intervals_and_order_invariance():
    fixture = cells()
    result = factorial_bootstrap(fixture, 100, 42)
    for contrast in result["comparisons"]:
        assert contrast["family_f1_ties"] == 2
        for metric in METRICS:
            assert contrast["metrics"][metric]["bonferroni_percentile_ci"] == [0, 0]
    shuffled = {key: {**value, "refog_records": value["refog_records"][::-1]} for key, value in reversed(list(fixture.items()))}
    assert factorial_bootstrap(shuffled, 100, 42) == result


def test_bootstrap_matches_independent_scalar_recomputation():
    fixture = cells()
    fixture["p1_c0_r0"]["refog_records"][1].update(true_positive=4, false_negative=2)
    result = factorial_bootstrap(fixture, 1000, 123)
    row = next(r for r in result["comparisons"] if r["on"] == "p1_c0_r0" and r["off"] == "p0_c0_r0")
    counts = np.random.Generator(np.random.PCG64(123)).multinomial(2, [0.5, 0.5], size=1000)
    def score(records, copies):
        ordered = sorted(records, key=lambda r: r["refog"])
        sums = [sum(n * r[key] / (r["genes"] - 1) for n, r in zip(copies, ordered))
                for key in ("true_positive", "false_positive", "false_negative")]
        tp, fp, fn = sums
        return [200 * tp / (2 * tp + fp + fn), 100 * tp / (tp + fp), 100 * tp / (tp + fn)]
    samples = np.array([[a - b for a, b in zip(score(fixture[row["on"]]["refog_records"], n),
                                               score(fixture[row["off"]]["refog_records"], n))] for n in counts])
    for i, metric in enumerate(METRICS):
        assert row["metrics"][metric]["paired_percentile_ci"] == pytest.approx(np.quantile(samples[:, i], [0.025, 0.975]))
        assert row["metrics"][metric]["bonferroni_percentile_ci"] == pytest.approx(np.quantile(samples[:, i], [0.05 / 72, 1 - 0.05 / 72]))
    observed = result["point_estimates_percent"]["p0_c0_r0"]["f_score"]
    assert observed == pytest.approx(score(fixture["p0_c0_r0"]["refog_records"], [1, 1])[0])
    assert observed != pytest.approx((100 + 200 / 9) / 2)


def test_failed_cell_removes_only_its_three_conditional_edges_without_zero_imputation():
    fixture = cells()
    fixture["p0_c0_r0"] = {"status": "failed", "reason": "native failure"}
    result = factorial_bootstrap(fixture, 100)
    assert sum(r["status"] == "unavailable" for r in result["comparisons"]) == 3
    assert result["multiplicity_endpoints"] == 36
    assert "p0_c0_r0" not in result["point_estimates_percent"]
    assert "| p0_c0_r0 | NA | NA | NA |" in render_report(result)


def test_all_failures_are_explicit_and_interactions_unavailable():
    result = factorial_bootstrap({cell: {"status": "failed", "reason": "failure"} for cell in CELLS}, 100)
    assert result["families"] == [] and result["point_estimates_percent"] == {}
    assert result["draws_generated"] == 0
    assert len(result["descriptive_interactions"]) == 6
    assert all(row["metrics"] is None for row in result["comparisons"])
    assert all(row["status"] == "unavailable" for row in result["descriptive_interactions"])


def test_interaction_sign_is_on_on_minus_on_off_minus_off_on_plus_off_off():
    fixture = cells()
    fixture["p1_c1_r1"]["refog_records"][1].update(true_positive=4, false_negative=2)
    result = factorial_bootstrap(fixture, 100)
    points = result["point_estimates_percent"]
    difference = points["p1_c1_r1"]["f_score"] - points["p0_c0_r0"]["f_score"]
    assert difference > 0
    for row in result["descriptive_interactions"]:
        expected = difference if next(iter(row["fixed"].values())) else 0
        assert row["difference_of_differences_percentage_points"]["f_score"] == pytest.approx(expected)
        assert "confidence_interval" not in row
    text = render_report(result)
    assert "Descriptive Family Counts" in text and "no interaction confidence intervals" in text


@pytest.mark.parametrize("change", ["missing", "pending", "failed_score", "family", "size"])
def test_invalid_panels_rejected(change):
    fixture = deepcopy(cells())
    if change == "missing":
        fixture.pop("p0_c0_r0")
    elif change == "pending":
        fixture["p0_c0_r0"]["status"] = "running"
    elif change == "failed_score":
        fixture["p0_c0_r0"].update(status="failed", reason="failure")
    elif change == "family":
        fixture["p0_c0_r0"]["refog_records"][0]["refog"] = "wrong"
    else:
        fixture["p0_c0_r0"]["refog_records"][0].update(genes=3, false_negative=2)
    with pytest.raises(ValueError):
        factorial_bootstrap(fixture, 100)
