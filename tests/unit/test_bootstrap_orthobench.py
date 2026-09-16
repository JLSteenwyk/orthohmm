import numpy as np
import pytest

from benchmark_tools.bootstrap_orthobench import paired_bootstrap, render_report, statistics, weighted_records
from benchmark_tools.score_orthobench_partition import score_partition


def records():
    return [
        dict(refog="a", genes=2, true_positive=1, false_positive=0, false_negative=0),
        dict(refog="b", genes=4, true_positive=0, false_positive=0, false_negative=6),
    ]


def test_recomputes_weighted_statistic_not_mean_family_f1():
    _, _, weights = weighted_records(records())
    assert statistics(weights.sum(axis=0))[0] == pytest.approx(50)
    assert statistics(weights)[:, 0].mean() == pytest.approx(50)
    # Unequal family TP/FN yields different macro F1 and pooled statistic.
    values = records()
    values[1].update(true_positive=1, false_negative=5)
    _, _, weights = weighted_records(values)
    assert statistics(weights.sum(axis=0))[0] != pytest.approx(statistics(weights)[:, 0].mean())


def test_identical_methods_have_zero_paired_uncertainty_and_order_is_irrelevant():
    result = paired_bootstrap({"a": records(), "b": records()[::-1]}, "a", 100, 123)
    assert result["comparisons"]["b"]["metrics"]["f_score"]["paired_percentile_ci"] == [0, 0]
    assert result["comparisons"]["b"]["family_f1_ties"] == 2
    assert result == paired_bootstrap({"b": records()[::-1], "a": records()}, "a", 100, 123)


def test_statistic_matches_scorer_with_low_certainty_exclusions():
    score = score_partition([{"a", "b", "x"}, {"c", "d"}, {"e"}],
                            {"one": {"a", "b"}, "two": {"c", "d", "e"}}, {"one": {"x"}})
    _, _, weights = weighted_records(score["refog_records"])
    assert statistics(weights.sum(axis=0)) == pytest.approx([score[k] for k in ("f_score", "precision", "recall")])


def test_zero_predictions_have_finite_zero_statistics():
    assert np.array_equal(statistics(np.array([0., 0., 10.])), np.zeros(3))


def test_report_includes_uncertainty_caveat_and_family_counts():
    report = render_report(paired_bootstrap({"a": records(), "b": records()}, "a", 100))
    assert "not independent confirmation" in report
    assert "| b | 0 | 2 | 0 |" in report
    assert "[0.000, 0.000]" in report


def test_incompatible_families_rejected():
    other = records()
    other[0]["refog"] = "wrong"
    with pytest.raises(ValueError, match="identical RefOG"):
        paired_bootstrap({"a": records(), "b": other}, "a", 100)


@pytest.mark.parametrize("field,value", [("genes", 1), ("false_positive", -1), ("true_positive", float("nan")), ("false_negative", 999)])
def test_invalid_statistics_rejected(field, value):
    values = records()
    values[0][field] = value
    with pytest.raises(ValueError):
        weighted_records(values)
