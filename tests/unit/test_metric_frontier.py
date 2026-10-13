import numpy as np
import pytest

from benchmark_tools.metric_frontier import ENDPOINTS, metric_frontier


def vector(value=.8):
    return dict.fromkeys(ENDPOINTS, value)


def test_tradeoffs_and_ties_retained_without_mean_ranking():
    a, b, c = vector(), vector(), vector(.7)
    b[ENDPOINTS[0]] = .9
    b[ENDPOINTS[1]] = .75
    report = metric_frontier({"a": a, "a_tie": a, "b": b, "c": c})
    assert report["nondominated_complete_labels"] == ["a", "a_tie", "b"]
    assert report["dominated_by"]["c"] == ["a", "a_tie", "b"]


def test_missing_metrics_are_not_zero_or_eligible():
    report = metric_frontier({"partial": {ENDPOINTS[0]: 1.}, "complete": vector(), "absent": {}})
    assert report["nondominated_complete_labels"] == ["complete"]
    assert len(report["excluded_incomplete"]["partial"]) == 7
    assert len(report["excluded_incomplete"]["absent"]) == 8


@pytest.mark.parametrize("value", [-.1, 1.01, np.nan, np.inf, True])
def test_invalid_values_rejected(value):
    scores = vector()
    scores[ENDPOINTS[0]] = value
    with pytest.raises(ValueError, match="Invalid"):
        metric_frontier({"bad": scores})


def test_wrong_eight_endpoint_names_cannot_qualify():
    scores = vector()
    scores["unknown"] = scores.pop(ENDPOINTS[0])
    with pytest.raises(ValueError, match="Unknown"):
        metric_frontier({"wrong": scores})


@pytest.mark.parametrize("endpoints", [[], ["same", "same"]])
def test_empty_or_repeated_endpoint_definitions_rejected(endpoints):
    with pytest.raises(ValueError, match="unique endpoints"):
        metric_frontier({}, endpoints)
