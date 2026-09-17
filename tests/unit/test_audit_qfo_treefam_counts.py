import pytest

from benchmark_tools.audit_qfo_treefam_counts import verify_counts


def endpoint():
    return {"native_participant": {"metric_x": .8, "metric_y": 2 / 3}, "score": 8 / 11}


def test_one_pooled_prior_not_one_prior_per_relation():
    scores = verify_counts({"TP": 6, "FP": 2, "FN": 0, "TN": 5}, endpoint())
    assert scores["TPR"] == .8
    assert scores["PPV"] == 2 / 3
    assert scores["F1"] == pytest.approx(8 / 11)


@pytest.mark.parametrize("field", ["metric_x", "metric_y", "score"])
def test_changed_native_metric_rejected(field):
    native = endpoint()
    if field == "score":
        native[field] = .5
    else:
        native["native_participant"][field] = .5
    with pytest.raises(ValueError):
        verify_counts({"TP": 6, "FP": 2, "FN": 0, "TN": 5}, native)
