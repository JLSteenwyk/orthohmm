import pytest

from benchmark_tools.audit_qfo_swiss_comparators import aggregate_verified
from benchmark_tools.audit_qfo_swiss_counts import statistics


def fixture():
    counts = {"A": {"TP": 8, "FP": 2, "FN": 20, "TN": 0},
              "B": {"TP": 2, "FP": 22, "FN": 0, "TN": 20}}
    values = [statistics(row) for row in counts.values()]
    p = sum(v["PPV"] for v in values) / 2
    r = sum(v["TPR"] for v in values) / 2
    return counts, {"metric_x": r, "metric_y": p}, 2 * p * r / (p + r)


def test_native_macro_harmonic_not_average_family_f1():
    counts, participant, f1 = fixture()
    values, means = aggregate_verified(counts, participant, f1)
    assert means["F1"] == pytest.approx(f1)
    assert means["F1"] != pytest.approx(sum(v["F1"] for v in values.values()) / 2)


@pytest.mark.parametrize("axis", ["metric_x", "metric_y"])
def test_native_mismatch_rejected(axis):
    counts, participant, f1 = fixture()
    participant[axis] += 0.001
    with pytest.raises(ValueError, match="endpoint"):
        aggregate_verified(counts, participant, f1)


def test_wrong_f1_rejected():
    counts, participant, f1 = fixture()
    with pytest.raises(ValueError, match="endpoint"):
        aggregate_verified(counts, participant, f1 + 0.01)
