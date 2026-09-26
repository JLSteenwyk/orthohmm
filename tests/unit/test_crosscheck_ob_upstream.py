import pytest

from benchmark_tools.crosscheck_ob_upstream import compare


def test_preserves_metric_order_and_full_precision():
    expected = dict(f_score=50.123456789, precision=40.234567891, recall=60.345678912)
    assert compare(list(expected.values()), expected) == dict(f_score=0, precision=0, recall=0)
    assert compare([expected["f_score"] + 1e-12, expected["precision"], expected["recall"]], expected)["f_score"] > 0


@pytest.mark.parametrize("values", [[50, 40], [50, 40, 60, 1], [float("nan"), 40, 60],
                                    [50, float("inf"), 60], [40, 50, 60], [50.0001, 40, 60]])
def test_rejects_missing_nonfinite_or_different_scores(values):
    with pytest.raises(ValueError):
        compare(values, dict(f_score=50, precision=40, recall=60))
