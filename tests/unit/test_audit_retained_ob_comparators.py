import pytest

from benchmark_tools.audit_retained_ob_comparators import agrees


def test_accepts_retained_rounding_not_material_difference():
    score = dict(f_score=55.065342292123, precision=59.078055333123,
                 recall=51.563064181123, exact_refogs=12)
    expected = {k + "_percent": round(v, 9) for k, v in score.items() if k != "exact_refogs"}
    expected["exact_refogs"] = 12
    assert agrees(score, expected)
    assert not agrees(dict(score, f_score=55.07), expected)


@pytest.mark.parametrize("key,value", [
    ("f_score", float("nan")), ("precision", float("inf")),
    ("recall", 0), ("exact_refogs", 13),
])
def test_rejects_metric_or_exact_group_disagreement(key, value):
    score = dict(f_score=50, precision=50, recall=50, exact_refogs=12)
    expected = dict(f_score_percent=50, precision_percent=50, recall_percent=50, exact_refogs=12)
    assert not agrees(dict(score, **{key: value}), expected)
