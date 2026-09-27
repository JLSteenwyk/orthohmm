import copy

import pytest

from benchmark_tools.export_qfo_canonical_comparison import compare, METRICS


def inputs():
    assessment = dict(participant="ohmm_qfo_corrected_factorial_p1_c1_r1",
                      secondary_six_metric_mean=0.5,
                      endpoints={m: dict(score=0.5, axes={}, score_semantics=m,
                                         native_participant={}) for m in METRICS})
    canonical = copy.deepcopy(assessment)
    canonical["participant"] = "ohmm_qfo_canonical_20260927"
    return (dict(status="historical_comparator_readmitted", assessment=assessment),
            dict(status="canonical_qfo_assessment_admitted", accuracy_admitted=True, assessment=canonical))


def test_equal_endpoints_preserve_fas_caveat():
    rows = compare(*inputs())
    assert len(rows) == 6
    assert all(r["difference"] == 0 for r in rows)
    assert [r["endpoint"] for r in rows if r["sampling_confounded"]] == ["FAS"]


@pytest.mark.parametrize("mutation", ["status", "admission", "participant", "missing", "nan", "mean", "axes", "semantics"])
def test_reject_invalid_comparison(mutation):
    old, new = inputs()
    a = new["assessment"]
    if mutation == "status":
        new["status"] = "running"
    elif mutation == "admission":
        new["accuracy_admitted"] = False
    elif mutation == "participant":
        a["participant"] = "other"
    elif mutation == "missing":
        del a["endpoints"]["FAS"]
    elif mutation == "nan":
        a["endpoints"]["GO"]["score"] = float("nan")
    elif mutation == "mean":
        a["secondary_six_metric_mean"] = 0.7
    elif mutation == "axes":
        a["endpoints"]["GO"]["axes"] = {"x_axis": "wrong"}
    else:
        a["endpoints"]["GO"]["score_semantics"] = "F1"
    with pytest.raises(ValueError):
        compare(old, new)
