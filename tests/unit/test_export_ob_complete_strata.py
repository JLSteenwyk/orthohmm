from copy import deepcopy

import pytest

from benchmark_tools.export_ob_complete_strata import summarize, CATEGORIES


def fixture():
    records = [dict(refog="a", genes=2, true_positive=1, false_positive=0, false_negative=0),
               dict(refog="b", genes=3, true_positive=2, false_positive=2, false_negative=1)]
    scores = {name: dict(refog_records=deepcopy(records)) for name in ("first", "second")}
    strata = {d + ":" + c: ["a", "b"] if i == 0 else [] for d, cs in CATEGORIES.items() for i, c in enumerate(cs)}
    return scores, strata


def test_weighted_recomputation_and_empty_bins():
    scores, strata = fixture()
    result = summarize(scores, strata)
    assert len(result) == 28
    for row in result:
        if row["family_count"]:
            assert row["weighted_counts"] == dict(tp=2., fp=1., fn=.5)
            assert row["metrics_percent"]["f_score"] == pytest.approx(100 * 4 / 5.5)
            assert row["status"] == "descriptive"
        else:
            assert row["metrics_percent"] is None
            assert row["weighted_counts"] is None
            assert row["status"] == "empty_nonestimable"


@pytest.mark.parametrize("change", ["missing_bin", "duplicate_family", "missing_family", "different_names", "different_sizes"])
def test_reject_inconsistent_inputs(change):
    scores, strata = fixture()
    label = next(k for k, v in strata.items() if v)
    if change == "missing_bin":
        strata.pop(label)
    elif change == "duplicate_family":
        strata[label].append("a")
    elif change == "missing_family":
        strata[label].pop()
    elif change == "different_names":
        scores["second"]["refog_records"][0]["refog"] = "c"
    else:
        scores["second"]["refog_records"][0].update(genes=3, false_negative=2)
    with pytest.raises(ValueError):
        summarize(scores, strata)
