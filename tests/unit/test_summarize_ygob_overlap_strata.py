from copy import deepcopy

import pytest

from benchmark_tools.score_ygob_groups import score_groups
from benchmark_tools.summarize_ygob_overlap_strata import METHODS, summarize


def fixture():
    score = score_groups({"x": ["a", "b", "c"], "y": ["d"]},
                         {"A": ["a", "b"], "B": ["c", "d"]}, ["a", "b", "c", "d"])
    return {m: deepcopy(score) for m in METHODS}


def test_cross_stratum_false_positives_retained():
    scores = fixture()
    result = summarize(scores, ["A"])
    a = result["screen_positive"][METHODS[0]]
    b = result["screen_negative"][METHODS[0]]
    assert a["counts"] == dict(tp=1, fp=1, fn=0)
    assert b["counts"] == dict(tp=0, fp=1, fn=1)
    assert a["metrics"]["f1"] == pytest.approx(2/3)
    assert a["reference_gene_coverage"] == b["reference_gene_coverage"] == 1
    assert a["counts"]["fp"] + b["counts"]["fp"] == scores[METHODS[0]]["counts"]["fp"]


def test_half_counts_and_undefined_ratios():
    score = score_groups({"x": ["a", "b"]}, {"A": ["a"], "B": ["b"]}, ["a", "b"])
    result = summarize({m: deepcopy(score) for m in METHODS}, ["A"])
    for methods in result.values():
        row = methods[METHODS[0]]
        assert row["counts"]["fp"] == .5
        assert row["zero_truth_pair_pillars"] == 1
        assert row["defined"] == dict(f1=True, precision=True, recall=False)
        assert row["difference_vs_orthofinder_percentage_points"]["recall"] is None


def test_empty_stratum_explicitly_undefined():
    result = summarize(fixture(), [])
    row = result["screen_positive"][METHODS[0]]
    assert row["reference_groups"] == 0
    assert not any(row["defined"].values())
    assert row["metrics"] == dict(f1=0, precision=0, recall=0)


@pytest.mark.parametrize("labels", [["A", "A"], ["C"], [1], "A"])
def test_invalid_labels(labels):
    with pytest.raises(ValueError):
        summarize(fixture(), labels)


@pytest.mark.parametrize("key,value", [
    ("tp", float("nan")), ("fp", -.5), ("fp", .25), ("tp", .5),
    ("genes", True), ("covered_genes", 3), ("covered_genes", 1),
    ("exact", True), ("exact", 1), ("fn", 100),
])
def test_invalid_pillar_statistics(key, value):
    scores = fixture()
    scores[METHODS[0]]["records"][0][key] = value
    with pytest.raises(ValueError):
        summarize(scores, ["A"])


@pytest.mark.parametrize("change", ["duplicate", "signature", "totals", "metrics", "missing_method"])
def test_rejects_inconsistent_sources(change):
    scores = fixture()
    score = scores[METHODS[0]]
    if change == "duplicate":
        score["records"].append(deepcopy(score["records"][0]))
    elif change == "signature":
        scores[METHODS[1]]["records"][1]["pillar"] = "C"
    elif change == "totals":
        score["counts"]["fp"] += 1
    elif change == "metrics":
        score["metrics"]["f1"] = 1
    else:
        scores.pop(METHODS[-1])
    with pytest.raises(ValueError):
        summarize(scores, ["A"])


def test_order_invariance_and_no_input_mutation():
    scores = fixture()
    saved = deepcopy(scores)
    expected = summarize(scores, ["A"])
    assert scores == saved
    for score in scores.values():
        score["records"].reverse()
    assert summarize(scores, ["A"]) == expected
