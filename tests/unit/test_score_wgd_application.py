from copy import deepcopy

import pytest

from benchmark_tools.score_wgd_application import score_pairs, summarize


def fixture():
    owners = {"a": "Scerevisiae", "b": "Scerevisiae", "c": "Smikatae",
              "d": "Suvarum", "foreign": "Suvarum", "unknown": "Smikatae"}
    reference = {"p": ["a", "b", "c", "d"], "q": ["foreign"]}
    pairs = [{"orf_pair": ["a", "b"], "split_eligible": True, "reference_eligible": True,
              "reference_pillar": "p", "available_pillar_members": reference["p"],
              "experimental_class": "Sparse"}]
    return pairs, reference, owners


def score(groups):
    pairs, reference, owners = fixture()
    return score_pairs(pairs, groups, reference, owners)[0]


def test_supported_split_separates_unknown_from_foreign():
    row = score({"one": ["a", "c", "foreign", "unknown"], "two": ["b", "d"]})
    assert row["separation_rate"] == row["supported_separation_rate"] == 1
    assert row["mean_non_scer_coverage"] == 1
    assert row["foreign_pillar_members"] == ["foreign"]
    assert row["unmapped_members"] == ["unknown"]
    assert row["pillar_native_group_count"] == 2


def test_singleton_splits_not_supported():
    row = score({"one": ["a"], "two": ["b"], "three": ["c", "d"]})
    assert row["separation_rate"] == 1 and row["supported_separation_rate"] == 0
    assert row["mean_non_scer_coverage"] == 0 and row["pillar_native_group_count"] == 3


def test_merged_group_has_coverage_but_no_separation():
    row = score({"one": ["a", "b", "c", "d"]})
    assert row["assignment_state"] == "merged"
    assert row["mean_non_scer_coverage"] == 1
    assert row["separation_rate"] == row["supported_separation_rate"] == 0
    assert row["union_size"] == 4


def test_one_missing_anchor_is_not_synthetic_singleton():
    row = score({"one": ["a", "c"]})
    assert row["assignment_state"] == "incomplete_assignment"
    assert row["anchor_groups"] == ["one", None]
    assert row["mean_non_scer_coverage"] == 0.5
    assert row["supported_separation_rate"] == row["separation_rate"] == 0


def test_both_missing_anchors_give_zero_not_missing_endpoint():
    row = score({"one": ["c", "d"]})
    assert row["mean_non_scer_coverage"] == 0
    assert row["separation_rate"] == 0


def test_each_anchor_requires_its_own_homolog_support():
    row = score({"one": ["a", "c", "d"], "two": ["b"]})
    assert row["homolog_support_by_anchor"] == [2, 0]
    assert row["supported_separation_rate"] == 0


@pytest.mark.parametrize("groups", [{"x": ["a", "a"]}, {"x": ["a"], "y": ["a"]}, {"x": ["alien"]}, {"x": []}])
def test_invalid_native_membership_rejected(groups):
    with pytest.raises(ValueError):
        score(groups)


def test_exclusions_remain_rows_with_unavailable_reference_endpoints():
    pairs, reference, owners = fixture()
    pairs[0]["reference_eligible"] = False
    row = score_pairs(pairs, {"x": ["a"], "y": ["b"]}, reference, owners)[0]
    assert row["separation_rate"] == 1 and row["supported_separation_rate"] is None
    pairs[0]["split_eligible"] = False
    rows = score_pairs(pairs, {}, reference, owners)
    assert rows[0]["separation_rate"] is None
    assert summarize(rows)["pairs"] == 1
    assert summarize(rows)["assignment_states"] == {"input_excluded": 1}


def test_reference_mismatch_and_duplicate_pairs_rejected():
    pairs, reference, owners = fixture()
    with pytest.raises(ValueError, match="duplicate"):
        score_pairs(pairs + deepcopy(pairs), {}, reference, owners)
    pairs[0]["available_pillar_members"] = ["a", "b"]
    with pytest.raises(ValueError, match="differs"):
        score_pairs(pairs, {}, reference, owners)


def test_summary_preserves_sparse_class_and_fixed_denominator():
    row = score({"x": ["a", "c"]})
    assert row["experimental_class"] == "Sparse"
    summary = summarize([row])
    assert summary["endpoints"]["separation_rate"] == {"pairs": 1, "mean": 0}
    assert summary["endpoints"]["mean_non_scer_coverage"]["mean"] == 0.5
