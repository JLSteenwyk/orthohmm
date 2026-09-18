from copy import deepcopy

import pytest

from benchmark_tools.summarize_ob_search_stage_trace import STAGES, summarize


@pytest.fixture
def fixture():
    row = {"refog": "family", "left": "a", "right": "b", "same_species": "False",
           "forward_normalized_hit": "2", "reverse_normalized_hit": "NA",
           **{stage: "True" for stage in STAGES}}
    row["root_hogs"] = "False"
    family = {"possible_pairs": 1, "search": {"both_direction_pairs": 0, "either_direction_pairs": 1},
              "stages": {stage: {"groups": [{"family_genes": ["a", "b"]}],
                                  "within_family_pairs": int(stage != "root_hogs")} for stage in STAGES}}
    return [row], {"family": family}, [("candidates", "root_hogs", "reconciliation")]


def test_hit_retained_but_grouping_lost(fixture):
    counts = summarize(*fixture)["aggregate_membership_counts"]
    assert counts["all/hits1/candidates_to_root_hogs/lost"] == 1
    assert counts["cross_species/hits1/pairs"] == 1


@pytest.mark.parametrize("field,value", [("same_species", "yes"), ("forward_normalized_hit", "nan"),
    ("forward_normalized_hit", "0"), ("root_hogs", "NA"), ("left", "b"), ("right", "unknown")])
def test_invalid_trace_rejected(fixture, field, value):
    fixture[0][0][field] = value
    with pytest.raises(ValueError):
        summarize(*fixture)


def test_missing_pair_rejected(fixture):
    fixture[0].clear()
    with pytest.raises(ValueError, match="Incomplete"):
        summarize(*fixture)


def test_duplicate_pair_rejected(fixture):
    fixture[0].append(deepcopy(fixture[0][0]))
    with pytest.raises(ValueError, match="duplicate"):
        summarize(*fixture)


def test_wrong_marginal_rejected(fixture):
    fixture[1]["family"]["stages"]["root_hogs"]["within_family_pairs"] = 1
    with pytest.raises(ValueError, match="marginals"):
        summarize(*fixture)


def test_no_direct_hit_but_grouped_is_retained(fixture):
    fixture[0][0]["forward_normalized_hit"] = "NA"
    fixture[1]["family"]["search"]["either_direction_pairs"] = 0
    counts = summarize(*fixture)["aggregate_membership_counts"]
    assert counts["all/hits0/candidates/together"] == 1
