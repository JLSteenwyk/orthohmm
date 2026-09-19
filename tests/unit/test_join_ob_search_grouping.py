import itertools

import pytest

from benchmark_tools.join_ob_search_grouping import category, join, DECISIONS


@pytest.mark.parametrize("left,right", list(itertools.product(sorted(DECISIONS), repeat=2)))
def test_category_symmetric(left, right):
    assert category(left, right) == category(right, left)
    assert (category(left, right) == "at_least_one_accepted_direction") == ("accepted" in (left, right))


def row():
    return dict(refog="f", left="a", right="b", forward_normalized_hit="NA",
                reverse_normalized_hit="NA", root_hogs="True", same_species="False")


def test_missing_direct_hits_can_be_grouped_and_overlaps_retained():
    first = row()
    second = {**first, "refog": "overlap"}
    values = {("a", "b"): "not_selected_by_prefilter", ("b", "a"): "scored_not_significant"}
    rows, counts, families = join([first, second], values)
    assert len(rows) == 2 and len(families) == 2
    assert counts == {"one_prefilter_excluded_one_scored_not_significant:grouped:cross_species": 2}


def test_reject_historical_mismatch():
    with pytest.raises(ValueError, match="disagrees"):
        join([row()], {("a", "b"): "accepted", ("b", "a"): "not_selected_by_prefilter"})


def test_reject_duplicate_pair():
    with pytest.raises(ValueError, match="Duplicate"):
        join([row(), row()], {("a", "b"): "scored_not_significant", ("b", "a"): "scored_not_significant"})


def test_reject_unknown_decision():
    with pytest.raises(ValueError):
        category("missing", "accepted")
