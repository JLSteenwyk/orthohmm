import pytest

from benchmark_tools.trace_wgd_cases import loss_stage, partition_index


@pytest.mark.parametrize("groups", [{"x": []}, {"x": ["a", "a"]}, {"x": ["a"], "y": ["a"]}])
def test_invalid_partition_rejected(groups):
    with pytest.raises(ValueError):
        partition_index(groups)


def test_partition_index_preserves_group_ids():
    assert partition_index({"x": ["a", "b"], "y": ["c"]}) == {"a": "x", "b": "x", "c": "y"}


@pytest.mark.parametrize("candidate,pre,final,expected", [
    ("c", "p", "f", "retained_in_anchor_group"),
    ("outside", "elsewhere", "elsewhere", "outside_anchor_candidates"),
    ("c", "elsewhere", "elsewhere", "root_lineage_split"),
    ("c", "p", "elsewhere", "satellite_constraint_split"),
])
def test_loss_categories_follow_execution_order(candidate, pre, final, expected):
    c = {"a": "c", "b": "d", "h": candidate}
    p = {"a": "p", "b": "q", "h": pre}
    f = {"a": "f", "b": "g", "h": final}
    assert loss_stage("h", ("a", "b"), c, p, f) == expected


def test_homolog_in_second_anchor_group_is_retained():
    assert loss_stage("h", ("a", "b"), {"a": 0, "b": 1, "h": 1},
                      {"a": 2, "b": 3, "h": 3}, {"a": 4, "b": 5, "h": 5}) == "retained_in_anchor_group"
