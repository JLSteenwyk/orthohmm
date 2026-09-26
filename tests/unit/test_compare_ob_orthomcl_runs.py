import pytest

from benchmark_tools.compare_ob_orthomcl_runs import differences, partition


def test_ignores_labels_and_order_not_membership():
    assert differences([("a", ["x", "y"]), ("b", ["z"])],
                       [("new", ["z"]), ("other", ["y", "x"])])["identical_groups"] == 2


def test_detects_split_beyond_missing_gene():
    result = differences([("a", ["a", "b", "c", "d"])], [("b", ["a", "b"]), ("c", ["c"])])
    assert result["april_only_genes"] == ["d"]
    assert result["july_only_groups"] == [["a", "b"], ["c"]]
    assert not result["equal_after_removing_april_only_genes"]


def test_detects_pure_gene_removal():
    assert differences([("a", ["x", "y"])], [("b", ["x"])])["equal_after_removing_april_only_genes"]


@pytest.mark.parametrize("rows", [[("a", [])], [("a", ["x", "x"])],
                                   [("a", ["x"]), ("b", ["x"])],
                                   [("a", ["x"]), ("a", ["y"])]])
def test_rejects_invalid_partition(rows):
    with pytest.raises(ValueError):
        partition(rows)
