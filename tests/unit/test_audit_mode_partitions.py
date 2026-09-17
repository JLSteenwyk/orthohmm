import pytest

from benchmark_tools.audit_mode_partitions import canonical


def test_group_labels_order_and_missing_genes_are_explicit():
    assert canonical({"x": ["b", "a"]}, {"a", "b", "c"}) == ({("a", "b")}, {"c"})
    assert canonical({"x": ["b", "a"]}, {"a", "b"}) == canonical({"y": ["a", "b"]}, {"a", "b"})


@pytest.mark.parametrize("groups", [{"x": ["a", "a"]}, {"x": ["a"], "y": ["a"]}, {"x": ["unknown"]}])
def test_invalid_memberships_rejected(groups):
    with pytest.raises(ValueError):
        canonical(groups, {"a", "b"})
