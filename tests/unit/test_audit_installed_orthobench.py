import pytest

from benchmark_tools.audit_installed_orthobench import compare_partitions, read_root_hogs


def test_strict_partition_and_label_invariance(tmp_path):
    path = tmp_path / "groups.tsv"
    path.write_text("root_hog\tsource_family\tgenes\nH1\tF1\ta,b\nH2\tF2\tc\n")
    groups = read_root_hogs(path, {"a", "b", "c"})
    assert compare_partitions(groups, groups[::-1])["label_invariant_equal"]
    split = [frozenset("a"), frozenset("b"), frozenset("c")]
    result = compare_partitions(groups, split)
    assert result["genes_in_changed_groups"] == 2
    assert result["identical_groups"] == 1


@pytest.mark.parametrize("body", [
    "H1\tF1\ta,a\nH2\tF2\tb\n", "H1\tF1\ta,b\nH2\tF2\ta\n",
    "H1\tF1\ta\n", "H1\tF1\ta,b,x\n", "H1\tF1\ta\nH1\tF1\tb\n",
    "H1\t\ta,b\n", "H1\tF1\ta,b,\n",
])
def test_invalid_memberships_rejected(tmp_path, body):
    path = tmp_path / "groups.tsv"
    path.write_text("root_hog\tsource_family\tgenes\n" + body)
    with pytest.raises(ValueError):
        read_root_hogs(path, {"a", "b"})
