import pytest

from benchmark_tools.audit_dgx_orthohmm_smoke import compare, pair_set, partition


def test_partition_ignores_group_labels(tmp_path):
    path = tmp_path / "groups"
    path.write_text("OG2: b a\nOG1: c\n")
    assert partition(path, "named_groups", {"a": "s1", "b": "s2", "c": "s3"}) == {("a", "b"), ("c",)}


@pytest.mark.parametrize("text", ["G: a b\n", "G: a b c\nH: a\n"])
def test_invalid_partition(tmp_path, text):
    path = tmp_path / "groups"
    path.write_text(text)
    with pytest.raises(ValueError):
        partition(path, "named_groups", {"a": "s1", "b": "s2", "c": "s3"})


def test_pair_duplicates_rejected(tmp_path):
    path = tmp_path / "pairs"
    path.write_text("gene_a\tspecies_a\tgene_b\tspecies_b\na\ts1\tb\ts2\na\ts1\tb\ts2\n")
    with pytest.raises(ValueError, match="unique"):
        pair_set(path, {"a": "s1", "b": "s2"})


def test_differences_not_hidden():
    assert compare({1, 2}, {2, 3}) == {"x86_count": 2, "arm_count": 2, "equal": False,
                                    "x86_only": 1, "arm_only": 1}
