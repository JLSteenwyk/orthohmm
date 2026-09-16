import pytest

from benchmark_tools.validate_factorial_partition import validate_partition


def fixture(tmp_path):
    candidate = tmp_path / "candidate.txt"
    root = tmp_path / "root.tsv"
    candidate.write_text("a b\nc\n")
    root.write_text("root_hog\tsource_family\tgenes\nRootHOG0000000\tFamily0000000\ta,b\nRootHOG0000001\tFamily0000001\tc\n")
    counts = {"candidate_families": 2, "root_hogs": 2, "reconciled_families": 1, "bypassed_families": 1}
    return candidate, root, counts


def test_preserves_singletons_and_all_genes(tmp_path):
    candidate, root, counts = fixture(tmp_path)
    groups, summary = validate_partition(candidate, root, {"a", "b", "c"}, counts)
    assert groups == [{"a", "b"}, {"c"}]
    assert summary["all_candidate_genes_preserved"]
    assert summary["cross_source_merges"] == 0


def test_valid_split_preserves_both_descendant_groups(tmp_path):
    candidate, root, counts = fixture(tmp_path)
    root.write_text("root_hog\tsource_family\tgenes\nRootHOG0000000\tFamily0000000\ta\n"
                    "RootHOG0000001\tFamily0000000\tb\nRootHOG0000002\tFamily0000001\tc\n")
    counts["root_hogs"] = 3
    groups, summary = validate_partition(candidate, root, {"a", "b", "c"}, counts)
    assert groups == [{"a"}, {"b"}, {"c"}]
    assert summary["split_source_families"] == 1


@pytest.mark.parametrize("change", ["duplicate", "missing", "unknown", "cross_source", "negative_source",
                                    "duplicate_hog", "candidate_duplicate", "count", "completion", "header"])
def test_rejects_invalid_partition(tmp_path, change):
    candidate, root, counts = fixture(tmp_path)
    text = root.read_text()
    if change == "duplicate":
        text = text.replace("a,b", "a,a,b")
    elif change == "missing":
        text = text.replace("a,b", "a")
    elif change == "unknown":
        text = text.replace("a,b", "a,b,x")
    elif change == "cross_source":
        text = text.replace("Family0000001", "Family0000000")
    elif change == "negative_source":
        text = text.replace("Family0000001", "Family-1")
    elif change == "duplicate_hog":
        text = text.replace("RootHOG0000001", "RootHOG0000000")
    elif change == "candidate_duplicate":
        candidate.write_text("a b\nb c\n")
    elif change == "count":
        counts["root_hogs"] = 3
    elif change == "completion":
        counts["bypassed_families"] = 0
    else:
        text = text.replace("source_family", "wrong")
    root.write_text(text)
    with pytest.raises(ValueError):
        validate_partition(candidate, root, {"a", "b", "c"}, counts)
