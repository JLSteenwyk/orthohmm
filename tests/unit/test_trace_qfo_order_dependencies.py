import pytest

from benchmark_tools.trace_qfo_order_dependencies import family_records, trace_changes


def row(source, target):
    return dict(source_genes=source, target_genes=target)


def test_numbering_and_membership_are_separate():
    left = [("F0", ("a",)), ("F1", ("b", "c")), ("F2", ("d",))]
    right = [("F0", ("d",)), ("F1", ("a", "b")), ("F2", ("c",))]
    result = trace_changes(left, right, [], [])
    assert result["shared_families"] == 1
    assert result["shared_families_with_changed_ids"] == 1
    assert len(result["changed_families"]["left"]) == 2


def test_multiset_distinguishes_reordering_repetition_and_direction():
    a, b = row(["a"], ["b"]), row(["b"], ["a"])
    result = trace_changes([], [], [a, b], [b, a])
    assert result["constraint_multiset_difference"] == dict(left_only=[], right_only=[])
    result = trace_changes([], [], [a, a], [b])
    assert result["constraint_multiset_difference"]["left_only"][0]["count"] == 2
    assert result["constraint_multiset_difference"]["right_only"][0]["source_genes"] == ["b"]


def test_family_parser_matches_native_sort_and_ignores_blank_lines(tmp_path):
    path = tmp_path / "partition"
    path.write_text("b a\n\nc\n")
    assert family_records(path) == [("Family0000000", ("a", "b")), ("Family0000001", ("c",))]


def test_family_parser_rejects_duplicates(tmp_path):
    path = tmp_path / "partition"
    path.write_text("a b\na\n")
    with pytest.raises(ValueError):
        family_records(path)
