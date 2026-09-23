import pytest

from benchmark_tools.inspect_swisstree_line_models import inspect_lines


def test_two_unterminated_candidates_not_concatenated():
    result = inspect_lines("((A,B)[&&NHX:D=Y],C)\n((B,A)[&&NHX:D=Y],C)\n")
    assert result["candidate_count"] == 2
    assert result["identical_leaf_sets"] is result["identical_stored_clades"] is result["identical_D_Y_clades"] is True
    assert result["chosen_candidate"] is None
    assert all(r["appended_terminator_for_parse"] for r in result["candidates"])
    assert result["candidates"][0]["explicit_D_Y_clades"] == [["A", "B"]]


def test_different_candidates_and_missing_labels_not_equated():
    result = inspect_lines("((A,B)[&&NHX:D=Y],C);\n((A,C),B);\n")
    assert result["identical_leaf_sets"] is True
    assert result["identical_stored_clades"] is result["identical_D_Y_clades"] is False
    assert not any(r["appended_terminator_for_parse"] for r in result["candidates"])


def test_empty_rejected():
    with pytest.raises(ValueError):
        inspect_lines("\n")
