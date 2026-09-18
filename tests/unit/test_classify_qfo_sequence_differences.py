import pytest

from benchmark_tools.classify_qfo_sequence_differences import classify


def test_noncanonical_to_x_only():
    result = classify("ACUOBZJ", "ACXXXXX")
    assert result["category"] == "noncanonical_to_X_only"
    assert result["changed_positions"] == 5
    assert sum(r["count"] for r in result["residue_changes"]) == 5


@pytest.mark.parametrize("original,native", [("ACU", "ACG"), ("AAA", "AAX"), ("ACD", "ACE"), ("AXA", "AAA")])
def test_other_changes_not_hidden(original, native):
    assert classify(original, native)["category"] == "other_same_length"


def test_unequal_lengths_not_zipped_or_aligned():
    result = classify("ACD", "ACDE")
    assert result["category"] == "different_length_unaligned"
    assert result["length_delta_native_minus_input"] == 1
    assert result["residue_changes"] is None


def test_identical_rejected():
    with pytest.raises(ValueError):
        classify("ACD", "ACD")
