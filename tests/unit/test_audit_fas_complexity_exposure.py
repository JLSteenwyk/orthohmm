import pytest

from benchmark_tools.audit_fas_complexity_exposure import summarize, verify_lookup


def test_both_flagged_counted_once():
    result = summarize([(("A", "B"), 0.5), (("A", "C"), 0.7), (("C", "D"), 0.2)],
                       {"A", "B"}, {"A", "B", "C", "D"})
    assert result["touches_flagged_protein"]["pairs"] == 2
    assert result["touches_flagged_protein"]["mean"] == 0.6
    assert result["neither_flagged"] == dict(pairs=1, mean=0.2)


def test_empty_category_is_not_zero_mean():
    result = summarize([(("A", "B"), 0.5)], set(), {"A", "B"})
    assert result["touches_flagged_protein"] == dict(pairs=0, mean=None)


def test_missing_annotation_endpoint_fails():
    with pytest.raises(ValueError, match="absent"):
        summarize([(("A", "B"), 0.5)], set(), {"A"})


@pytest.mark.parametrize("score", [float("nan"), float("inf"), -0.1, 1.1])
def test_invalid_score_fails(score):
    with pytest.raises(ValueError, match="Invalid"):
        summarize([(("A", "B"), score)], set(), {"A", "B"})


def test_lookup_checks_reversed_pair_and_numeric_mean():
    result = verify_lookup({("A", "B"): 0.5, ("C", "D"): None},
                           [("B_A", ["0.2", "0.8"]), ("X_Y", ["NA", "NA"])])
    assert result == dict(entries_scanned=2, unique_saved_precomputed_pairs=1, unique_saved_new_pairs=1,
                         relevant_canonical_overwrites=0)


@pytest.mark.parametrize("entries", [[], [("A_B", [0.1, 0.2])], [("A_B", ["NA", "NA"])],
    [("A_B", [0.5, 0.5]), ("A_B", [0.5, 0.5])], [("A_B", [float("nan"), 0.5])]])
def test_bad_lookup_rejected(entries):
    with pytest.raises(ValueError):
        verify_lookup({("A", "B"): 0.5}, entries)


def test_new_pair_cannot_have_valid_precomputed_score():
    with pytest.raises(ValueError, match="Logged new"):
        verify_lookup({("A", "B"): None}, [("A_B", [0.5, 0.5])])


def test_last_valid_canonical_orientation_wins_like_native_loader():
    result = verify_lookup({("A", "B"): 0.5}, [("A_B", [0.1, 0.3]), ("B_A", [0.5, 0.5])])
    assert result["relevant_canonical_overwrites"] == 1


def test_later_invalid_orientation_does_not_erase_valid_score():
    result = verify_lookup({("A", "B"): 0.5}, [("A_B", [0.5, 0.5]), ("B_A", ["NA", "NA"])])
    assert result["unique_saved_precomputed_pairs"] == 1
