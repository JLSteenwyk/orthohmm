import json

import numpy as np
import pytest

from orthohmm.search.profile import ProfileHMM
from orthohmm.search.profile_spans import ProfileAlignmentSpan, profile_alignment_span
from orthohmm.search.viterbi import _viterbi_score_one, batch_viterbi_c, is_c_available


def profile(symbols, match=5, mismatch=-20, gap_open=-12, gap_extend=-3, gap_close=-1, insert=-1):
    emissions = np.full((len(symbols), 20), mismatch, dtype=np.int8)
    for position, symbol in enumerate(symbols):
        emissions[position, symbol] = match
    return ProfileHMM(len(symbols), emissions, np.full(20, insert, dtype=np.int8),
                      np.array([0, gap_open, gap_open, gap_close, gap_extend, gap_close, gap_extend], dtype=np.int32))


def target(symbols):
    return np.array(symbols, dtype=np.uint8)


def test_exact_full_length_coordinates_and_coverage():
    result = profile_alignment_span(profile([0, 1, 2]), target([0, 1, 2]), expected_score=15)
    assert result == ProfileAlignmentSpan(15, 0, 3, 0, 3, 3, 3, 3)
    assert result.passes_guard and result.category == "both_sides_covered"
    record = result.as_record()
    assert record["aligned_model_coverage"] == record["aligned_target_coverage"] == 1.
    json.dumps(record, allow_nan=False)


def test_local_restart_discards_unrelated_flanks():
    result = profile_alignment_span(profile([1, 1, 1, 0, 0]), target([2] * 6 + [0, 0]))
    assert result == ProfileAlignmentSpan(10, 3, 5, 6, 8, 2, 5, 8)
    assert not result.passes_guard and result.category == "partial_shared_region"


def test_fragment_and_fusion_compatible_asymmetry_is_not_automatically_rejected():
    result = profile_alignment_span(profile([0, 1, 2, 3, 4]), target([1, 2, 3]))
    assert result == ProfileAlignmentSpan(15, 1, 4, 0, 3, 3, 5, 3)
    assert result.passes_guard and result.category == "asymmetric_coverage"
    assert "not orthology truth" in result.as_record()["interpretation"]


def test_insertions_expand_bounding_span_without_counting_as_alignment():
    result = profile_alignment_span(profile([0, 1], gap_open=-1, gap_extend=0, gap_close=0, insert=0),
                                    target([0, 2, 2, 2, 2, 1]))
    assert result == ProfileAlignmentSpan(9, 0, 2, 0, 6, 2, 2, 6)
    record = result.as_record()
    assert record["bounding_target_coverage"] == 1
    assert record["aligned_target_coverage"] == pytest.approx(1 / 3)
    assert result.passes_guard and result.category == "asymmetric_coverage"


def test_deletions_expand_model_span_without_counting_as_alignment():
    result = profile_alignment_span(profile([0, 2, 2, 2, 2, 1], gap_open=-1, gap_extend=0, gap_close=0),
                                    target([0, 1]))
    assert result == ProfileAlignmentSpan(9, 0, 6, 0, 2, 2, 6, 2)
    record = result.as_record()
    assert record["bounding_model_coverage"] == 1
    assert record["aligned_model_coverage"] == pytest.approx(1 / 3)


def test_unknown_match_states_are_not_counted_as_known_alignment_evidence():
    result = profile_alignment_span(profile([0, 0, 0]), target([0, 20, 0]))
    assert result == ProfileAlignmentSpan(10, 0, 3, 0, 3, 2, 3, 3)
    assert result.as_record()["bounding_model_coverage"] == 1
    assert not result.passes_guard and result.category == "partial_shared_region"


def test_repeated_equal_optima_use_first_terminal_match_not_a_selected_favorable_span():
    result = profile_alignment_span(profile([0]), target([0, 0, 0]))
    assert result == ProfileAlignmentSpan(5, 0, 1, 0, 1, 1, 1, 3)


@pytest.mark.parametrize("symbols,sequence", [([], []), ([], [0]), ([0], []), ([0], [20]), ([0], [255])])
def test_empty_or_unknown_only_inputs_have_zero_spans(symbols, sequence):
    result = profile_alignment_span(profile(symbols), target(sequence), expected_score=0)
    assert result.score == result.model_start == result.model_end == result.target_start == result.target_end == 0
    assert result.matched_known_residues == 0 and not result.passes_guard
    assert result.category == "no_positive_alignment"
    json.dumps(result.as_record(), allow_nan=False)


def test_score_mismatch_fails_instead_of_falling_back_to_another_alignment():
    with pytest.raises(ValueError, match="score mismatch"):
        profile_alignment_span(profile([0]), target([0]), expected_score=4)


@pytest.mark.parametrize("expected", [True, np.nan, np.inf, -1, 5.5, "5"])
def test_invalid_expected_scores_rejected(expected):
    with pytest.raises(ValueError, match="Expected score"):
        profile_alignment_span(profile([0]), target([0]), expected_score=expected)


@pytest.mark.parametrize("band", [True, 1.5, None, "64", 2 ** 32])
def test_invalid_band_types_rejected(band):
    with pytest.raises(ValueError, match="Band width"):
        profile_alignment_span(profile([0]), target([0]), band_width=band)


@pytest.mark.parametrize("sequence", [[[0]], [-1], [256], [0.5]])
def test_invalid_target_shape_or_values_rejected(sequence):
    with pytest.raises(ValueError, match="Target|target"):
        profile_alignment_span(profile([0]), sequence)


@pytest.mark.parametrize("field,value", [
    ("length", True), ("length", -1), ("length", 1.5),
    ("match_emissions", np.zeros((1, 21), dtype=np.int8)),
    ("match_emissions", np.full((1, 20), 128, dtype=np.int16)),
    ("insert_emissions", np.zeros(19, dtype=np.int8)),
    ("transitions", np.ones(7, dtype=np.int32)),
])
def test_invalid_profile_arrays_and_lengths_rejected(field, value):
    model = profile([0])
    setattr(model, field, value)
    with pytest.raises(ValueError, match="Invalid"):
        profile_alignment_span(model, target([0]))


@pytest.mark.parametrize("band", [0, 1, 8, 64])
def test_extended_recurrence_score_matches_original_numba_and_installed_c(band):
    if not is_c_available():
        pytest.skip("Requires original native HMM library for exact-score comparison")
    rng = np.random.default_rng(18102026)
    lengths = [0, 1, 7, 49, 50, 51, 97]
    target_lengths = [0, 1, 3, 49, 50, 51, 89, 201]
    for length in lengths:
        emissions = rng.integers(-10, 12, size=(length, 20), dtype=np.int8)
        model = ProfileHMM(length, emissions, np.full(20, -1, dtype=np.int8),
                           np.array([0, -12, -12, -1, -3, -1, -3], dtype=np.int32))
        for size in target_lengths:
            sequence = rng.integers(0, 21, size=size, dtype=np.uint8)
            before = [array.tobytes() for array in (emissions, model.insert_emissions, model.transitions, sequence)]
            reference = _viterbi_score_one(emissions, model.insert_emissions, model.transitions, sequence, size, length, band)
            c_score = batch_viterbi_c(emissions, model.insert_emissions, model.transitions,
                np.array([0], dtype=np.int64), np.array([length], dtype=np.int32), sequence,
                np.array([0], dtype=np.int64), np.array([size], dtype=np.int32),
                np.array([[0, 0]], dtype=np.int32), band, n_threads=1)[0]
            result = profile_alignment_span(model, sequence, band_width=band, expected_score=c_score)
            assert result.score == reference
            assert [array.tobytes() for array in (emissions, model.insert_emissions, model.transitions, sequence)] == before


def test_guard_boundary_is_exact_integer_arithmetic():
    assert ProfileAlignmentSpan(1, 0, 5, 0, 5, 4, 5, 5).passes_guard
    assert not ProfileAlignmentSpan(1, 0, 5, 0, 5, 3, 5, 5).passes_guard
