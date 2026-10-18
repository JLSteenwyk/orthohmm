"""Selective span evidence from OrthoHMM's existing integer Viterbi model.

Coordinates describe one deterministic best path, not an alignment posterior
or a biological domain annotation. Insertions/deletions widen spans but do
not increase the count of aligned, known amino-acid match states.
"""

from dataclasses import asdict, dataclass
from numbers import Integral, Real

import numpy as np
from numba import njit

from .viterbi import NEG_INF


PROFILE_ALIGNMENT_POLICIES = ("off", "strict_span_guard")


def validate_profile_alignment_policy(policy):
    if not isinstance(policy, str) or policy not in PROFILE_ALIGNMENT_POLICIES:
        raise ValueError("Unknown profile alignment policy; expected off or strict_span_guard")
    return policy


@njit(cache=True, nogil=True)
def _span_one(match_emit, insert_emit, transitions, target, band_width):
    length, target_length = len(match_emit), len(target)
    result = np.zeros(6, dtype=np.int64)
    if length == 0 or target_length == 0:
        return result
    width = max(length, target_length)
    if band_width > 0 and width > 50:
        width = band_width
    # State records contain score, model start, target start and known matches.
    previous = np.zeros((3, target_length + 1, 4), dtype=np.int32)
    current = np.zeros_like(previous)
    previous[:, :, 0] = NEG_INF
    current[:, :, 0] = NEG_INF
    for i in range(1, length + 1):
        center = (i * target_length) // length
        lo, hi = max(1, center - width), min(target_length, center + width)
        for j in range(lo - 1, min(target_length, hi + 1) + 1):
            for state in range(3):
                current[state, j, 0] = NEG_INF
        for j in range(lo, hi + 1):
            amino_acid = target[j - 1]
            known = int(amino_acid < 20)
            emission = int(match_emit[i - 1, amino_acid]) if known else 0
            insertion = int(insert_emit[amino_acid]) if known else -1

            # Fresh start wins a zero tie, followed by M, I, then D.
            best, predecessor = 0, -1
            for state, transition in ((0, 0), (1, 3), (2, 5)):
                value = previous[state, j - 1, 0] + transitions[transition]
                if value > best:
                    best, predecessor = value, state
            current[0, j, 0] = best + emission
            if predecessor < 0:
                current[0, j, 1] = i - 1
                current[0, j, 2] = j - 1
                current[0, j, 3] = known
            else:
                for field in range(1, 4):
                    current[0, j, field] = previous[predecessor, j - 1, field]
                current[0, j, 3] += known

            best, predecessor = NEG_INF, 0
            for state, transition in ((0, 1), (1, 4)):
                value = current[state, j - 1, 0] + transitions[transition]
                if value > best:
                    best, predecessor = value, state
            current[1, j, 0] = best + insertion if best > NEG_INF // 2 else NEG_INF
            for field in range(1, 4):
                current[1, j, field] = current[predecessor, j - 1, field]

            best, predecessor = NEG_INF, 0
            for state, transition in ((0, 2), (2, 6)):
                value = previous[state, j, 0] + transitions[transition]
                if value > best:
                    best, predecessor = value, state
            current[2, j, 0] = best
            for field in range(1, 4):
                current[2, j, field] = previous[predecessor, j, field]

            if current[0, j, 0] > result[0]:
                result[0] = current[0, j, 0]
                result[1], result[2] = current[0, j, 1], i
                result[3], result[4] = current[0, j, 2], j
                result[5] = current[0, j, 3]
        previous, current = current, previous
    return result


@dataclass(frozen=True)
class ProfileAlignmentSpan:
    score: int
    model_start: int
    model_end: int
    target_start: int
    target_end: int
    matched_known_residues: int
    model_length: int
    target_length: int

    @property
    def passes_guard(self):
        return (self.score > 0 and min(self.model_length, self.target_length) > 0
                and 5 * self.matched_known_residues >= 4 * min(self.model_length, self.target_length))

    @property
    def category(self):
        if self.score <= 0:
            return "no_positive_alignment"
        model = 5 * self.matched_known_residues >= 4 * self.model_length
        target = 5 * self.matched_known_residues >= 4 * self.target_length
        if model and target:
            return "both_sides_covered"
        if model or target:
            return "asymmetric_coverage"
        return "partial_shared_region"

    def as_record(self):
        return {
            **asdict(self), "passes_guard": self.passes_guard, "category": self.category,
            "aligned_model_coverage": self.matched_known_residues / self.model_length if self.model_length else 0.,
            "aligned_target_coverage": self.matched_known_residues / self.target_length if self.target_length else 0.,
            "bounding_model_coverage": (self.model_end - self.model_start) / self.model_length if self.model_length else 0.,
            "bounding_target_coverage": (self.target_end - self.target_start) / self.target_length if self.target_length else 0.,
            "coordinate_system": "zero_based_half_open",
            "interpretation": "Representative Viterbi path evidence, not orthology truth or domain architecture",
        }


def _integer_array(value, shape, lower, upper, name, dtype):
    array = np.asarray(value)
    if (array.shape != shape or array.dtype.kind not in "iu"
            or np.any(array < lower) or np.any(array > upper)):
        raise ValueError("Invalid " + name + " shape or integer values")
    return np.ascontiguousarray(array, dtype=dtype)


def profile_alignment_span(profile, target, band_width=64, expected_score=None):
    """Trace a selected profile/target pair, optionally requiring score equality."""
    if (isinstance(band_width, (bool, np.bool_)) or not isinstance(band_width, Integral)
            or not np.iinfo(np.int32).min <= band_width <= np.iinfo(np.int32).max):
        raise ValueError("Band width must be an integer")
    length = profile.length
    if (isinstance(length, (bool, np.bool_)) or not isinstance(length, Integral)
            or not 0 <= length < np.iinfo(np.int32).max):
        raise ValueError("Invalid profile length")
    target = np.asarray(target)
    if target.ndim != 1:
        raise ValueError("Target must be a one-dimensional encoded sequence")
    emissions = _integer_array(profile.match_emissions, (length, 20), -128, 127, "match emissions", np.int8)
    insertions = _integer_array(profile.insert_emissions, (20,), -128, 127, "insert emissions", np.int8)
    transitions = _integer_array(profile.transitions, (7,), -10000, 0, "transitions", np.int32)
    target = _integer_array(target, target.shape, 0, 255, "target", np.uint8)
    if expected_score is not None and (
        isinstance(expected_score, (bool, np.bool_)) or not isinstance(expected_score, Real)
        or not np.isfinite(expected_score) or expected_score < 0 or int(expected_score) != expected_score
    ):
        raise ValueError("Expected score must be a finite nonnegative integer value")
    raw = _span_one(emissions, insertions, transitions, target, int(band_width))
    result = ProfileAlignmentSpan(*[int(value) for value in raw], int(length), len(target))
    if expected_score is not None and result.score != expected_score:
        raise ValueError(f"Profile span score mismatch: expected {expected_score}, observed {result.score}")
    if result.score > 0 and not (
        0 <= result.model_start < result.model_end <= length
        and 0 <= result.target_start < result.target_end <= len(target)
        and 0 <= result.matched_known_residues <= min(
            result.model_end - result.model_start, result.target_end - result.target_start)
    ):
        raise ValueError("Invalid representative Viterbi span coordinates/count")
    return result
