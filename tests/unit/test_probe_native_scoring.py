from copy import deepcopy

import numpy as np
import pytest

from benchmark_tools.probe_native_scoring import BANDS, LENGTHS, fixtures, flat, compare


def report():
    return {"status": "checked", "fixture_sha256": "fixture", "probe_sha256": "probe",
            "commit": "commit", "pairs": 2, "packages": {"numpy": "same"}, "backend": "scalar_C",
            "bands": [{"band": band, "scores": [1, 2], "scalar_scores": [1, 2], "jit_scores": [1, 2],
                       "evalues": [1., 1e-6], "normalized": [.1, .2],
                       "decisions": {str(t): [False, True] for t in (1e-3, 1e-4, 1e-5)},
                       "scalar_mismatches": 0, "jit_mismatches": 0} for band in BANDS]}


def test_fixtures_reproducible_and_cover_lengths():
    a, b = fixtures(), fixtures()
    for left, right in zip(a[:2], b[:2]):
        for x, y in zip(left, right):
            np.testing.assert_array_equal(x, y)
    np.testing.assert_array_equal(a[2], b[2])
    assert [len(q) for q in a[0]] == list(LENGTHS)
    assert len(a[2]) == 405
    assert any(20 in q for q in a[0]) and any(20 in t for t in a[1])
    values, offsets, lengths = flat(a[1])
    for sequence, offset, length in zip(a[1], offsets, lengths):
        np.testing.assert_array_equal(values[offset:offset + length], sequence)


def test_matching_scores():
    a, b = report(), report()
    b["backend"] = "multipair_AVX2"
    b["bands"][0]["evalues"][0] += 1e-14
    assert compare(a, b)["end_to_end_equivalence"] is False


@pytest.mark.parametrize("field", ["fixture_sha256", "probe_sha256", "commit", "pairs", "packages"])
def test_provenance_mismatch(field):
    a, b = report(), report()
    b[field] = "different"
    with pytest.raises(ValueError, match="provenance"):
        compare(a, b)


@pytest.mark.parametrize("problem", ["integer", "jit", "floating", "nan", "decisions", "bands", "status", "counter", "length"])
def test_invalid_evidence(problem):
    a, b = report(), deepcopy(report())
    row = b["bands"][0]
    if problem == "integer":
        row["scores"][0] += 1
    elif problem == "jit":
        row["jit_scores"][0] += 1
    elif problem == "floating":
        row["normalized"][0] += .01
    elif problem == "nan":
        row["normalized"][0] = float("nan")
    elif problem == "decisions":
        row["decisions"] = {}
    elif problem == "bands":
        b["bands"].pop()
    elif problem == "status":
        b["status"] = "mismatch"
    elif problem == "counter":
        row["jit_mismatches"] = 1
    else:
        row["normalized"].pop()
    with pytest.raises(ValueError):
        compare(a, b)
