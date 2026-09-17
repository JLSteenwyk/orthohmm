from copy import deepcopy

import pytest

from benchmark_tools.probe_native_scoring import BANDS, fixtures, digest
from benchmark_tools.summarize_native_scoring import summarize


def evidence():
    queries, targets, pairs = fixtures()
    fixture = {"queries": [q.tolist() for q in queries], "targets": [t.tolist() for t in targets],
               "pairs": pairs.tolist(), "bands": list(BANDS)}
    n = len(pairs)
    return {"status": "checked", "fixture_sha256": digest(fixture), "probe_sha256": "same",
            "pairs": n, "commit": "same", "packages": {}, "backend": "scalar_C",
            "bands": [{"band": b, "scores": [0] * n, "scalar_scores": [0] * n, "jit_scores": [0] * n,
                       "evalues": [1.] * n, "normalized": [0.] * n, "scalar_mismatches": 0,
                       "jit_mismatches": 0, "decisions": {str(t): [False] * n for t in (1e-3, 1e-4, 1e-5)}}
                      for b in BANDS]}


def test_diagnostics_never_admit_timing():
    a = evidence()
    result = summarize(a, deepcopy(a))
    assert result["admitted"] is False
    assert "gate_result" in result


def test_mismatch_retained_with_identity():
    a = evidence()
    b = deepcopy(a)
    b["status"] = "mismatch"
    b["bands"][0]["scores"][0] = 4
    result = summarize(a, b)
    assert result["gate_error"] == "Local backend discrepancy"
    assert result["bands"][0]["native_score_differences"][0]["pair_index"] == 0


@pytest.mark.parametrize("field,value", [("fixture_sha256", "bad"), ("probe_sha256", "bad"),
                                        ("pairs", 2), ("bands", [])])
def test_changed_diagnostic_evidence(field, value):
    a, b = evidence(), evidence()
    b[field] = value
    with pytest.raises(ValueError):
        summarize(a, b)
