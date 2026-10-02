"""Read back retained synthetic scores; never rerun native inference in CI."""

import gzip
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import beta

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "benchmark_tools/results"
REGIMES = ("blosum_background", "uniform", "half_glutamine")
LENGTHS = (50, 150, 400)
EXPECTED_HITS = {
    ("blosum_background", 50): (0, 0),
    ("blosum_background", 150): (2, 1),
    ("blosum_background", 400): (0, 0),
    ("uniform", 50): (1, 1),
    ("uniform", 150): (7, 5),
    ("uniform", 400): (12, 3),
    ("half_glutamine", 50): (8983, 8983),
    ("half_glutamine", 150): (10000, 10000),
    ("half_glutamine", 400): (10000, 10000),
}


@pytest.fixture(scope="module")
def retained():
    receipt = json.loads((RESULTS / "frozen_null_score_receipt_20261002.json").read_bytes())
    pins = {p["path"]: p for p in receipt["files"]}
    name = "benchmark_tools/results/frozen_null_score_observations_20261002.json.gz"
    compressed = (ROOT / name).read_bytes()
    assert len(compressed) == pins[name]["bytes"] == 172690
    assert hashlib.sha256(compressed).hexdigest() == pins[name]["sha256"]
    raw = gzip.decompress(compressed)
    original = pins["benchmarks/work/frozen_null_score_panel_20261002/result.json"]
    assert len(raw) == original["bytes"] == 2657291
    assert hashlib.sha256(raw).hexdigest() == original["sha256"]
    return json.loads(raw), receipt


def test_retained_design_provenance_and_scope(retained):
    report, receipt = retained
    assert report["status"] == "prespecified_frozen_null_panel_completed"
    assert report["revision"] == receipt["scientific_revision"]
    assert receipt["protocol_commit"] == "186bd5e47ecaa9361a86064c33d7a57c70ce1ce3"
    assert report["independent_pairs"] == 90000
    assert report["native_score_evaluations"] == 180000
    assert report["tail_endpoints"] == 90
    expected = {(r, length, seed) for r in REGIMES for length in LENGTHS for seed in range(10)}
    rows = {(r["regime"], r["length"], r["seed"]) for r in report["rows"]}
    assert len(report["rows"]) == 90 and rows == expected
    assert len(report["summaries"]) == 9
    assert {(s["regime"], s["length"]) for s in report["summaries"]} == set(EXPECTED_HITS)
    for source, name in (("source", "benchmark_tools/probe_frozen_null_scores.py"),
                         ("protocol", "benchmark_tools/results/FROZEN_NULL_SCORE_PROTOCOL_20261002.md")):
        raw = (ROOT / name).read_bytes()
        assert len(raw) == report[source]["bytes"]
        assert hashlib.sha256(raw).hexdigest() == report[source]["sha256"]
    for key in ("native_pipeline_rerun", "prefilter_executed", "calibration_established",
                "benchmark_scores_or_defaults_changed", "publication_ready"):
        assert report[key] is False
    audit = receipt["audit"]
    assert audit["reference_python_scores_rechecked"] == 180
    assert audit["sequence_digests_regenerated"] == 180
    assert audit["all_endpoint_counts_and_intervals_recomputed"] is True
    assert audit["all_native_scores_recomputed"] is False
    assert audit["calibration_established"] is False


@pytest.mark.parametrize("regime,length", EXPECTED_HITS)
def test_all_retained_endpoints_against_raw_scores(retained, regime, length):
    report, _ = retained
    rows = sorted((r for r in report["rows"] if (r["regime"], r["length"]) == (regime, length)),
                  key=lambda r: r["seed"])
    summary = next(s for s in report["summaries"] if (s["regime"], s["length"]) == (regime, length))
    scores = {}
    for band in ("0", "64"):
        for row in rows:
            assert set(row["scores"]) == {"0", "64"}
            assert len(row["scores"][band]) == 1000
            assert all(type(s) is int and s >= 0 for s in row["scores"][band])
        scores[band] = np.array([s for row in rows for s in row["scores"][band]])
        tails = summary["bands"][band]
        assert [t["threshold"] for t in tails] == [1, 0.1, 0.01, 0.001, 0.0001]
        for tail in tails:
            boundary = 1
            while 0.134 * length * length * math.exp(-0.3176 * boundary) >= tail["threshold"]:
                boundary += 1
            assert tail["minimum_integer_score"] == boundary
            k = int(np.count_nonzero(scores[band] >= boundary))
            assert tail["hits"] == k and tail["trials"] == 10000
            assert tail["fraction"] == k / 10000
            e = 0.134 * length * length * math.exp(-0.3176 * boundary)
            assert tail["boundary_approximate_e"] == pytest.approx(e, rel=1e-12, abs=0)
            model = -math.expm1(-e)
            assert tail["poisson_model_tail_reference"] == pytest.approx(model, rel=1e-12, abs=0)
            assert tail["observed_to_model_tail_ratio"] == pytest.approx(k / 10000 / model)
            for field, alpha in (("nominal_clopper_pearson", 0.05),
                                 ("bonferroni_clopper_pearson", 0.05 / 90)):
                low = 0.0 if k == 0 else beta.ppf(alpha / 2, k, 10001 - k)
                high = 1.0 if k == 10000 else beta.ppf(1 - alpha / 2, k + 1, 10000 - k)
                assert tail[field] == pytest.approx([low, high], rel=1e-12, abs=1e-10)
        assert tails[-1]["hits"] == EXPECTED_HITS[regime, length][int(band == "64")]
    assert (scores["64"] <= scores["0"]).all()
    assert summary["band_changed_scores"] == int(np.count_nonzero(scores["64"] != scores["0"]))
    for tail in summary["bands"]["0"]:
        gate = tail["minimum_integer_score"]
        lost = np.count_nonzero((scores["0"] >= gate) & (scores["64"] < gate))
        assert summary["band_lost_gate_hits"][str(tail["threshold"])] == int(lost)
