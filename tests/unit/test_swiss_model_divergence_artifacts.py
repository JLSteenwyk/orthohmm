"""Retained complete evidence, with no rerun of inference or raw benchmark scoring."""

import hashlib
import json
from pathlib import Path
import tarfile

import pytest

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "benchmark_tools/results"
RECEIPT = json.loads((RESULTS / "swiss_model_divergence_execution_23932_v1.json").read_text())


def sha(data):
    return hashlib.sha256(data).hexdigest()


@pytest.mark.parametrize("key", ["feature_readback", "projection", "projection_readback", "evidence_archive"])
def test_actual_artifact_hashes_and_sizes_match_execution_receipt(key):
    ref = RECEIPT[key]
    data = (ROOT / ref["path"]).read_bytes()
    assert len(data) == ref["bytes"] and sha(data) == ref["sha256"]


def test_archive_contains_every_selected_alignment_and_new_attempt_payload_exactly():
    readback = json.loads((ROOT / RECEIPT["feature_readback"]["path"]).read_text())
    marker = "/benchmarks/results/"
    expected = {r["path"].split(marker, 1)[1]: r for r in readback["checked_inputs"] if marker in r["path"]}
    assert len(expected) == 183
    with tarfile.open(ROOT / RECEIPT["evidence_archive"]["path"], "r:gz") as archive:
        files = [m for m in archive.getmembers() if m.isfile()]
        assert len(files) == len(set(m.name for m in files)) == len(expected)
        assert {m.name for m in files} == set(expected)
        assert not any(m.issym() or m.islnk() or m.name.startswith("/") or ".." in Path(m.name).parts
                       for m in archive.getmembers())
        for member in files:
            data = archive.extractfile(member).read()
            ref = expected[member.name]
            assert len(data) == ref["bytes"] and sha(data) == ref["sha256"]
        feature_path = "swiss_model_divergence_20261007_v1/report.json"
        feature_bytes = archive.extractfile(feature_path).read()
        assert sha(feature_bytes) == RECEIPT["feature_report"]["sha256"]
        feature = json.loads(feature_bytes)
    assert feature["failed_families"] == []
    assert len(feature["runs"]) == 18
    assert all(r["exit_code"] == 0 and r["timed_out"] is False for r in feature["runs"])
    assert sum(r["features"]["pairs"] for r in feature["runs"]) == 10765
    assert feature["memberships"] == json.loads((ROOT / RECEIPT["projection"]["path"]).read_text())["memberships"]


def test_reader_covers_full_population_and_both_bins_without_invented_scope():
    feature = json.loads((ROOT / RECEIPT["feature_readback"]["path"]).read_text())
    projection = json.loads((ROOT / RECEIPT["projection_readback"]["path"]).read_text())
    assert feature["status"] == "features_verified"
    assert (feature["families_checked"], feature["proteins_checked"], feature["pairs_checked"]) == (18, 563, 10765)
    assert len(feature["checked_inputs"]) == 190
    assert len(feature["strata"]["lower_or_equal_median"]) == len(feature["strata"]["higher_than_median"]) == 9
    assert projection["status"] == "projection_verified"
    assert (projection["family_rows_checked"], projection["score_rows_checked"], projection["differences_checked"]) == (54, 9, 6)
    for doc in (feature, projection, RECEIPT):
        assert doc["new_bootstrap_draws"] == 0
        for flag in ("new_uncertainty", "independent_confirmation", "publication_ready",
                     "new_accuracy_or_resource_admission", "scientific_timings_admitted"):
            assert doc[flag] is False


def test_every_retained_table_binding_is_current_and_semantics_unchanged():
    report = json.loads((ROOT / RECEIPT["projection"]["path"]).read_text())
    for ref in report["outputs"].values():
        path = ROOT / "benchmark_tools/results/swiss_model_divergence_strata_20261007_v1" / Path(ref["path"]).name
        data = path.read_bytes()
        assert len(data) == ref["bytes"] and sha(data) == ref["sha256"]
    assert [r["stratum"] for r in report["rows"]] == ["all"] * 3 + ["lower_or_equal_median"] * 3 + ["higher_than_median"] * 3
    assert report["cells"][1]["timing_eligible"] is False
    assert report["cells"][1]["timing_admitted"] is False
    assert report["rows"][1]["prediction_semantics"] == "native_pair"
    assert all(r["prediction_semantics"] == "group_clique" for r in report["rows"] if r["cell"] != "p0_c0_r1")


def test_preexecution_sources_and_protocol_have_not_changed():
    feature_reader = json.loads((ROOT / RECEIPT["feature_readback"]["path"]).read_text())
    projection = json.loads((ROOT / RECEIPT["projection"]["path"]).read_text())
    refs = [r for r in feature_reader["checked_inputs"] if "/benchmark_tools/" in r["path"]
            and Path(r["path"]).name in ("prepare_swiss_model_divergence.py", "swiss_model_divergence_batch_20261007.sh",
                                         "SWISS_MODEL_DIVERGENCE_PROTOCOL_20261007.md")]
    refs += [feature_reader["source"], projection["source"]]
    for ref in refs:
        relative = ref["path"].split("/benchmark_tools/", 1)[1]
        data = (ROOT / "benchmark_tools" / relative).read_bytes()
        assert sha(data) == ref["sha256"] and len(data) == ref["bytes"]


def test_incomplete_auxiliary_resource_accounting_is_not_zero_or_native_timing():
    assert RECEIPT["scheduler"]["total_cpu_seconds"] is None
    assert RECEIPT["scheduler"]["maximum_rss_bytes"] is None
    assert RECEIPT["scientific_timings_admitted"] is False
    assert RECEIPT["existing_native_jobs_unchanged"] == [23902, 23910]
    for key in ("feature_readback", "projection", "projection_readback"):
        assert RECEIPT[key]["exit_code"] == 0 and RECEIPT[key]["elapsed_seconds"] > 0
