"""Actual printed content and exact scope of the separately observed visual review."""

import hashlib
import json
from pathlib import Path

import fitz

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "benchmark_tools/results"


def load(name):
    return json.loads((BASE / name).read_text())


def check(ref):
    path = Path(ref["path"])
    if not path.is_absolute():
        path = ROOT / path
    assert path.stat().st_size == ref["bytes"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == ref["sha256"]


def test_new_render_and_review_use_stable_assets_and_preserve_frozen_sources():
    render = load("publication_main_render_20261007_v3.json")
    assert render["status"] == "manuscript_review_rendered"
    assert render["local_occurrences"] == 117 and render["unique_targets"] == 114
    assert render["untracked_targets"] == [] and "iqtree3_2026" in render["citation_ids"]
    assert all(not r["path"].endswith("PUBLICATION_PROGRESS.md") for r in render["targets"])
    assert any(r["path"].endswith("PUBLICATION_MAIN_TEXT_20261007_v2.md") for r in render["targets"])
    for ref in [render["html"], *render["sources"], *render["targets"]]:
        check(ref)
    review = load("publication_main_pdf_review_20261007_v3/report.json")
    assert review["page_count"] == 21 and review["bounds_violations"] == []
    assert review["visual_review_complete"] is False
    assert review["publication_ready"] is False
    for ref in review["checked_records"]:
        check(ref)


def test_pdf_retains_all_six_distance_rows_methods_limits_and_citation():
    printed = load("publication_main_print_20261007_v3/print.json")
    assert printed["status"] == "verified_html_printed" and printed["page_count"] == 21
    check(printed["pdf"])
    with fitz.open(printed["pdf"]["path"]) as document:
        text = " ".join(" ".join(p.get_text().split()) for p in document)
        for row in (
            "All families 18 R at P0/C0 +10.039 +30.521 -6.533",
            "All families 18 C at P0/R0 -0.347 -3.523 +4.375",
            "Lower/equal distance 9 R at P0/C0 +11.945 +32.417 -7.950",
            "Lower/equal distance 9 C at P0/R0 -1.245 -4.547 +5.233",
            "Higher distance 9 R at P0/C0 +8.427 +28.625 -5.117",
            "Higher distance 9 C at P0/R0 +0.324 -2.500 +3.517",
        ):
            assert row in text
        for phrase in ("fixed WAG+G4", "all 10,765 distances", "2.24582463005 substitutions/site",
                       "not calibrated biological time", "unknown and potentially tool-dependent impact",
                       "IQ-TREE 3: Phylogenomic Inference Software Using Complex",
                       "10.1093/molbev/msag117"):
            assert phrase in text
        for page in (document[9], document[10]):
            assert "Delta Precision" in page.get_text() and "Delta Recall" in page.get_text()


def test_visual_receipt_binds_all_actual_viewed_pages_without_readiness_claim():
    manual = load("publication_main_pdf_visual_review_20261007_v3.json")
    automated = load("publication_main_pdf_review_20261007_v3/report.json")
    assert manual["manual_pages_actually_viewed"] == list(range(1, 22))
    assert manual["viewed_rasters"] == automated["rendered_pages"]
    assert len(manual["viewed_rasters"]) == 21
    assert manual["visual_review_complete"] is True
    for ref in [*manual["records"], manual["pdf"], *manual["viewed_rasters"]]:
        check(ref)
    for flag in ("publication_ready", "new_accuracy_or_resource_admission",
                 "scientific_inference_reexecuted", "old_render_or_archive_rebuilt"):
        assert manual[flag] is False
