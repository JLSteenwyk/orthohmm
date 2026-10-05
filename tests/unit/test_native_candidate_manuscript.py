"""Keep the new fixed-boundary result distinct from earlier recovery diagnostics."""

import json
from pathlib import Path


BASE = Path(__file__).resolve().parents[2] / "benchmark_tools/results"


def section():
    text = (BASE / "PUBLICATION_MANUSCRIPT_DRAFT_20260916.md").read_text()
    return text.split("A separate [full-native candidate-only diagnostic]", 1)[1].split(
        "A subsequent [fresh recovery installation]", 1)[0]


def test_candidate_result_numbers_are_bound_to_actual_report():
    report = json.loads((BASE / "native_candidate_factorial_22429/report.json").read_text())
    summary = json.loads((BASE / "native_candidate_factorial_22429/summary.json").read_text())
    paragraph = section()
    assert report["controls"] == {"historical": True, "fresh_full": True}
    assert report["genes"] == 251378
    assert len(report["rows"]) == 5
    assert all(row["groups"] == 54745 for row in report["rows"])
    assert [r["genes_changed_vs_historical"] for r in summary["rows"]] == [0, 0, 58, 58, 58]
    assert [r["genes_changed_vs_native"] for r in summary["rows"]] == [58, 58, 0, 0, 0]
    for literal in ("54,745", "251,378", "58-gene", "P0C1R0", "ten pairwise comparisons"):
        assert literal in paragraph


def test_new_case_does_not_replace_old_case_or_promote_defaults():
    paragraph = section()
    for literal in ("93-gene dependency-control", "profile expansion disabled",
                    "reconciliation disabled", "fixed-boundary effect", "No default was changed",
                    "not comparative timing", "unknown and potentially tool-dependent"):
        assert literal in paragraph
    claims = (BASE / "PUBLICATION_CLAIMS_20260916.md").read_text()
    assert "remaining 93-gene candidate discrepancy" in claims
    assert "58-gene difference" in claims
    assert "original upstream-order provenance or a new accuracy advantage" in claims


def test_new_local_evidence_links_exist():
    for name in ("NATIVE_CANDIDATE_FACTORIAL_RESULT_22429.md", "native_candidate_factorial_22429/summary.md"):
        assert (BASE / name).is_file()
        assert name in section()
        assert name in (BASE / "PUBLICATION_CLAIMS_20260916.md").read_text()


def test_inherited_unchanged_accuracy_is_not_new_scoring():
    score = json.loads((BASE / "native_factorial_orthobench_score_22429.json").read_text())
    report = json.loads((BASE / "native_candidate_factorial_22429/report.json").read_text())
    assert all(v == 0 for v in score["difference_from_original_percentage_points"].values())
    assert score["canonical_partition_comparison"]["reference_touching_groups_equal"]
    assert score["canonical_partition_comparison"]["reference_genes_in_changed_groups"] == []
    assert report["accuracy_scored"] is False
    assert "separately evaluated OrthoBench statistic was unchanged" in section()
