"""Read back actual native-evidence review artifacts, not fixture inference."""

import hashlib
import json
from pathlib import Path
import subprocess

import fitz


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "benchmark_tools/results"


def read(name):
    return json.loads((BASE / name).read_text())


def identity(data):
    return len(data), hashlib.sha256(data).hexdigest()


def retained(reference):
    path = Path(reference["path"])
    assert identity(path.read_bytes()) == (reference["bytes"], reference["sha256"])
    return path


def pdf_text():
    review = read("publication_main_native_visual_review_20261006.json")
    with fitz.open(retained(review["pdf"])) as document:
        return " ".join(" ".join(page.get_text() for page in document).split())


def test_actual_render_sources_match_source_commit_before_later_ledger_updates():
    review = read("publication_main_native_visual_review_20261006.json")
    render = read("publication_main_render_20261006_v1.json")
    revision = review["source_revision"]
    assert revision == "47b7f14ba0c2e8c5b7f070f5b11620dd356a1eb9"
    for reference in [*render["sources"], *render["targets"]]:
        path = Path(reference["path"])
        if path.is_relative_to(ROOT):
            data = subprocess.check_output([
                "git", "show", revision + ":" + str(path.relative_to(ROOT))], cwd=ROOT)
            assert identity(data) == (reference["bytes"], reference["sha256"])
        else:
            retained(reference)
    assert len(render["citation_ids"]) == 18
    assert render["unique_targets"] == 85 and render["local_occurrences"] == 88
    assert render["untracked_targets"] == []
    assert render["stderr"] == {"parse": "", "render": ""}


def test_print_bounds_and_manual_scope_remain_distinct():
    review = read("publication_main_native_visual_review_20261006.json")
    for key in ("render", "print", "bounds", "pdf", "source", "original_main", "original_combined_review"):
        retained(review[key])
    printing = read("publication_main_print_20261006_v1/print.json")
    bounds = read("publication_main_pdf_review_20261006_v1/report.json")
    assert printing["returncode"] == 0 and printing["page_count"] == bounds["page_count"] == 18
    assert bounds["bounds_violations"] == []
    pages = [3, 4, 5, 6, 7, 8, 10, 11, 12, 14, 15, 16, 17, 18]
    assert review["main_pages_actually_viewed"] == pages
    assert review["rendered_page_records"] == bounds["rendered_pages"]
    for page, ref in zip(pages, review["rendered_page_records"]):
        assert retained(ref).name == f"page_{page:03d}.png"
    for flag in ("observed_clipping_or_incoherent_overlap", "all_eighteen_main_pages_fresh_manual_review_claimed",
                 "linked_figures_fresh_manual_review_claimed", "combined_figure_assembly_created",
                 "native_inference_or_scoring_repeated", "old_archive_or_pdf_changed", "publication_ready"):
        assert review[flag] is False
    assert printing["visual_review_complete"] is bounds["visual_review_complete"] is False
    assert printing["publication_ready"] is bounds["publication_ready"] is False
    assert review["new_bootstrap_draws"] == 0


def test_all_native_point_tables_and_conditional_effects_survive_pdf_rendering():
    text = pdf_text()
    ob = read("native_factorial_progress_20261005_v6/report.json")
    for row in ob["rows"]:
        if row["dataset"] == "orthobench":
            assert row["cell"].upper().replace("_", "/") in text
            for metric in ("f1", "precision", "recall"):
                assert f'{100 * row[metric]:.4f}' in text
    qfo = read("native_qfo_scientific_scores_20261006_v1/report.json")
    for row in qfo["rows"]:
        if row["accuracy_admitted"]:
            for value in row["scores"].values():
                assert f'{value:.8f}' in text
            assert f'{row["secondary_mean"]:.8f}' in text
            assert f'{100 * row["relation_coverage"]:.4f}%' in text
    binding = read("native_factorial_uncertainty_binding_20261005_v4.json")
    for row in binding["contrasts"]:
        if row["metrics"] is not None:
            metric = row["metrics"]["f_score"]
            assert f'{metric["difference_percentage_points"]:+.4f}' in text
            low, high = metric["bonferroni_percentile_ci"]
            assert f'[{low:.4f}, {high:.4f}]' in text
    binding = read("native_qfo_swiss_uncertainty_binding_22449_20261006.json")
    for row in binding["contrasts"]:
        if row["metrics"] is not None:
            for metric in row["metrics"].values():
                assert f'{100 * metric["difference"]:+.4f}' in text
                low, high = metric["bonferroni_percentile_ci"]
                assert f'[{100 * low:.4f}, {100 * high:.4f}]' in text


def test_native_table_page_break_retains_headers_and_all_rows():
    review = read("publication_main_native_visual_review_20261006.json")
    with fitz.open(retained(review["pdf"])) as document:
        first, second = (" ".join(document[i].get_text().split()) for i in (5, 6))
        for text in (first, second):
            assert "Native Cell Group-Recovery F1 (%) Precision (%) Recall (%)" in text
        assert "P0/C0/R0 69.7634" in first
        for label in ("P0/C0/R1", "P0/C1/R0", "P0/C1/R1", "P1/C0/R1", "P1/C1/R0"):
            assert label in second
        for endpoint in ("VGNC", "SwissTrees", "TreeFam-A", "GO", "EC", "FAS"):
            assert endpoint in second
    assert review["phrase_matches"]["Fresh Native QfO Trade-Off And Pair Composition"] == [7]


def test_broad_science_failure_and_shared_host_scope_survive_pdf():
    text = pdf_text()
    for phrase in (
        "Both intervals include zero", "not the total HMM contribution",
        "five remain unavailable", "other 13 planned contrasts remain unavailable",
        "failed 22437 timing", "not successful timing", "F1 interval includes zero",
        "not the selected-default all-tool comparison", "14/1/3", "817,432 rows",
        "0.3989%", "2.6358%", "Intersection-only means would change the endpoint",
        "paired uncertainty remains unfinished", "unknown and potentially tool-dependent impact",
        "not estimates of isolated performance", "No background overhead is subtracted",
        "27 of 27 reviewed attempts", "24 eligible observations", "not 27 eligible measurements",
        "Gene-Tree And Event-History Controls Localize Errors",
        "Transfer And Biological Recovery Reveal Trade-Offs",
        "Original TreeFam-A family mappings", "not untouched publication validation",
        "complete cached-execution costs remain unknown", "No submission-ready release or archival DOI is claimed",
        "neither archive includes this 6 October native-evidence revision",
        "did not consistently outperform full OrthoFinder", "Not submission-ready",
    ):
        assert phrase in text
    section = (BASE / "threadripper_shared_resource_section_20261004_v27/resource_section.md").read_text()
    compact = "".join(text.split())
    for line in section.splitlines():
        if line.startswith("| Ortho"):
            for cell in line.split("|")[1:-1]:
                assert "".join(cell.split()) in compact
    assert text.count("Unavailable") == 9
