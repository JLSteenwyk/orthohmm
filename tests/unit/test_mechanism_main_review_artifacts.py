"""Read back the revised presentation without repeating scientific inference."""

import hashlib
import json
from pathlib import Path
import subprocess

import fitz


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "benchmark_tools/results"
REVISION = "45af9c000485ea7003a046f74ef582cf7200237e"


def read(name):
    return json.loads((RESULTS / name).read_text())


def identity(data):
    return len(data), hashlib.sha256(data).hexdigest()


def retained(reference):
    path = Path(reference["path"])
    if not path.is_absolute():
        path = ROOT / path
    assert identity(path.read_bytes()) == (reference["bytes"], reference["sha256"])
    return path


def test_render_time_inputs_match_the_actual_commit_not_future_ledger_edits():
    review = read("publication_main_mechanism_visual_review_20261004.json")
    render = read("publication_main_render_20261004_v2.json")
    refs = render["sources"] + render["targets"]
    expected = {str(Path(r["path"]).relative_to(ROOT)): r for r in refs
                if Path(r["path"]).is_relative_to(ROOT)}
    observed = {r["path"]: r for r in review["checked_committed_direct_inputs"]}
    assert expected.keys() == observed.keys()
    assert review["render_time_source_revision"] == REVISION
    for path, ref in observed.items():
        assert ref["revision"] == REVISION
        data = subprocess.check_output(["git", "show", REVISION + ":" + path], cwd=ROOT)
        assert identity(data) == (ref["bytes"], ref["sha256"])
        assert identity(data) == (expected[path]["bytes"], expected[path]["sha256"])
        blob = subprocess.check_output(["git", "rev-parse", REVISION + ":" + path],
                                       cwd=ROOT, text=True).strip()
        assert blob == ref["git_blob"]


def test_actual_render_print_bounds_and_review_scope_are_bound():
    review = read("publication_main_mechanism_visual_review_20261004.json")
    for key in ("render", "print", "bounds", "selection", "assembly", "pdf", "prior_figure_review"):
        retained(review[key])
    render = read("publication_main_render_20261004_v2.json")
    printing = read("publication_main_print_20261004_v2/print.json")
    bounds = read("publication_main_pdf_review_20261004_v2/report.json")
    assert len(render["citation_ids"]) == 18 and render["untracked_targets"] == []
    retained(render["html"])
    retained(printing["pdf"])
    assert printing["returncode"] == 0 and printing["page_count"] == 14
    assert bounds["page_count"] == 14 and bounds["bounds_violations"] == []
    assert review["main_pages_actually_viewed"] == list(range(1, 15))
    assert review["new_guide_pages_actually_viewed"] == [15, 16, 17, 18]
    assert review["observed_clipping_or_overlap"] is False
    assert review["all_sixteen_figures_fresh_manual_review_claimed"] is False
    assert review["rc3_archive_changed"] is review["publication_ready"] is False
    assert bounds["visual_review_complete"] is printing["visual_review_complete"] is False
    assert review["original_failed_bounds_selector"]["output_created"] is False


def test_all_thirty_actual_source_pages_retain_text_geometry_and_pixels():
    assembly = read("publication_main_with_figures_20261004_v2/assembly.json")
    assert [assembly[k] for k in ("main_pages", "guide_pages", "figure_pages", "total_pages")] == [14, 4, 16, 34]
    assert len(assembly["preserved_source_pages"]) == 30 and assembly["mupdf_warnings"] == ""
    with fitz.open(retained(assembly["pdf"])) as combined:
        assert len(combined) == 34
        for ref in assembly["preserved_source_pages"]:
            with fitz.open(retained(ref["source"])) as source:
                before = source[ref["original_page"] - 1]
                after = combined[ref["combined_page"] - 1]
                assert before.rect == after.rect and before.get_text("words") == after.get_text("words")
                a, b = before.get_pixmap(alpha=False), after.get_pixmap(alpha=False)
                assert (a.width, a.height, a.n, a.samples) == (b.width, b.height, b.n, b.samples)
                assert hashlib.sha256(a.samples).hexdigest() == ref["pixel_sha256"]


def test_figures_captions_and_existing_review_pdf_are_not_replaced():
    old = read("publication_figure_selection_20261004.json")
    new = read("publication_figure_selection_20261004_v2.json")
    assert new["figures"] == old["figures"] and len(new["figures"]) == 16
    assert new["provenance"][:len(old["provenance"])] == old["provenance"]
    retained(read("publication_main_with_figures_20261004/assembly.json")["pdf"])
    for ref in [new["main"], *new["provenance"], *[f["pdf"] for f in new["figures"]]]:
        retained(ref)


def test_bookmarks_and_guide_and_main_links_reach_correct_figure_pages():
    assembly = read("publication_main_with_figures_20261004_v2/assembly.json")
    with fitz.open(retained(assembly["pdf"])) as combined:
        expected = [[1, "Main Text", 1], [1, "Figure Guide", 15]] + [
            [1, "A%d. %s" % (f["number"], f["title"]), 18 + f["number"]]
            for f in assembly["figures"]]
        assert combined.get_toc() == expected
        assert len(assembly["redirected_main_figure_links"]) == 7
        assert len(assembly["restored_original_file_uri_actions"]) == 63
        for group in range(4):
            page = combined[14 + group]
            assert "Main text: pages 1-14." in page.get_text()
            assert [r["page"] for r in page.get_links() if r["kind"] == fitz.LINK_GOTO] == list(
                range(18 + group * 4, 22 + group * 4))
        for index in range(14):
            assert sorted(r["page"] + 1 for r in combined[index].get_links()
                          if r["kind"] == fitz.LINK_GOTO) == sorted(
                r["combined_page"] for r in assembly["redirected_main_figure_links"]
                if r["main_page"] == index + 1)


def test_new_mechanism_claims_and_all_resource_tables_survive_rendering():
    with fitz.open(RESULTS / "publication_main_print_20261004_v2/document.pdf") as pdf:
        text = " ".join(" ".join(p.get_text() for p in pdf).split())
    replay = read("simulation_mechanism_reporting_replay_20261004.json")["summary"]
    conditions = {c["condition"]: c for c in replay["conditions"]}
    for key in ("divergent", "divergent_turnover"):
        cell = conditions[key]
        assert format(cell["cross_candidate_fn"], ",") in text
        assert format(cell["oracle_fn"], ",") in text
        assert "%.3f%%" % (100 * cell["cross_candidate_fn"] / cell["oracle_fn"]) in text
        for metric in ("different_graph_components", "connected_but_separated", "direct_graph_edge_but_separated"):
            assert format(cell["upstream_totals"][metric], ",") in text
    for phrase in ("380 TP, 292 TN, 62 FN and 66 FP", "All 62 false negatives",
                   "All 66 have overlap evidence", "Contention distortion is unknown and potentially method dependent",
                   "not a comparison against OrthoFinder", "References",
                   "All 27 planned identities are terminal-reviewed", "not 27 eligible measurements",
                   "No defaults or endpoints changed"):
        assert phrase in text
    assert text.count("Unavailable") == 9
    section = (RESULTS / "threadripper_shared_resource_section_20261004_v27/resource_section.md").read_text()
    compact = "".join(text.split())
    table_lines = [line for line in section.splitlines() if line.startswith("|")]
    numeric_cells = [cell.strip() for line in table_lines for cell in line.split("|")
                     if cell.strip() and cell.strip()[0].isdigit()]
    assert numeric_cells
    for cell in numeric_cells:
        assert "".join(cell.split()) in compact
