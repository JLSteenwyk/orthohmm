import copy
import csv
import json
from pathlib import Path
import shutil
import sys

from PIL import Image
import pytest

from benchmark_tools import plot_native_qfo_three_cell_scores as plot
from benchmark_tools import review_native_qfo_three_cell_figure as reviewer
from benchmark_tools.prepare_ob_candidate_neighborhood import record

RESULTS = Path(__file__).resolve().parents[2] / "benchmark_tools/results"


def inputs():
    return [json.loads((RESULTS / name).read_text()) for name in (
        "native_qfo_scientific_scores_20261006_v2/report.json",
        "native_qfo_candidate_swiss_uncertainty_20261006_v1.json",
        "native_qfo_candidate_swiss_readback_20261006_v1.json")]


def write(path, value):
    path.write_text(json.dumps(value, sort_keys=True))
    return record(path)


def test_three_cell_points_keep_endpoints_coverage_and_missing_values_distinct():
    snapshot, binding, readback = inputs()
    rows, scores, coverage, intervals = plot.figure_data(snapshot, binding, readback)
    assert len(scores) == 18 and len(coverage) == 3 and len(intervals) == 6
    assert sum(r["statistic"] == "F1" for r in scores) == 9
    assert not any("mean" in r["endpoint"].lower() for r in scores)
    assert rows["p0_c0_r1"]["resources"] is None
    assert all(r["adjusted_low_pp"] < 0 < r["adjusted_high_pp"] for r in intervals if r["metric"] == "F1")
    assert coverage[2]["relation_coverage"] == 585610/984137
    assert all(v is None for v in rows["p0_c1_r1"]["scores"].values())


@pytest.mark.parametrize("value", [True, float("nan"), float("inf"), "0.5", -.1, 1.1])
def test_invalid_candidate_score_refused(value):
    snapshot, binding, readback = inputs()
    snapshot["rows"][2]["scores"]["FAS"] = value
    with pytest.raises(ValueError):
        plot.figure_data(snapshot, binding, readback)


@pytest.mark.parametrize("change", ["missing_endpoint", "wrong_f1", "missing_as_zero", "duplicate_cell",
    "extra_admission", "wrong_denominator", "wrong_fraction", "boolean_count", "failed_timing_relabelled",
    "wrong_adjustment", "unmatched_candidate", "duplicate_contrast", "unmatched_as_interval",
    "invalid_interval", "nonnested_interval", "wrong_effect", "wrong_readback", "new_draws", "ready"])
def test_scientific_and_partial_scope_guards(change):
    snapshot, binding, readback = inputs()
    effect = next(r for r in binding["contrasts"] if r["name"] == "C_at_P0_R0")
    if change == "missing_endpoint": snapshot["rows"][2]["scores"].pop("GO")
    elif change == "wrong_f1": snapshot["rows"][2]["scores"]["SwissTrees"] = .123
    elif change == "missing_as_zero": snapshot["rows"][3]["scores"]["VGNC"] = 0
    elif change == "duplicate_cell": snapshot["rows"][3]["cell"] = plot.CELLS[0]
    elif change == "extra_admission": snapshot["rows"][3]["accuracy_admitted"] = True
    elif change == "wrong_denominator": snapshot["rows"][2]["input_accessions"] -= 1
    elif change == "wrong_fraction": snapshot["rows"][2]["relation_coverage"] = .9
    elif change == "boolean_count": snapshot["rows"][2]["submitted_pairs"] = True
    elif change == "failed_timing_relabelled": snapshot["rows"][1]["timing_eligible"] = True
    elif change == "wrong_adjustment": binding["multiplicity_endpoints"] = 6
    elif change == "unmatched_candidate": effect["status"] = "native_records_unavailable"
    elif change == "duplicate_contrast": binding["contrasts"][-1] = copy.deepcopy(effect)
    elif change == "unmatched_as_interval": binding["contrasts"][-1]["metrics"] = effect["metrics"]
    elif change == "invalid_interval": effect["metrics"]["F1"]["bonferroni_percentile_ci"][0] = float("nan")
    elif change == "nonnested_interval": effect["metrics"]["F1"]["bonferroni_percentile_ci"] = [.1, .2]
    elif change == "wrong_effect": effect["metrics"]["F1"]["difference"] = .9
    elif change == "wrong_readback": readback["candidate_pair_labels_matched"] = 100
    elif change == "new_draws": readback["new_bootstrap_draws"] = 1
    else: snapshot["publication_ready"] = True
    with pytest.raises(ValueError):
        plot.figure_data(snapshot, binding, readback)


def test_existing_output_refused_before_reading_inputs(tmp_path):
    with pytest.raises(ValueError, match="Output already exists"):
        plot.run(Path("absent"), "unused", Path("absent"), "unused", Path("absent"), "unused", tmp_path, Path("absent"))


@pytest.mark.parametrize("change", ["snapshot", "readback_binding", "readback_source"])
def test_cross_input_binding_refused_before_render(tmp_path, change):
    snapshot, binding, readback = inputs()
    s = write(tmp_path / "snapshot.json", snapshot)
    if change != "snapshot": binding["snapshot"] = s
    b = write(tmp_path / "binding.json", binding)
    if change != "readback_binding": readback["binding"] = b
    if change == "readback_source": readback["source"] = {}
    r = write(tmp_path / "readback.json", readback)
    with pytest.raises(ValueError):
        plot.run(Path(s["path"]), s["sha256"], Path(b["path"]), b["sha256"], Path(r["path"]),
                 r["sha256"], tmp_path / "figure", Path("absent"))
    assert not (tmp_path / "figure").exists()


@pytest.fixture(scope="module")
def rendered(tmp_path_factory):
    root = tmp_path_factory.mktemp("three_cell_figure")
    snapshot, binding, readback = inputs()
    s = write(root / "snapshot.json", snapshot)
    binding["snapshot"] = s
    b = write(root / "binding.json", binding)
    readback["binding"] = b
    r = write(root / "readback.json", readback)
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(plot.previous, "replay_binding", lambda *args: dict(python=record(sys.executable), fixture=True))
        manifest = plot.run(Path(s["path"]), s["sha256"], Path(b["path"]), b["sha256"],
            Path(r["path"]), r["sha256"], root / "figure", Path(sys.executable))
    return root / "figure", manifest


def copy_figure(tmp_path, rendered):
    source, original = rendered
    root = tmp_path / "figure"
    shutil.copytree(source, root)
    manifest = copy.deepcopy(original)
    manifest["evidence"] = []
    manifest["outputs"] = [record(root / Path(r["path"]).name) for r in manifest["outputs"]]
    (root / "manifest.json").write_text(json.dumps(manifest))
    return root, manifest


def test_fixture_renders_and_independent_readback_decodes_actual_pdf(tmp_path, rendered):
    root, manifest = copy_figure(tmp_path, rendered)
    result = reviewer.review(root, tmp_path / "preview.png")
    assert result["score_endpoints_checked"] == 18 and result["coverage_endpoints_checked"] == 3
    assert result["contrast_endpoints_checked"] == result["precision_recall_endpoints_checked"] == 6
    assert result["png"]["size"] == [2760, 1960]
    assert all(n > 100 for n in result["png"]["color_pixels"].values())
    assert result["exact_table_readback"] is result["svg_and_decoded_pdf_required_labels_present"] is True
    assert result["automatic_visual_certification"] is result["publication_ready"] is False
    assert manifest["new_scoring_or_admission"] is manifest["scientific_timings_admitted"] is False


@pytest.mark.parametrize("change", ["score", "statistic", "duplicate_score", "missing_score", "coverage",
    "coverage_count", "semantics", "precision", "missing_coverage", "interval", "interval_name",
    "missing_interval", "svg_label", "blank_png", "pdf_content", "output_hash", "publication_ready",
    "new_bootstrap_draws", "source"])
def test_readback_refuses_resealed_content_changes(tmp_path, rendered, change):
    root, manifest = copy_figure(tmp_path, rendered)
    if change in ("score", "statistic", "duplicate_score", "missing_score"):
        path = root / "scores.tsv"
        with path.open() as stream: rows = list(csv.reader(stream, delimiter="\t"))
        if change == "score": rows[1][3] = str(float(rows[1][3]) + 1e-10)
        elif change == "statistic": rows[1][2] = "not F1"
        elif change == "duplicate_score": rows[-1] = rows[1]
        else: rows.pop()
        with path.open("w") as stream: csv.writer(stream, delimiter="\t").writerows(rows)
    elif change in ("coverage", "coverage_count", "semantics", "precision", "missing_coverage"):
        path = root / "coverage.tsv"
        with path.open() as stream: rows = list(csv.reader(stream, delimiter="\t"))
        if change == "coverage": rows[1][4] = ".9"
        elif change == "coverage_count": rows[1][1] = "100"
        elif change == "semantics": rows[1][-1] = "incorrect semantics"
        elif change == "precision": rows[1][5] = ".9"
        else: rows.pop()
        with path.open("w") as stream: csv.writer(stream, delimiter="\t").writerows(rows)
    elif change in ("interval", "interval_name", "missing_interval"):
        path = root / "swiss_intervals.tsv"
        with path.open() as stream: rows = list(csv.reader(stream, delimiter="\t"))
        if change == "interval": rows[1][4] = "1.0"
        elif change == "interval_name": rows[1][0] = "wrong_contrast"
        else: rows.pop()
        with path.open("w") as stream: csv.writer(stream, delimiter="\t").writerows(rows)
    elif change == "svg_label":
        path = root / (plot.STEM + ".svg")
        path.write_text(path.read_text().replace("Coverage is not accuracy", "Changed interpretation"))
    elif change == "blank_png": Image.new("RGB", (2760, 1960), "white").save(root / (plot.STEM + ".png"))
    elif change == "pdf_content": (root / (plot.STEM + ".pdf")).write_bytes(b"%PDF-not-valid")
    elif change == "new_bootstrap_draws": manifest[change] = 1
    elif change == "source": manifest["source"] = {}
    elif change == "publication_ready": manifest[change] = True
    manifest["outputs"] = [record(root / Path(r["path"]).name) for r in manifest["outputs"]]
    if change == "output_hash": manifest["outputs"][0]["sha256"] = "0" * 64
    (root / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises((ValueError, fitz_exception())):
        reviewer.review(root, tmp_path / "preview.png")


def fitz_exception():
    return reviewer.fitz.FileDataError


def test_existing_preview_refused_before_input_reads(tmp_path):
    preview = tmp_path / "preview.png"
    preview.write_text("retain me")
    with pytest.raises(ValueError, match="PDF preview already exists"):
        reviewer.review(Path("absent"), preview)
    assert preview.read_text() == "retain me"


def test_actual_three_cell_manifest_and_review_remain_bound():
    root = RESULTS / "native_qfo_three_cell_figure_20261006_v1"
    manifest = json.loads((root / "manifest.json").read_text())
    readback = json.loads((RESULTS / "native_qfo_three_cell_figure_readback_20261006_v1.json").read_text())
    assert manifest["source"] == record(plot.__file__)
    assert readback["source"] == record(reviewer.__file__)
    assert readback["manifest"] == record(root / "manifest.json")
    assert manifest["validation"]["exact_binding_replay"] is True
    assert manifest["validation"]["versions"]["python"] == "3.10.13"
    assert manifest["validation"]["invocation"].endswith("native_factorial_review_py310_20261004/bin/python")
    assert readback["score_endpoints_checked"] == 18
    assert readback["coverage_endpoints_checked"] == 3
    assert readback["contrast_endpoints_checked"] == 6
    assert readback["automatic_visual_certification"] is manifest["visual_review_complete"] is False
    for ref in manifest["outputs"]:
        assert record(root / Path(ref["path"]).name) == ref


def test_current_caption_and_reproduction_guide_keep_partial_scope():
    manuscript = (RESULTS / "PUBLICATION_MANUSCRIPT_DRAFT_20260916.md").read_text()
    section = manuscript.split("### Third Native QfO Cell And Candidate-Expansion Trade-Off", 1)[1].split(
        "### Native Functional-Pair Composition", 1)[0]
    assert "(native_qfo_three_cell_figure_20261006_v1/native_qfo_three_cell.png)" in section
    for phrase in ("three admitted cells and four unavailable", "Both adjusted F1 intervals include zero",
                   "coverage, not accuracy", "are not plotted as zero", "No other-endpoint uncertainty",
                   "visually inspected separately", "not estimates of isolated performance"):
        assert phrase in section
    guide = (RESULTS.parent / "PUBLICATION_REPRODUCTION.md").read_text().split("## Current Entry Point", 1)[0]
    assert "contains three admitted cells and four unavailable" in guide
    assert "plot_native_qfo_three_cell_scores" in guide and "review_native_qfo_three_cell_figure" in guide
    assert "original\nPython3.10 venv entry point" in guide
    assert "Native22444 remains live" not in guide
