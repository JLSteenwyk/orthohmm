"""Four-cell presentation guards and synthetic-path rendering/readback."""

from copy import deepcopy
import csv
import json
from pathlib import Path
import shutil
import sys

from PIL import Image
import pytest

from benchmark_tools import plot_native_qfo_four_cell_scores as plot
from benchmark_tools import review_native_qfo_four_cell_figure as reviewer
from benchmark_tools.prepare_ob_candidate_neighborhood import record

RESULTS = Path(__file__).resolve().parents[2] / "benchmark_tools/results"


def inputs():
    return [json.loads((RESULTS / name).read_text()) for name in (
        "native_qfo_scientific_scores_20261007_v3/report.json",
        "native_qfo_four_cell_swiss_uncertainty_20261007_v1.json",
        "native_qfo_profile_swiss_readback_20261007_v1.json")]


def write(path, value):
    path.write_text(json.dumps(value, sort_keys=True))
    return record(path)


def test_four_cell_data_retains_conditional_endpoints_and_missing_cells():
    rows, scores, coverage, intervals = plot.figure_data(*inputs())
    assert len(scores) == 24 and len(coverage) == 4 and len(intervals) == 9
    assert sum(r["statistic"] == "F1" for r in scores) == 12
    assert not any("mean" in r["endpoint"].lower() for r in scores)
    assert coverage[-1]["relation_coverage"] == 542336/984137
    assert rows["p0_c0_r1"]["resources"] is None
    assert all(v is None for v in rows["p1_c1_r0"]["scores"].values())
    profile = next(r for r in intervals if r["contrast"] == "P_at_C0_R1" and r["metric"] == "F1")
    assert profile["difference_pp"] == -.32021523014126307
    assert all(r["adjusted_low_pp"] < 0 < r["adjusted_high_pp"] for r in intervals if r["metric"] == "F1")


@pytest.mark.parametrize("value", [True, float("nan"), float("inf"), "0.5", -.1, 1.1])
def test_invalid_new_profile_score_refused(value):
    snapshot, binding, readback = inputs()
    snapshot["rows"][4]["scores"]["FAS"] = value
    with pytest.raises(ValueError):
        plot.figure_data(snapshot, binding, readback)


@pytest.mark.parametrize("change", ["schema", "extra_admission", "duplicate_cell", "missing_score", "wrong_f1",
    "wrong_denominator", "wrong_fraction", "boolean_count", "failed_timing", "allocated_timing", "missing_as_zero",
    "adjustment", "new_draws", "unmatched", "extra_interval", "duplicate_contrast", "nan_interval",
    "nonnested", "wrong_effect", "wrong_reference", "wrong_readback", "changed_prior", "ready"])
def test_scientific_and_partial_scope_refused(change):
    snapshot, binding, readback = inputs()
    row = snapshot["rows"][4]
    effect = next(r for r in binding["contrasts"] if r["name"] == "P_at_C0_R1")
    if change == "schema": snapshot["schema"] = "native_qfo_scientific_reporting_snapshot_v1"
    elif change == "extra_admission": snapshot["rows"][5]["accuracy_admitted"] = True
    elif change == "duplicate_cell": snapshot["rows"][5]["cell"] = plot.CELLS[0]
    elif change == "missing_score": row["scores"].pop("GO")
    elif change == "wrong_f1": row["scores"]["SwissTrees"] = .5
    elif change == "wrong_denominator": row["input_accessions"] -= 1
    elif change == "wrong_fraction": row["relation_coverage"] = .9
    elif change == "boolean_count": row["submitted_pairs"] = True
    elif change == "failed_timing": snapshot["rows"][1]["timing_eligible"] = True
    elif change == "allocated_timing": row["scientific_timings_admitted"] = True
    elif change == "missing_as_zero": snapshot["rows"][5]["scores"]["VGNC"] = 0
    elif change == "adjustment": binding["multiplicity_endpoints"] = 9
    elif change == "new_draws": readback["new_bootstrap_draws"] = 1
    elif change == "unmatched": effect["status"] = "native_records_unavailable"
    elif change == "extra_interval": binding["contrasts"][-1]["metrics"] = effect["metrics"]
    elif change == "duplicate_contrast": binding["contrasts"][-1] = deepcopy(effect)
    elif change == "nan_interval": effect["metrics"]["F1"]["bonferroni_percentile_ci"][0] = float("nan")
    elif change == "nonnested": effect["metrics"]["F1"]["bonferroni_percentile_ci"] = [.1, .2]
    elif change == "wrong_effect": effect["metrics"]["F1"]["difference"] = .9
    elif change == "wrong_reference": effect["reference"] = plot.CELLS[0]
    elif change == "wrong_readback": readback["profile_pair_labels_matched"] = 1
    elif change == "changed_prior": readback["prior_matched_contrasts_unchanged"] = False
    else: snapshot["publication_ready"] = True
    with pytest.raises(ValueError):
        plot.figure_data(snapshot, binding, readback)


def test_existing_output_refused_before_reads(tmp_path):
    with pytest.raises(ValueError, match="already exists"):
        plot.run(Path("absent"), "unused", Path("absent"), "unused", Path("absent"), "unused", tmp_path, Path("absent"))


@pytest.mark.parametrize("change", ["snapshot", "binding", "source"])
def test_cross_input_binding_refused(tmp_path, change):
    snapshot, binding, readback = inputs()
    s = write(tmp_path / "snapshot.json", snapshot)
    if change != "snapshot": binding["snapshot"] = s
    b = write(tmp_path / "binding.json", binding)
    if change != "binding": readback["binding"] = b
    if change == "source": readback["source"] = {}
    r = write(tmp_path / "readback.json", readback)
    with pytest.raises(ValueError):
        plot.run(Path(s["path"]), s["sha256"], Path(b["path"]), b["sha256"], Path(r["path"]),
                 r["sha256"], tmp_path / "figure", Path("absent"))
    assert not (tmp_path / "figure").exists()


@pytest.fixture(scope="module")
def rendered(tmp_path_factory):
    root = tmp_path_factory.mktemp("four_cell_fixture")
    snapshot, binding, readback = inputs()
    # A presentation fixture: no claims of replaying production provenance.
    snapshot["evidence"], snapshot["outputs"] = [], []
    binding["evidence"], binding["helpers"], readback["checked_inputs"] = [], [], []
    s = write(root / "snapshot.json", snapshot)
    binding["snapshot"] = s
    b = write(root / "binding.json", binding)
    readback["binding"] = b
    r = write(root / "readback.json", readback)
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(plot, "replay_snapshot", lambda *args: dict(python=record(sys.executable),
            replay_helper=record(plot.__file__), fixture=True))
        manifest = plot.run(Path(s["path"]), s["sha256"], Path(b["path"]), b["sha256"], Path(r["path"]),
            r["sha256"], root / "figure", Path(sys.executable))
    return root / "figure", manifest


def copy_figure(tmp_path, rendered):
    source, original = rendered
    root = tmp_path / "figure"
    shutil.copytree(source, root)
    manifest = deepcopy(original)
    manifest["outputs"] = [record(root / Path(r["path"]).name) for r in manifest["outputs"]]
    write(root / "manifest.json", manifest)
    return root, manifest


def test_fixture_render_and_separate_asset_readback(tmp_path, rendered):
    root, manifest = copy_figure(tmp_path, rendered)
    result = reviewer.review(root, tmp_path / "preview.png")
    assert result["score_endpoints_checked"] == 24 and result["coverage_endpoints_checked"] == 4
    assert result["contrast_endpoints_checked"] == 9 and result["precision_recall_endpoints_checked"] == 8
    assert result["png"]["size"] == [2960, 2120]
    assert result["exact_table_readback"] and result["pdf_text_within_page_bounds"]
    assert all(n > 100 for n in result["png"]["color_pixels"].values())
    assert result["automatic_visual_certification"] is result["publication_ready"] is False
    assert manifest["new_scoring_or_admission"] is manifest["scientific_timings_admitted"] is False


@pytest.mark.parametrize("change", ["score", "statistic", "missing_score", "duplicate_score", "coverage",
    "count", "semantics", "precision", "missing_coverage", "interval", "interval_name", "missing_interval",
    "svg_label", "blank_png", "pdf", "hash", "scope", "source"])
def test_resealed_table_asset_or_scope_change_refused(tmp_path, rendered, change):
    root, manifest = copy_figure(tmp_path, rendered)
    if change in ("score", "statistic", "missing_score", "duplicate_score"):
        path = root / "scores.tsv"
        with path.open() as stream: rows = list(csv.reader(stream, delimiter="\t"))
        if change == "score": rows[-1][3] = str(float(rows[-1][3]) + 1e-10)
        elif change == "statistic": rows[-1][2] = "F1"
        elif change == "duplicate_score": rows[-1] = rows[1]
        else: rows.pop()
        with path.open("w") as stream: csv.writer(stream, delimiter="\t").writerows(rows)
    elif change in ("coverage", "count", "semantics", "precision", "missing_coverage"):
        path = root / "coverage.tsv"
        with path.open() as stream: rows = list(csv.reader(stream, delimiter="\t"))
        if change == "coverage": rows[-1][4] = ".9"
        elif change == "count": rows[-1][1] = "100"
        elif change == "semantics": rows[-1][-1] = "wrong semantics"
        elif change == "precision": rows[-1][5] = ".9"
        else: rows.pop()
        with path.open("w") as stream: csv.writer(stream, delimiter="\t").writerows(rows)
    elif change in ("interval", "interval_name", "missing_interval"):
        path = root / "swiss_intervals.tsv"
        with path.open() as stream: rows = list(csv.reader(stream, delimiter="\t"))
        if change == "interval": rows[-1][4] = "1.0"
        elif change == "interval_name": rows[-1][0] = "wrong contrast"
        else: rows.pop()
        with path.open("w") as stream: csv.writer(stream, delimiter="\t").writerows(rows)
    elif change == "svg_label":
        path = root / (plot.STEM + ".svg")
        path.write_text(path.read_text().replace("Coverage is not accuracy", "Changed interpretation"))
    elif change == "blank_png": Image.new("RGB", (2960, 2120), "white").save(root / (plot.STEM + ".png"))
    elif change == "pdf": (root / (plot.STEM + ".pdf")).write_bytes(b"%PDF-invalid")
    elif change == "scope": manifest["publication_ready"] = True
    elif change == "source": manifest["source"] = {}
    manifest["outputs"] = [record(root / Path(r["path"]).name) for r in manifest["outputs"]]
    if change == "hash": manifest["outputs"][0]["sha256"] = "0" * 64
    write(root / "manifest.json", manifest)
    with pytest.raises((ValueError, reviewer.fitz.FileDataError)):
        reviewer.review(root, tmp_path / "preview.png")


def test_existing_preview_refused(tmp_path):
    preview = tmp_path / "retain.png"
    preview.write_text("retain\n")
    with pytest.raises(ValueError, match="already exists"):
        reviewer.review(Path("absent"), preview)
    assert preview.read_text() == "retain\n"
