import csv
import json
from pathlib import Path
import shutil

from PIL import Image
import pytest

from benchmark_tools import review_native_qfo_figure as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "benchmark_tools/results"
FIGURE = RESULTS / "native_qfo_p0c0_figure_20261006_v1"


def fixture(tmp_path):
    root = tmp_path / "figure"
    shutil.copytree(FIGURE, root)
    manifest = json.loads((root / "manifest.json").read_text())
    manifest["evidence"] = []
    manifest["source"] = record(ROOT / "benchmark_tools/plot_native_qfo_scientific_scores.py")
    manifest["snapshot"] = record(RESULTS / "native_qfo_scientific_scores_20261006_v1/report.json")
    binding = json.loads((RESULTS / "native_qfo_swiss_uncertainty_binding_22449_20261006.json").read_text())
    binding["snapshot"] = manifest["snapshot"]
    (tmp_path / "binding.json").write_text(json.dumps(binding))
    manifest["swiss_binding"] = record(tmp_path / "binding.json")
    manifest["outputs"] = [record(root / Path(ref["path"]).name) for ref in manifest["outputs"]]
    (root / "manifest.json").write_text(json.dumps(manifest))
    return root, manifest


def test_readback_fixture_checks_tables_assets_without_new_science(tmp_path):
    root, _ = fixture(tmp_path)
    result = module.review(root, RESULTS / "native_qfo_p0c0_pdf_preview_20261006.png")
    assert result["exact_score_and_interval_readback"] is True
    assert result["score_endpoints_checked"] == 12
    assert result["contrast_endpoints_checked"] == 3
    assert result["png"]["size"] == [2500, 1600]
    assert result["publication_ready"] is result["automatic_visual_certification"] is False
    assert result["new_scoring_or_admission"] is False


@pytest.mark.parametrize("change", ["score", "statistic", "duplicate_score", "missing_score",
    "interval", "missing_interval", "missing_label", "blank_png", "pdf_header", "wrong_hash",
    "publication_ready", "visual_review_complete", "new_scoring_or_admission", "new_bootstrap_draws"])
def test_readback_refuses_changed_content_even_with_refreshed_output_hashes(tmp_path, change):
    root, manifest = fixture(tmp_path)
    if change in ("score", "statistic", "duplicate_score", "missing_score"):
        path = root / "scores.tsv"
        with path.open() as stream:
            rows = list(csv.reader(stream, delimiter="\t"))
        if change == "score": rows[1][3] = str(float(rows[1][3]) + 1e-10)
        elif change == "statistic": rows[1][2] = "incorrect statistic"
        elif change == "duplicate_score": rows[-1] = rows[1]
        else: rows.pop()
        with path.open("w") as stream:
            csv.writer(stream, delimiter="\t", lineterminator="\n").writerows(rows)
    elif change in ("interval", "missing_interval"):
        path = root / "swiss_intervals.tsv"
        with path.open() as stream:
            rows = list(csv.reader(stream, delimiter="\t"))
        if change == "interval": rows[1][4] = "1.0"
        else: rows.pop()
        with path.open("w") as stream:
            csv.writer(stream, delimiter="\t", lineterminator="\n").writerows(rows)
    elif change == "missing_label":
        path = root / "native_qfo_p0c0.svg"
        path.write_text(path.read_text().replace("Orthology F1", "Incorrect label"))
    elif change == "blank_png":
        Image.new("RGB", (2500, 1600), "white").save(root / "native_qfo_p0c0.png")
    elif change == "pdf_header": (root / "native_qfo_p0c0.pdf").write_bytes(b"not a PDF")
    elif change == "new_bootstrap_draws": manifest[change] = 1
    elif change != "wrong_hash": manifest[change] = True
    manifest["outputs"] = [record(ref["path"]) for ref in manifest["outputs"]]
    if change == "wrong_hash": manifest["outputs"][0]["sha256"] = "0" * 64
    (root / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError):
        module.review(root, RESULTS / "native_qfo_p0c0_pdf_preview_20261006.png")


def test_actual_committed_figure_and_readback_hashes_remain_exact():
    manifest = json.loads((FIGURE / "manifest.json").read_text())
    review = json.loads((RESULTS / "native_qfo_p0c0_figure_readback_20261006.json").read_text())
    for ref in manifest["outputs"]:
        actual = record(FIGURE / Path(ref["path"]).name)
        assert (actual["sha256"], actual["bytes"]) == (ref["sha256"], ref["bytes"])
    assert record(FIGURE / "manifest.json")["sha256"] == review["manifest"]["sha256"]
    assert record(RESULTS / "native_qfo_p0c0_pdf_preview_20261006.png")["sha256"] == review["pdf_preview_record"]["sha256"]
    assert record(ROOT / "benchmark_tools/review_native_qfo_figure.py")["sha256"] == review["source"]["sha256"]
    assert review["checked_evidence_records"] == 114
    assert review["automatic_visual_certification"] is manifest["visual_review_complete"] is False
    assert manifest["validation"]["exact_binding_replay"] is True
    assert manifest["validation"]["versions"]["python"] == "3.10.13"
    assert manifest["validation"]["invocation"].endswith("native_factorial_review_py310_20261004/bin/python")


def test_current_manuscript_caption_keeps_partial_scope_and_actual_numbers():
    manuscript = (RESULTS / "PUBLICATION_MANUSCRIPT_DRAFT_20260916.md").read_text()
    section = manuscript.split("### Second Native QfO Cell And Conditional Reconciliation Contrast", 1)[1].split(
        "The original-release QfO factorial", 1)[0]
    snapshot = json.loads((RESULTS / "native_qfo_scientific_scores_20261006_v1/report.json").read_text())
    for row in snapshot["rows"][:2]:
        for endpoint, value in row["scores"].items():
            assert format(value, ".8f" if endpoint in ("GO", "EC") else ".10f") in section
    for text in ("two admitted cells and five", "timing remains failed and ineligible",
                 "not a uniform endpoint improvement", "not a non-HMM control",
                 "not official QfO", "10.0390", "[-4.6418, 24.4199]", "[14.3252, 46.5473]",
                 "[-18.0413, -0.0436]", "all 42", "18 development-exposed families",
                 "pair-IID SEM is not a", "unknown, tool-dependent", "does not establish isolated performance"):
        assert text in section
    assert "(native_qfo_p0c0_figure_20261006_v1/native_qfo_p0c0.png)" in section
    assert "(NATIVE_QFO_FIGURE_RESULT_20261006.md)" in section
    claims = (RESULTS / "PUBLICATION_CLAIMS_20260916.md").read_text()
    assert "Two of seven fresh cells are admitted; five remain unavailable" in claims
    assert "F1 +10.0390 pp, adjusted interval [-4.6418, 24.4199] includes zero" in claims
    assert "| The package is publication-ready | All sections below | Not achieved |" in claims
    guide = (ROOT / "benchmark_tools/PUBLICATION_REPRODUCTION.md").read_text()
    assert "--validation-python" in guide and "original Python3.10 venv entry point" in guide
    assert "not portable fresh-install" in guide
