"""Fixed selection, descriptive scope and fixture rendering of stratum panels."""

from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools import plot_native_qfo_three_cell_strata as plot

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "benchmark_tools/results"


def inputs():
    return (json.loads((RESULTS / name).read_text()) for name, _ in plot.PINS.values())


def test_all_nonempty_nonredundant_bins_and_exact_values_are_selected():
    report, readback = inputs()
    points, excluded = plot.figure_data(report, readback)
    assert len(points) == 66 and len(excluded) == 9
    assert sum(r["reason"] == "empty_bin" for r in excluded) == 5
    assert sum(r["reason"] == "redundant_all_family_copy" for r in excluded) == 4
    assert {r["families"] for r in points} == {3, 6, 7, 9, 11, 12, 15, 18}
    for point in points:
        row = next(r for r in report["differences"] if all(
            point[k] == r[k] for k in ("suite", "stratum", "contrast")))
        assert point["difference_pp"] == 100 * row[point["metric"]]
    assert "not confidence intervals" in plot.CAPTION


@pytest.mark.parametrize("change", ["scope", "missing", "duplicate", "nan", "unit", "unplotted", "empty"])
def test_invalid_selection_and_scope_rejected(change):
    report, readback = inputs()
    report = deepcopy(report)
    if change == "scope":
        report["new_uncertainty"] = True
    elif change == "missing":
        report["differences"].pop()
    elif change == "duplicate":
        report["differences"][1] = deepcopy(report["differences"][0])
    elif change == "nan":
        next(r for r in report["differences"] if r["stratum"] == "higher_entropy")["PPV"] = float("nan")
    elif change == "unit":
        report["differences"][0]["reference"] = "wrong"
    elif change == "unplotted":
        report["bins"]["sequence"]["missing"] = ["APP"]
    else:
        readback["score_rows_checked"] = 0
    with pytest.raises((ValueError, KeyError)):
        plot.figure_data(report, readback)


def test_fixture_render_is_nonblank_and_has_all_formats(tmp_path):
    from PIL import Image
    import fitz
    report, readback = inputs()
    # Render invented fixture differences, not the selected numeric figure.
    for row in report["differences"]:
        for k in plot.METRICS:
            if row[k] is not None:
                row[k] = .01 if k == "TPR" else -.01
    points, _ = plot.figure_data(report, readback)
    plot.render(points, tmp_path)
    with Image.open(tmp_path / "fixed_stratum_contrasts.png") as image:
        assert image.size == (2240, 1280)
        assert len(image.convert("RGB").getcolors(image.width * image.height)) > 50
    with fitz.open(tmp_path / "fixed_stratum_contrasts.pdf") as doc:
        assert len(doc) == 1
        text = doc[0].get_text()
        assert "Candidate expansion minus baseline" in text
        assert "not confidence intervals" in text
    assert (tmp_path / "fixed_stratum_contrasts.svg").stat().st_size > 10000


def test_existing_namespace_refused_before_loading(tmp_path):
    with pytest.raises(ValueError, match="already exists"):
        plot.build(tmp_path / "missing_repo", tmp_path)


def test_actual_figure_current_sources_inputs_outputs_and_plotted_tsv():
    import csv
    root = RESULTS / "native_qfo_three_cell_strata_figure_20261007_v1"
    manifest = json.loads((root / "manifest.json").read_text())
    assert plot.record(root / "manifest.json")["sha256"] == (
        "e0157657d24ed39da8f8e39a41590311ec480537ac5aa35e06483d758ef85b02"
    )
    assert manifest["source"] == plot.record(plot.__file__)
    for ref in [*manifest["checked_inputs"], *manifest["outputs"]]:
        assert plot.record(ref["path"]) == ref
    report, readback = inputs()
    points, excluded = plot.figure_data(report, readback)
    assert manifest["plotted_rows"] == points and manifest["excluded_bins"] == excluded
    with (root / "plotted_rows.tsv").open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        assert reader.fieldnames == list(plot.FIELDS)
        rows = list(reader)
    assert len(rows) == 66
    for row, point in zip(rows, points):
        assert all(row[k] == str(point[k]) for k in plot.FIELDS)
    assert manifest["new_uncertainty"] is manifest["publication_ready"] is False


def test_actual_png_pdf_svg_labels_pixels_and_text_bounds():
    from PIL import Image
    import fitz
    import xml.etree.ElementTree as ET
    root = RESULTS / "native_qfo_three_cell_strata_figure_20261007_v1"
    for name in ("fixed_stratum_contrasts.png", "pdf_page_1.png"):
        with Image.open(root / name) as image:
            assert image.size == (2240, 1280)
            rgb = image.convert("RGB")
            for box in ((672, 256, 1320, 1024), (1480, 256, 2190, 1024)):
                colors = rgb.crop(box).getcolors(image.width * image.height)
                assert len(colors) > 50
                assert sum(count for count, color in colors if min(color) < 150) > 2000
    with fitz.open(root / "fixed_stratum_contrasts.pdf") as doc:
        assert len(doc) == 1
        page = doc[0]
        text = page.get_text()
        assert "not confidence intervals" in text
        assert "Axes use different ranges" in text
        assert "Reconciliation minus baseline" in text and "Candidate expansion minus baseline" in text
        report, readback = inputs()
        points, _ = plot.figure_data(report, readback)
        for p in points[:33:3]:
            assert p["label"] in text
        for block in page.get_text("dict")["blocks"]:
            for line in block.get("lines", []):
                for span in line.get("spans", []):
                    x0, y0, x1, y1 = span["bbox"]
                    assert 0 <= x0 <= x1 <= page.rect.width
                    assert 0 <= y0 <= y1 <= page.rect.height
    svg = ET.parse(root / "fixed_stratum_contrasts.svg")
    words = " ".join("".join(e.itertext()) for e in svg.iter() if e.tag.endswith("}text"))
    assert "Native SwissTrees fixed-stratum contrasts" in words
    assert "not confidence intervals" in words


def test_manual_review_and_source_commit_are_bound():
    import subprocess
    review = json.loads((RESULTS / "native_qfo_three_cell_strata_figure_review_20261007_v1.json").read_text())
    assert plot.record(ROOT / review["manifest_path"])["sha256"] == review["manifest_sha256"]
    ref = review["pdf_rasterization"]
    current = plot.record(ROOT / ref["path"])
    assert current["bytes"] == ref["bytes"] and current["sha256"] == ref["sha256"]
    assert review["full_png_actually_viewed"] is review["full_pdf_page_raster_actually_viewed"] is True
    assert review["old_figures_or_archives_rebuilt"] is review["publication_ready"] is False
    source = subprocess.run(["git", "show", review["source_commit_before_render"] +
                             ":benchmark_tools/plot_native_qfo_three_cell_strata.py"],
                            cwd=ROOT, capture_output=True, check=True)
    assert source.stdout == Path(plot.__file__).read_bytes()
