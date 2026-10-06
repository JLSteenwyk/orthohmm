"""Actual figure readback, bounded panel pixels and deliberate corruption refusal."""

import json
from pathlib import Path
import shutil

import pytest
from PIL import Image

from benchmark_tools import review_native_qfo_swiss_mechanism as reviewer

RESULTS = Path(reviewer.__file__).parent / "results"
ROOT = RESULTS / "native_qfo_swiss_mechanism_figure_20261006_v3"
PREVIEW = RESULTS / "native_qfo_swiss_mechanism_pdf_preview_20261006_v3.png"


def test_actual_complete_figure_readback():
    result = reviewer.review(ROOT, PREVIEW)
    assert result["plotted_table_rows_checked"] == 30 and result["panels_checked"] == 4
    assert result["original_raster"]["size"] == [2640, 1760]
    assert all(value > 100 for value in result["original_raster"]["colors"].values())
    assert all(value > 1000 for value in result["pdf_raster"]["panel_nonwhite_pixels"])
    assert result["human_visual_review_complete"] is result["publication_ready"] is False


def test_actual_readback_receipt_matches_current_sources():
    receipt = json.loads((RESULTS / "native_qfo_swiss_mechanism_figure_readback_20261006_v3.json").read_text())
    assert receipt["source"] == reviewer.record(reviewer.__file__)
    assert receipt["manifest"] == reviewer.record(ROOT / "manifest.json")
    assert receipt["pdf_preview"] == reviewer.record(PREVIEW)
    assert receipt["plotted_table_rows_checked"] == 30 and receipt["panels_checked"] == 4
    assert receipt["scientific_timings_admitted"] is receipt["publication_ready"] is False


def test_layout_corrections_preserve_original_plotted_bytes():
    versions = [RESULTS / ("native_qfo_swiss_mechanism_figure_20261006_v" + str(v)) for v in (1, 2, 3)]
    for name in ("stages.tsv", "strata.tsv", "differences.tsv"):
        assert len({(directory / name).read_bytes() for directory in versions}) == 1


@pytest.mark.parametrize("fault", ("scope", "source", "count", "stage_value", "stratum_value", "difference",
                                 "table_header", "svg_label", "png_blank"))
def test_corrupted_actual_figure_refuses(tmp_path, fault):
    directory = tmp_path / "figure"
    shutil.copytree(ROOT, directory)
    manifest = json.loads((directory / "manifest.json").read_text())
    manifest["outputs"] = [reviewer.record(directory / Path(r["path"]).name) for r in manifest["outputs"]]
    if fault == "scope":
        manifest["scientific_timings_admitted"] = True
    elif fault == "source":
        manifest["source"] = reviewer.record(__file__)
    elif fault == "count":
        manifest["changed_pairs"] += 1
    elif fault == "png_blank":
        Image.new("RGB", (2640, 1760), "white").save(directory / "native_swiss_mechanism.png")
    else:
        name = {"stage_value": "stages.tsv", "stratum_value": "strata.tsv", "difference": "differences.tsv",
                "table_header": "strata.tsv", "svg_label": "native_swiss_mechanism.svg"}[fault]
        path = directory / name
        data = path.read_text()
        if fault == "stage_value":
            data = data.replace("\t60\t334\t", "\t61\t334\t", 1)
        elif fault == "stratum_value":
            data = data.replace("0.689183963152562", "0.1", 1)
        elif fault == "difference":
            data = data.replace("30.52118719287774", "0", 1)
        elif fault == "table_header":
            data = data.replace("PPV", "wrong", 1)
        elif fault == "svg_label":
            data = data.replace("No direct hit", "wrong", 1)
        path.write_text(data)
    manifest["outputs"] = [reviewer.record(directory / Path(r["path"]).name) for r in manifest["outputs"]]
    (directory / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError):
        reviewer.review(directory, PREVIEW)


def test_blank_image_refused(tmp_path):
    path = tmp_path / "blank.png"
    Image.new("RGB", (100, 100), "white").save(path)
    with pytest.raises(ValueError, match="colors"):
        reviewer.raster(path)
