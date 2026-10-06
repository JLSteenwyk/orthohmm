"""Independently read back native QfO figure tables, hashes and raster/vector assets."""

import argparse
import csv
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

import numpy as np
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.export_native_factorial_progress import require


def pixels(path):
    with Image.open(path) as picture:
        array = np.asarray(picture.convert("RGB"))
        result = dict(size=list(picture.size), nonwhite_pixels=int(np.count_nonzero(np.any(array < 240, axis=2))),
            teal_pixels=int(np.count_nonzero(np.all(array == [18, 120, 111], axis=2))),
            brick_pixels=int(np.count_nonzero(np.all(array == [182, 75, 55], axis=2))))
    require(result["nonwhite_pixels"] > 5000 and result["teal_pixels"] > 100
            and result["brick_pixels"] > 100, "Blank or missing method colors")
    return result


def review(root, pdf_preview):
    manifest_path = root / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    require(manifest["schema"] == "native_qfo_p0c0_figure_v1"
            and manifest["publication_ready"] is False
            and manifest["visual_review_complete"] is False
            and manifest["new_bootstrap_draws"] == 0
            and manifest["new_scoring_or_admission"] is False
            and manifest["scientific_timings_admitted"] is False, "Figure scope changed")
    for ref in [manifest["source"], *manifest["outputs"], *manifest["evidence"]]:
        check(ref)
    require({Path(ref["path"]).name for ref in manifest["outputs"]} == {
        "native_qfo_p0c0.png", "native_qfo_p0c0.pdf", "native_qfo_p0c0.svg", "scores.tsv", "swiss_intervals.tsv"},
        "Figure output inventory differs")
    inputs = []
    for name in ("snapshot", "swiss_binding"):
        ref = manifest[name]
        check(ref)
        inputs.append(json.loads(Path(ref["path"]).read_text()))
    snapshot, binding = inputs
    require(binding["snapshot"] == manifest["snapshot"], "Cross-snapshot figure binding")
    cells = ("p0_c0_r0", "p0_c0_r1")
    native = {row["cell"]: row for row in snapshot["rows"]}
    expected = {(cell, endpoint): value for cell in cells for endpoint, value in native[cell]["scores"].items()}
    with (root / "scores.tsv").open() as stream:
        scores = list(csv.DictReader(stream, delimiter="\t"))
    require(len(scores) == 12 and len({(r["cell"], r["endpoint"]) for r in scores}) == 12,
            "Duplicate or missing plotted endpoint")
    require({(r["cell"], r["endpoint"]): float(r["value"]) for r in scores} == expected,
            "Figure score readback differs")
    require(all(r["statistic"] == ("F1" if r["endpoint"] in ("VGNC", "SwissTrees", "TreeFam-A")
                else "similarity") for r in scores), "Figure conflates F1 and similarity")
    contrast = next(row for row in binding["contrasts"] if row["name"] == "R_at_P0_C0")
    require(contrast["status"] == "native_records_matched", "Unmatched plotted uncertainty")
    with (root / "swiss_intervals.tsv").open() as stream:
        intervals = list(csv.DictReader(stream, delimiter="\t"))
    require([r["metric"] for r in intervals] == ["F1", "PPV", "TPR"], "Interval metric inventory differs")
    for row in intervals:
        source = contrast["metrics"][row["metric"]]
        expected_interval = dict(difference_pp=100*source["difference"],
            nominal_low_pp=100*source["paired_percentile_ci"][0],
            nominal_high_pp=100*source["paired_percentile_ci"][1],
            adjusted_low_pp=100*source["bonferroni_percentile_ci"][0],
            adjusted_high_pp=100*source["bonferroni_percentile_ci"][1])
        require({k: float(v) for k, v in row.items() if k != "metric"} == expected_interval,
                "Figure interval readback differs")
    require(float(intervals[0]["adjusted_low_pp"]) < 0 < float(intervals[0]["adjusted_high_pp"]),
            "Adjusted F1 interval no longer crosses zero")
    svg = ET.parse(root / "native_qfo_p0c0.svg")
    labels = " ".join(node.text or "" for node in svg.iter() if node.tag.endswith("}text"))
    for text in ("Orthology F1", "Similarity endpoints (not F1)", "SwissTrees precision-recall",
                 "R-on minus R-off", "2/7 fresh cells admitted", "42-endpoint adjusted interval",
                 "18 development-exposed SwissTrees families", "F1 adjusted interval includes zero",
                 "Failed R-on inference timing excluded", "Not a selected-default comparison"):
        require(text in labels, "Missing figure label: " + text)
    require((root / "native_qfo_p0c0.pdf").read_bytes().startswith(b"%PDF-"), "Invalid PDF header")
    png = pixels(root / "native_qfo_p0c0.png")
    require(png["size"] == [2500, 1600], "Changed raster dimensions")
    return dict(schema="native_qfo_figure_content_readback_v1", manifest=record(manifest_path),
        snapshot=manifest["snapshot"], swiss_binding=manifest["swiss_binding"],
        source=record(__file__), checked_evidence_records=len(manifest["evidence"]),
        score_endpoints_checked=12, contrast_endpoints_checked=3,
        exact_score_and_interval_readback=True, png=png, pdf_preview=pixels(pdf_preview),
        pdf_preview_record=record(pdf_preview), svg_required_labels_present=True,
        automatic_visual_certification=False, new_scoring_or_admission=False, publication_ready=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("figure-directory", "pdf-preview", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()
    require(not args.output.exists() and not args.output.is_symlink(), "Output already exists")
    result = review(args.figure_directory, args.pdf_preview)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps(dict(score_endpoints_checked=12, contrast_endpoints_checked=3,
                         exact_score_and_interval_readback=result["exact_score_and_interval_readback"])))
