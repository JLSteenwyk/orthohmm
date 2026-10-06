"""Independent figure-table readback and raster/vector checks, not human review."""

import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

import numpy as np
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.export_native_factorial_progress import require
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

CELLS = ("p0_c0_r0", "p0_c0_r1")
BINS = ("all", "lower_entropy", "higher_entropy", "short_relative", "not_short_relative")


def raster(path):
    with Image.open(path) as picture:
        pixels = np.asarray(picture.convert("RGB"))
        colors = {name: int(np.count_nonzero(np.all(pixels == color, axis=2))) for name, color in
                  (("r0_teal", [18, 120, 111]), ("r1_brick", [182, 75, 55]),
                   ("bidirectional_blue", [38, 133, 154]), ("same_root_gray", [87, 103, 116]))}
        require(all(value > 100 for value in colors.values()), "Blank or missing panel colors")
        height, width = pixels.shape[:2]
        panels = [int(np.count_nonzero(np.any(pixels[y0:y1, x0:x1] < 240, axis=2)))
                  for y0, y1 in ((int(.18 * height), int(.50 * height)), (int(.60 * height), int(.86 * height)))
                  for x0, x1 in ((int(.15 * width), int(.48 * width)), (int(.68 * width), int(.98 * width)))]
        require(all(value > 1000 for value in panels), "Blank panel crop")
        return dict(size=list(picture.size), colors=colors, panel_nonwhite_pixels=panels)


def review(root, preview):
    manifest_ref = record(root / "manifest.json")
    manifest = json.loads((root / "manifest.json").read_text())
    require(manifest["schema"] == "native_qfo_swiss_mechanism_figure_v1" and manifest["panels"] == 4
            and manifest["source"] == record(Path(__file__).with_name("plot_native_qfo_swiss_mechanism.py"))
            and all(manifest[k] is False for k in ("scientific_evidence_replayed", "raw_search_rescanned", "new_uncertainty",
                "new_scoring_or_admission", "scientific_timings_admitted", "publication_ready", "visual_review_complete")),
            "Wrong figure source/scope")
    refs = [manifest_ref, manifest["source"], *manifest["evidence"], *manifest["outputs"]]
    require(len(manifest["outputs"]) == 6 and {Path(r["path"]).name for r in manifest["outputs"]} == {
        "native_swiss_mechanism.png", "native_swiss_mechanism.pdf", "native_swiss_mechanism.svg",
        "stages.tsv", "strata.tsv", "differences.tsv"}, "Wrong figure artifact inventory")
    for ref in refs:
        check(ref)
    source = {}
    for kind, paired in manifest["inputs"].items():
        report_ref, readback_ref = paired
        check(report_ref)
        check(readback_ref)
        refs.extend(paired)
        report = json.loads(Path(report_ref["path"]).read_text())
        readback = json.loads(Path(readback_ref["path"]).read_text())
        require(readback["report"] == report_ref, "Wrong diagnostic readback binding")
        source[kind] = report
    require(set(source) == {"search", "reconciliation", "strata"}, "Wrong diagnostic inventory")
    expected = []
    for label in ("TP", "FP"):
        selected = [r for r in source["search"]["cases"] if r["before"] == label]
        count = Counter(r["direct_search"][CELLS[0]]["support"] for r in selected)
        for kind in ("no_direct_hit", "one_direction", "both_directions"):
            expected.append(dict(panel="A", label=label, category=kind, count=count[kind], total=len(selected),
                                 percentage=100 * count[kind] / len(selected)))
        root_counts = Counter(r["same_root_hog"] for r in selected)
        for flag in (True, False):
            expected.append(dict(panel="B", label=label, category="same_root" if flag else "different_root",
                count=root_counts[flag], total=len(selected), percentage=100 * root_counts[flag] / len(selected)))
    originals = {(r["cell"], r["stratum"]): r for r in source["strata"]["rows"]}
    score_rows = [{k: originals[cell, name][k] for k in ("cell", "stratum", "families", "F1", "PPV", "TPR")}
                  for name in BINS for cell in CELLS]
    differences = [dict(stratum=name, metric=metric, difference_pp=100 *
                   (originals[CELLS[1], name][metric] - originals[CELLS[0], name][metric]))
                   for name in BINS for metric in ("PPV", "TPR")]
    for name, wanted in (("stages.tsv", expected), ("strata.tsv", score_rows), ("differences.tsv", differences)):
        with (root / name).open(newline="") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            require(reader.fieldnames == list(wanted[0]), "Wrong plotted-table columns")
            observed = list(reader)
        require(observed == [{k: str(v) for k, v in row.items()} for row in wanted], "Plotted-table values differ")
    require(manifest["changed_pairs"] == len(source["search"]["cases"]) == 2023
            and all(manifest[k] == 10 for k in ("stage_rows", "strata_rows", "difference_rows")), "Wrong plotted counts")
    svg = ET.parse(root / "native_swiss_mechanism.svg")
    labels = " ".join(node.text or "" for node in svg.iter() if node.tag.endswith("}text"))
    for text in ("Direct significant search support", "Final membership of excluded pairs", "Frozen sequence-strata F1",
                 "Precision-recall trade-offs", "No direct hit", "One direction", "Both directions", "Same RootHOG",
                 "Different RootHOG", "Lower entropy (9)", "Higher entropy (9)", "Short-relative (7)",
                 "R-off: group-clique", "R-on: resolved pairs", "no subgroup intervals or significance claims",
                 "homology support is not orthology confidence", "Failed R-on native timing remains ineligible"):
        require(text in labels, "Missing SVG scope/data label: " + text)
    original_pixels, preview_pixels = raster(root / "native_swiss_mechanism.png"), raster(preview)
    require(original_pixels["size"] == [2640, 1760], "Unexpected original raster dimensions")
    preview_ref = record(preview)
    for ref in refs:
        check(ref)
    return dict(schema="native_qfo_swiss_mechanism_figure_readback_v1", manifest=manifest_ref, source=record(__file__),
        checked_inputs=refs, pdf_preview=preview_ref, original_raster=original_pixels, pdf_raster=preview_pixels,
        plotted_table_rows_checked=30, panels_checked=4, new_uncertainty=False, new_scoring_or_admission=False,
        scientific_timings_admitted=False, publication_ready=False, human_visual_review_complete=False,
        limitations=["Exact table/report readback and image/label checks, not mathematical validation of raster geometry.",
                     "Human PNG/PDF layout inspection must be recorded separately; automatic checks do not prove nonoverlap.",
                     "Retained scientific diagnostic receipts reused, not raw inference/search/tree/scoring replay."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "pdf-preview", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists() and not args.output.is_symlink(), "Output already exists")
    result = review(args.root, args.pdf_preview)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({k: result[k] for k in ("plotted_table_rows_checked", "panels_checked")}))
