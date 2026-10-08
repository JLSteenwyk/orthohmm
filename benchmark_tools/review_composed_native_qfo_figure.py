"""Independently check composed figure tables, scope and decoded assets."""

import argparse
import csv
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import fitz
import numpy as np
from PIL import Image

from benchmark_tools.export_native_factorial_progress import require
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

STEM = "composed_native_qfo"
CELLS = ("p0_c0_r0", "p0_c0_r1", "p0_c1_r0", "p1_c0_r1", "p1_c1_r1")
LABELS = ("Orthology F1", "Similarities (not F1)", "All-input relation coverage",
    "Conditional SwissTrees contrasts", "5/7 native cells admitted", "two scores unavailable",
    "Initial HMM search remains on", "42-endpoint adjusted intervals",
    "18 development-exposed SwissTrees families", "Coverage is not accuracy",
    "secondary mean is not plotted", "Failed timing remains excluded", "No isolated timing")


def table(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def pixels(path):
    with Image.open(path) as picture:
        data = np.asarray(picture.convert("RGB"))
        result = dict(size=list(picture.size), nonwhite_pixels=int(np.count_nonzero(np.any(data < 240, axis=2))),
            color_pixels=[int(np.count_nonzero(np.all(data == rgb, axis=2))) for rgb in (
                [18, 120, 111], [182, 75, 55], [70, 83, 166], [154, 108, 0], [125, 67, 141])])
    require(result["nonwhite_pixels"] > 5000 and all(count > 100 for count in result["color_pixels"]),
        "Blank asset or missing admitted-cell color")
    return result


def review(root, preview):
    require(not preview.exists() and not preview.is_symlink(), "Preview already exists")
    manifest = json.loads((root / "manifest.json").read_text())
    require(manifest["schema"] == "composed_native_qfo_figure_v1" and manifest["plotted_cells"] == list(CELLS)
        and manifest["plotted_score_endpoints"] == 30 and manifest["plotted_coverage_endpoints"] == 5
        and manifest["table_score_rows"] == 42 and manifest["table_coverage_rows"] == 7
        and manifest["table_status_rows"] == 7 and type(manifest["new_bootstrap_draws"]) is int
        and manifest["new_bootstrap_draws"] == 0
        and manifest["consumed_native_bootstrap_draws"] == 100000
        and all(manifest[key] is False for key in ("new_scoring_or_admission", "scientific_timings_admitted",
            "independent_confirmation", "publication_ready", "visual_review_complete")), "Changed figure scope")
    require(manifest["source"] == record(Path(__file__).with_name("plot_composed_native_qfo_scores.py")),
        "Changed figure source")
    expected_names = {f"{STEM}.{suffix}" for suffix in ("png", "pdf", "svg")}
    expected_names.update(("scores.tsv", "coverage.tsv", "status.tsv", "swiss_intervals.tsv"))
    require(len(manifest["outputs"]) == len(expected_names)
        and {Path(ref["path"]).name for ref in manifest["outputs"]} == expected_names
        and all(record(root / Path(ref["path"]).name) == ref for ref in manifest["outputs"]),
        "Changed asset/table inventory")
    evidence = [manifest["source"], *manifest["outputs"], *manifest["evidence"],
        manifest["snapshot"], manifest["native_intervals"]]
    for ref in evidence: check(ref)
    snapshot, uncertainty = [json.loads(Path(manifest[key]["path"]).read_text())
        for key in ("snapshot", "native_intervals")]
    require(snapshot["schema"] == "composed_native_qfo_scientific_reporting_v1"
        and snapshot["source"] == record(Path(__file__).with_name("export_composed_native_qfo_scientific_scores.py"))
        and uncertainty["schema"] == "composed_native_qfo_swiss_bootstrap_v1"
        and uncertainty["source"] == record(Path(__file__).with_name("bootstrap_composed_native_qfo_swiss.py"))
        and uncertainty["snapshot"] == manifest["snapshot"], "Changed source or cross-snapshot interval")
    rows = {row["cell"]: row for row in snapshot["rows"]}
    require(len(rows) == 7 and {cell for cell, row in rows.items() if row["accuracy_admitted"]} == set(CELLS)
        and manifest["unavailable_score_cells"] == [row["cell"] for row in snapshot["rows"] if not row["accuracy_admitted"]]
        and rows["p1_c1_r0"]["scoring_status"] == "OUT_OF_MEMORY"
        and rows["p0_c0_r1"]["resources"] is None and rows["p0_c0_r1"]["timing_eligible"] is False,
        "Changed missing/failure/timing scope")
    scores = table(root / "scores.tsv")
    expected = {(cell, endpoint): value for cell, row in rows.items() for endpoint, value in row["scores"].items()}
    observed = {(row["cell"], row["endpoint"]): float(row["value"]) if row["value"] else None for row in scores}
    require(len(scores) == len(observed) == 42 and observed == expected
        and all(row["statistic"] == ("F1" if row["endpoint"] in ("VGNC", "SwissTrees", "TreeFam-A")
            else "similarity") and row["accuracy_admitted"] == str(rows[row["cell"]]["accuracy_admitted"])
            for row in scores), "Changed score, statistic or null outcome")
    coverage, statuses = table(root / "coverage.tsv"), table(root / "status.tsv")
    require(len(coverage) == len(statuses) == 7
        and [row["cell"] for row in coverage] == [row["cell"] for row in snapshot["rows"]]
        and [row["cell"] for row in statuses] == [row["cell"] for row in snapshot["rows"]],
        "Changed coverage/status inventory")
    for row, status in zip(coverage, statuses):
        original = rows[row["cell"]]
        for key in ("input_accessions", "relation_accessions", "submitted_pairs"):
            require((int(row[key]) if row[key] else None) == original[key], "Changed coverage count")
        require((float(row["relation_coverage"]) if row["relation_coverage"] else None) == original["relation_coverage"]
            and row["prediction_semantics"] == original["prediction_semantics"]
            and row["accuracy_admitted"] == status["accuracy_admitted"] == str(original["accuracy_admitted"])
            and int(status["index"]) == original["index"] and status["status"] == original["status"]
            and status["scoring_status"] == original.get("scoring_status", ""), "Changed coverage or failed status")
    require(uncertainty["observed_cells"] == list(CELLS) and uncertainty["replicates"] == 100000
        and uncertainty["seed"] == 20260922 and uncertainty["multiplicity_endpoints"] == 42
        and uncertainty["retained_intervals_reused"] is False and uncertainty["unobserved_cells_imputed"] is False,
        "Changed interval scope")
    expected_intervals = []
    for effect in uncertainty["comparisons"]:
        if effect["metrics"] is None: continue
        for metric in ("F1", "PPV", "TPR"):
            value = effect["metrics"][metric]
            expected_intervals.append(dict(contrast=effect["name"], metric=metric,
                difference_pp=100*value["difference"], nominal_low_pp=100*value["paired_percentile_ci"][0],
                nominal_high_pp=100*value["paired_percentile_ci"][1],
                adjusted_low_pp=100*value["bonferroni_percentile_ci"][0],
                adjusted_high_pp=100*value["bonferroni_percentile_ci"][1]))
    actual = table(root / "swiss_intervals.tsv")
    require(len(actual) == len(expected_intervals) == manifest["plotted_swiss_contrast_endpoints"],
        "Changed interval table inventory")
    require([dict(row, **{key: float(value) for key, value in row.items()
        if key not in ("contrast", "metric")}) for row in actual] == expected_intervals, "Changed interval table value")
    svg = ET.parse(root / f"{STEM}.svg")
    text = " ".join(node.text or "" for node in svg.iter() if node.tag.endswith("}text"))
    require(all(label in text for label in LABELS), "Missing SVG scope label")
    png = pixels(root / f"{STEM}.png")
    require(png["size"] == [3000, 2400], "Changed raster dimensions")
    with fitz.open(root / f"{STEM}.pdf") as document:
        require(len(document) == 1, "Changed PDF page count")
        page = document[0]
        text = " ".join(page.get_text().split())
        require(all(label in text for label in LABELS), "Missing PDF scope label")
        spans = [span for block in page.get_text("dict")["blocks"] if "lines" in block
            for line in block["lines"] for span in line["spans"]]
        require(all(span["bbox"][0] >= -.5 and span["bbox"][1] >= -.5
            and span["bbox"][2] <= page.rect.width+.5 and span["bbox"][3] <= page.rect.height+.5
            for span in spans), "PDF text outside page bounds")
        page.get_pixmap(matrix=fitz.Matrix(120/72, 120/72)).save(preview)
    decoded = pixels(preview)
    return dict(schema="composed_native_qfo_figure_readback_v1", source=record(__file__),
        manifest=record(root / "manifest.json"), snapshot=manifest["snapshot"], native_intervals=manifest["native_intervals"],
        score_rows=42, coverage_rows=7, status_rows=7, interval_rows=len(actual),
        png=png, pdf_preview=record(preview), decoded_pdf=decoded, checked_inputs=evidence,
        scientific_tables_matched=True, assets_decoded=True, visual_review_complete=False,
        new_scoring_or_admission=False, independent_confirmation=False, publication_ready=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("figure-root", "pdf-preview", "output"):
        parser.add_argument("--"+name, type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists() and not args.output.is_symlink(), "Output already exists")
    result = review(args.figure_root, args.pdf_preview)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps(dict(score_rows=42, interval_rows=result["interval_rows"], assets_decoded=True)))


if __name__ == "__main__":
    main()
