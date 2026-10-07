"""Read back three-cell QfO figure tables and independently decode its assets."""

import argparse
import csv
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

import fitz
import numpy as np
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.export_native_factorial_progress import require
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

CELLS = ("p0_c0_r0", "p0_c0_r1", "p0_c1_r0")
CONTRASTS = ("C_at_P0_R0", "R_at_P0_C0")
STEM = "native_qfo_three_cell"
LABELS = ("Orthology F1", "Similarities (not F1)", "All-input relation coverage",
    "SwissTrees: candidate expansion", "SwissTrees: reconciliation", "SwissTrees precision-recall",
    "3/7 fresh cells admitted", "four scores unavailable", "Initial HMM search on",
    "42-endpoint adjusted interval", "Both adjusted F1 intervals include zero",
    "18 development-exposed SwissTrees families", "Coverage is not accuracy",
    "secondary mean is not plotted", "Failed R-on timing remains excluded", "No timing comparison")


def pixels(path):
    with Image.open(path) as picture:
        data = np.asarray(picture.convert("RGB"))
        result = dict(size=list(picture.size), nonwhite_pixels=int(np.count_nonzero(np.any(data < 240, axis=2))),
            color_pixels={name: int(np.count_nonzero(np.all(data == rgb, axis=2)))
                for name, rgb in (("baseline", [18, 120, 111]), ("reconciliation", [182, 75, 55]),
                                  ("candidate", [70, 83, 166]))})
    require(result["nonwhite_pixels"] > 5000 and all(n > 100 for n in result["color_pixels"].values()),
            "Blank asset or missing cell color")
    return result


def table(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def review(root, pdf_preview):
    require(not pdf_preview.exists() and not pdf_preview.is_symlink(), "PDF preview already exists")
    manifest_path = root / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    require(manifest["schema"] == "native_qfo_three_cell_figure_v1"
            and manifest["plotted_cells"] == list(CELLS)
            and manifest["admitted_cells_in_snapshot"] == 3
            and manifest["plotted_score_endpoints"] == 18
            and manifest["plotted_coverage_endpoints"] == 3
            and manifest["plotted_swiss_contrast_endpoints"] == 6
            and type(manifest["new_bootstrap_draws"]) is int and manifest["new_bootstrap_draws"] == 0
            and all(manifest[k] is False for k in ("publication_ready", "visual_review_complete",
                "new_scoring_or_admission", "scientific_timings_admitted")), "Changed figure scope")
    require(manifest["source"] == record(Path(__file__).with_name("plot_native_qfo_three_cell_scores.py")),
            "Changed figure source")
    checked = [manifest["source"], *manifest["outputs"], *manifest["evidence"],
        manifest["snapshot"], manifest["swiss_binding"], manifest["candidate_readback"]]
    for ref in checked:
        check(ref)
    names = {f"{STEM}.{suffix}" for suffix in ("png", "pdf", "svg")}
    names.update(("scores.tsv", "coverage.tsv", "swiss_intervals.tsv"))
    require(len(manifest["outputs"]) == len(names)
            and {Path(r["path"]).name for r in manifest["outputs"]} == names
            and all(record(root / Path(r["path"]).name) == r for r in manifest["outputs"]),
            "Changed figure output inventory")
    snapshot = json.loads(Path(manifest["snapshot"]["path"]).read_text())
    binding = json.loads(Path(manifest["swiss_binding"]["path"]).read_text())
    readback = json.loads(Path(manifest["candidate_readback"]["path"]).read_text())
    require(binding["snapshot"] == manifest["snapshot"] and readback["binding"] == manifest["swiss_binding"],
            "Cross-snapshot figure binding")
    native = {r["cell"]: r for r in snapshot["rows"]}
    require(sum(r["accuracy_admitted"] is True for r in native.values()) == 3
            and manifest["unavailable_score_cells"] == [r["cell"] for r in snapshot["rows"] if not r["accuracy_admitted"]]
            and len(manifest["unavailable_score_cells"]) == 4, "Changed missing-score scope")
    require(native[CELLS[1]]["resources"] is None and native[CELLS[1]]["timing_eligible"] is False
            and native[CELLS[1]]["timing_admitted"] is False, "Recovered timing relabeled")
    scores = table(root / "scores.tsv")
    expected = {(c, e): v for c in CELLS for e, v in native[c]["scores"].items()}
    require(len(scores) == len(expected) == 18 and len({(r["cell"], r["endpoint"]) for r in scores}) == 18
            and {(r["cell"], r["endpoint"]): float(r["value"]) for r in scores} == expected
            and all(r["statistic"] == ("F1" if r["endpoint"] in ("VGNC", "SwissTrees", "TreeFam-A")
                else "similarity") for r in scores), "Changed plotted score/statistic")
    coverage = table(root / "coverage.tsv")
    require([r["cell"] for r in coverage] == list(CELLS), "Changed coverage row inventory")
    for row in coverage:
        original = native[row["cell"]]
        for name in ("input_accessions", "relation_accessions", "submitted_pairs"):
            require(int(row[name]) == original[name], "Changed coverage count")
        require(float(row["relation_coverage"]) == original["relation_accessions"]/original["input_accessions"]
                and row["prediction_semantics"] == original["prediction_semantics"]
                and all(float(row[k]) == original["endpoint_details"]["SwissTrees"][k] for k in ("precision", "recall")),
                "Changed coverage/precision/recall value or semantics")
    intervals = table(root / "swiss_intervals.tsv")
    require([(r["contrast"], r["metric"]) for r in intervals] ==
            [(c, m) for c in CONTRASTS for m in ("F1", "PPV", "TPR")], "Changed plotted interval inventory")
    require(binding["multiplicity_endpoints"] == 42 and binding["replicates_reused"] == 100000
            and binding["new_bootstrap_draws"] == 0, "Changed interval adjustment")
    for row in intervals:
        selected = [r for r in binding["contrasts"] if r["name"] == row["contrast"]]
        require(len(selected) == 1 and selected[0]["status"] == "native_records_matched", "Unmatched plotted interval")
        value = selected[0]["metrics"][row["metric"]]
        expected = dict(difference_pp=100*value["difference"], nominal_low_pp=100*value["paired_percentile_ci"][0],
            nominal_high_pp=100*value["paired_percentile_ci"][1], adjusted_low_pp=100*value["bonferroni_percentile_ci"][0],
            adjusted_high_pp=100*value["bonferroni_percentile_ci"][1])
        require({k: float(v) for k, v in row.items() if k not in ("contrast", "metric")} == expected,
                "Changed plotted interval value")
        if row["metric"] == "F1":
            require(expected["adjusted_low_pp"] < 0 < expected["adjusted_high_pp"], "Changed F1 interpretation")
    svg = ET.parse(root / f"{STEM}.svg")
    labels = " ".join(node.text or "" for node in svg.iter() if node.tag.endswith("}text"))
    require(all(label in labels for label in LABELS), "Missing SVG scope label")
    png = pixels(root / f"{STEM}.png")
    require(png["size"] == [2760, 1960], "Changed raster dimensions")
    with fitz.open(root / f"{STEM}.pdf") as document:
        require(len(document) == 1, "Changed PDF page inventory")
        text = " ".join(document[0].get_text().split())
        require(all(label in text for label in LABELS), "Missing PDF scope label")
        document[0].get_pixmap(matrix=fitz.Matrix(120/72, 120/72)).save(pdf_preview)
    for ref in checked:
        check(ref)
    return dict(schema="native_qfo_three_cell_figure_readback_v1", manifest=record(manifest_path),
        snapshot=manifest["snapshot"], swiss_binding=manifest["swiss_binding"], source=record(__file__),
        checked_evidence_records=len(checked), score_endpoints_checked=18, coverage_endpoints_checked=3,
        precision_recall_endpoints_checked=6, contrast_endpoints_checked=6, exact_table_readback=True,
        png=png, pdf_preview=pixels(pdf_preview), pdf_preview_record=record(pdf_preview),
        svg_and_decoded_pdf_required_labels_present=True, automatic_visual_certification=False,
        new_scoring_or_admission=False, scientific_timings_admitted=False, publication_ready=False,
        limitations=["Presentation readback of admitted snapshots, not new raw scientific admission or uncertainty certification.",
            "Decoded PDF/pixel/text checks do not replace separate visual inspection.",
            "Three of seven fresh cells; no missing-score imputation, timing repair or independent confirmation."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("figure-directory", "pdf-preview", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists() and not args.output.is_symlink(), "Output already exists")
    result = review(args.figure_directory, args.pdf_preview)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({k: result[k] for k in ("score_endpoints_checked", "coverage_endpoints_checked", "contrast_endpoints_checked")}))
