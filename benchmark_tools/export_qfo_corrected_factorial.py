"""Export corrected factorial admissions without substituting historical scores."""

import argparse
import csv
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.export_qfo_corrected_comparison import ENDPOINTS
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_qfo_corrected_factorial_assessment import CELLS, CONVERTERS, validate_stage
from benchmark_tools.run_simulation_methods import read_frozen


def extract(report, conversion):
    if (report.get("status") != "corrected_factorial_assessment_admitted"
            or report.get("accuracy_admitted") is not True or report.get("publication_ready") is not False):
        raise ValueError("Require admitted corrected factorial scores")
    index = report["index"]
    if type(index) is not int or not 0 <= index < 8 or report["cell"] != CELLS[index]:
        raise ValueError("Wrong corrected cell identity")
    if report["conversion"] != conversion:
        raise ValueError("Embedded conversion differs from bound manifest")
    validate_stage(conversion, index, report["conversion_scheduler"])
    assessment = report["assessment"]
    participant = conversion["participant"]
    if assessment["participant"] != participant or set(assessment["endpoints"]) != set(ENDPOINTS):
        raise ValueError("Wrong participant or incomplete endpoints")
    scores, details = {}, {}
    for endpoint in ENDPOINTS:
        result = assessment["endpoints"][endpoint]
        native = result["native_participant"]
        if native["participant_id"] != participant:
            raise ValueError("Mixed endpoint participants")
        x, y = native["metric_x"], native["metric_y"]
        if any(type(v) not in (int, float) or not math.isfinite(v) or v < 0 for v in (x, y)) or y > 1:
            raise ValueError("Invalid native endpoint")
        if endpoint in ("VGNC", "SwissTrees", "TreeFam-A"):
            if x > 1:
                raise ValueError("Invalid recall")
            value = 2 * x * y / (x + y) if x + y else 0.
            details[endpoint] = {"precision": y, "recall": x, "statistic": "F1"}
        else:
            value = y
            details[endpoint] = {"assessed_relations": x, "statistic": result["score_semantics"]}
        if type(result["score"]) not in (int, float) or not math.isclose(result["score"], value, rel_tol=0, abs_tol=1e-12):
            raise ValueError("Native score arithmetic differs")
        scores[endpoint] = value
    mean = sum(scores.values()) / len(ENDPOINTS)
    supplied = assessment["secondary_six_metric_mean"]
    if type(supplied) not in (int, float) or not math.isclose(supplied, mean, rel_tol=0, abs_tol=1e-12):
        raise ValueError("Secondary mean differs")
    return {"index": index, "cell": CELLS[index], "status": "admitted", "scores": scores,
            "details": details, "secondary_mean": mean, "participant": participant,
            "prediction_semantics": conversion["semantics"],
            **{key: conversion[key] for key in ("total_pairs", "retained_pairs", "removed_mapping_pairs")}}


def export(sources, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    rows, checked = {}, []
    for path, digest in sources:
        path = Path(path).resolve()
        report = read_frozen(path, digest)
        admission = record(path)
        pair_record = report["pairs_manifest"]
        check(pair_record)
        conversion = json.loads(Path(pair_record["path"]).read_text())
        row = extract(report, conversion)
        if row["index"] in rows:
            raise ValueError("Duplicate cell admission")
        rows[row["index"]] = {**row, "admission": admission, "conversion": pair_record}
        checked.extend([admission, pair_record])
    table = [rows.get(i, {"index": i, "cell": cell, "status": "not_admitted",
             "scores": {e: None for e in ENDPOINTS}, "secondary_mean": None,
             "total_pairs": None, "retained_pairs": None, "removed_mapping_pairs": None,
             "prediction_semantics": CONVERTERS[i % 2][2]}) for i, cell in enumerate(CELLS)]
    for item in checked:
        check(item)
    output.mkdir(parents=True, exist_ok=False)
    headers = ["Cell", "Status", "GO similarity", "EC similarity", "VGNC F1", "SwissTrees F1",
               "TreeFam-A F1", "FAS", "Secondary mean", "Submitted pairs", "Retained pairs", "Mapping losses", "Prediction semantics"]
    values = [[r["cell"], r["status"], *[r["scores"][e] for e in ENDPOINTS], r["secondary_mean"],
               r["total_pairs"], r["retained_pairs"], r["removed_mapping_pairs"], r["prediction_semantics"]] for r in table]
    with (output / "scores.tsv").open("x") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(headers)
        writer.writerows(values)
    lines = ["# Corrected-Release QfO Factorial", "", "| " + " | ".join(headers) + " |",
             "| " + " | ".join(["---"] * len(headers)) + " |"]
    for row in values:
        lines.append("| " + " | ".join("not admitted" if v is None else f"{v:.6f}" if isinstance(v, float) else str(v) for v in row) + " |")
    limitations = ["Only supplied corrected-release admissions are included; missing is not zero or a scheduler-state claim.",
        "GO/EC similarity and FAS are not F1; the six-metric mean is a project-defined secondary summary.",
        "P-off retains initial HMM search. R changes group-derived pairs to native phylogenetic pairs.",
        "Pair totals measure prediction volume, not protein coverage. No paired uncertainty or ranking is established here.",
        "This export checks admission hashes, conversion binding and score arithmetic, not the complete inference/scorer workflow anew."]
    (output / "scores.md").write_text("\n".join(lines) + "\n\n" + "\n\n".join(limitations) + "\n")
    result = {"status": "corrected_qfo_factorial_exported", "admitted_cells": len(rows), "cells": table,
              "checked_records": checked, "source": record(__file__), "publication_ready": False,
              "limitations": limitations, "outputs": [record(output / n) for n in ("scores.tsv", "scores.md")],
              "helpers": [record(Path(__file__).with_name(n)) for n in (
                  "run_qfo_corrected_factorial_assessment.py", "bootstrap_qfo_factorial.py", "export_qfo_corrected_comparison.py")]}
    (output / "manifest.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assessment", nargs=2, action="append", default=[], metavar=("PATH", "SHA256"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    export(args.assessment, args.output.absolute())
