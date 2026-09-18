"""Export all six original-release QfO factorial endpoints and pair counts."""

import argparse
import csv
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_qfo_factorial_swiss import selected_admission
from benchmark_tools.bootstrap_qfo_factorial import CELLS
from benchmark_tools.export_qfo_corrected_comparison import ENDPOINTS
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen


def extract(report, index):
    native = selected_admission(report, index)
    assessment = native["assessment"]
    conversion = report["stage"] if index in (0, 4) else report["conversion"]
    if (conversion["status"] != "cell_pairs_prepared_unscored" or conversion["cell"] != CELLS[index]
            or conversion["index"] != index or set(assessment["endpoints"]) != set(ENDPOINTS)):
        raise ValueError("Wrong conversion cell or incomplete endpoints")
    semantics = "native phylogenetically inferred pairs" if index % 2 else "cross-species group-derived clique pairs"
    if conversion["semantics"] != semantics:
        raise ValueError("Wrong prediction semantics")
    total, retained, removed = (conversion[k] for k in ("total_pairs", "retained_pairs", "removed_mapping_pairs"))
    if (any(type(v) is not int or v < 0 for v in (total, retained, removed))
            or total == 0 or total != retained + removed):
        raise ValueError("Invalid pair accounting")
    scores, details = {}, {}
    for endpoint in ENDPOINTS:
        value = assessment["endpoints"][endpoint]
        participant = value["native_participant"]
        x, y = participant["metric_x"], participant["metric_y"]
        if participant["participant_id"] != assessment["participant"]:
            raise ValueError("Mixed native participants")
        if any(type(v) not in (int, float) or not math.isfinite(v) or v < 0 for v in (x, y)) or y > 1:
            raise ValueError("Invalid native statistic")
        if endpoint in ("VGNC", "SwissTrees", "TreeFam-A"):
            if x > 1:
                raise ValueError("Invalid recall")
            score = 2 * x * y / (x + y) if x + y else 0.
            details[endpoint] = {"precision": y, "recall": x}
        else:
            score = y
            details[endpoint] = {"assessed_relations": x, "statistic": value["score_semantics"]}
        if not math.isclose(value["score"], score, rel_tol=0, abs_tol=1e-12):
            raise ValueError("Native score arithmetic differs")
        scores[endpoint] = score
    mean = sum(scores.values()) / 6
    if not math.isclose(mean, assessment["secondary_six_metric_mean"], rel_tol=0, abs_tol=1e-12):
        raise ValueError("Secondary mean differs")
    return {"cell": CELLS[index], "scores": scores, "details": details, "secondary_mean": mean,
            "total_pairs": total, "retained_pairs": retained, "removed_mapping_pairs": removed,
            "prediction_semantics": semantics, "participant": assessment["participant"]}


def export(inventory, sha, output):
    if output.exists():
        raise FileExistsError(output)
    manifest = read_frozen(inventory, sha)
    if [r["cell"] for r in manifest["cells"]] != list(CELLS):
        raise ValueError("Require eight ordered cell admissions")
    checked, rows = [record(inventory)], []
    for index, entry in enumerate(manifest["cells"]):
        check(entry["admission"])
        report = json.loads(Path(entry["admission"]["path"]).read_text())
        pair_record = report["pairs_manifest"]
        check(pair_record)
        conversion = json.loads(Path(pair_record["path"]).read_text())
        if conversion != (report["stage"] if index in (0, 4) else report["conversion"]):
            raise ValueError("Embedded conversion differs from bound file")
        rows.append({**extract(report, index), "admission": entry["admission"], "conversion": pair_record})
        checked.extend([entry["admission"], pair_record])
    for item in checked:
        check(item)
    output.mkdir(parents=True)
    headers = ["Cell", "GO similarity", "EC similarity", "VGNC F1", "SwissTrees F1", "TreeFam-A F1", "FAS",
               "Secondary mean", "Submitted pairs", "Retained pairs", "Mapping losses", "Prediction semantics"]
    values = [[r["cell"], *[r["scores"][e] for e in ENDPOINTS], r["secondary_mean"], r["total_pairs"],
               r["retained_pairs"], r["removed_mapping_pairs"], r["prediction_semantics"]] for r in rows]
    with (output / "scores.tsv").open("x") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(headers)
        writer.writerows(values)
    lines = ["# Original-Release QfO Factorial", "", "| " + " | ".join(headers) + " |",
             "| " + " | ".join(["---"] * len(headers)) + " |"]
    for row in values:
        lines.append("| " + " | ".join(f"{v:.6f}" if isinstance(v, float) else str(v) for v in row) + " |")
    limitations = ["Original-release 976,504-gene inputs; not corrected-input results or independent validation.",
        "GO/EC similarity and FAS are not F1; the six-metric mean is a project-defined secondary summary.",
        "P-off retains initial HMM search. R changes group-derived pairs to native inferred pairs.",
        "Pair totals measure prediction volume, not the fraction of proteins recovered or benchmark-assessed relations.",
        "No confidence intervals for other QfO endpoints or the secondary mean are supplied by this export.",
        "This validates admitted report identities, conversion binding and score arithmetic, not the full inference workflow anew."]
    (output / "scores.md").write_text("\n".join(lines) + "\n\n" + "\n\n".join(limitations) + "\n")
    result = {"status": "original_qfo_factorial_exported", "cells": rows, "checked_records": checked,
              "source": record(__file__), "helpers": [record(Path(__file__).with_name(n)) for n in
                  ("audit_qfo_factorial_swiss.py", "bootstrap_qfo_factorial.py", "export_qfo_corrected_comparison.py")],
              "publication_ready": False, "limitations": limitations,
              "outputs": [record(output / name) for name in ("scores.tsv", "scores.md")]}
    (output / "manifest.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    export(args.inventory.resolve(), args.sha256, args.output.resolve())
