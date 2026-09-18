"""Export admitted corrected-release QfO scores with explicit missing rows."""

import argparse
import csv
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.publication_comparison import METHODS
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.prepare_qfo_corrected_orthofinder_pairs import SEMANTICS as OF_SEMANTICS

ENDPOINTS = ("GO", "EC", "VGNC", "SwissTrees", "TreeFam-A", "FAS")
METHOD_KEYS = {"proteinortho": "proteinortho_6_3_6", "sonic": "sonicparanoid_2_0_9",
               "orthofinder_full": "orthofinder_3_1_5_full",
               "orthofinder_sequence_only": "orthofinder_3_1_5_sequence_only"}


def extract(report, conversion):
    if (report.get("status") != "corrected_comparator_assessment_admitted"
            or report.get("accuracy_admitted") is not True or report.get("publication_ready") is not False):
        raise ValueError("Require admitted corrected comparator assessment")
    method = report["method"]
    if method not in METHOD_KEYS:
        raise ValueError("Method requires its own audited admission adapter")
    participant = "qfo_corrected_" + method
    assessment = report["assessment"]
    if assessment["participant"] != participant or set(assessment["endpoints"]) != set(ENDPOINTS):
        raise ValueError("Wrong participant or incomplete endpoint set")
    status = ("corrected_orthofinder_pairs_prepared_unscored" if method in OF_SEMANTICS
              else "corrected_comparator_pairs_prepared_unscored")
    if (conversion["status"] != status
            or conversion["method"] != method or conversion["participant"] != participant):
        raise ValueError("Wrong corrected conversion binding")
    if method in OF_SEMANTICS and conversion["semantics"] != OF_SEMANTICS[method]:
        raise ValueError("Wrong OrthoFinder prediction semantics")
    total, retained, removed = (conversion[k] for k in ("total_pairs", "retained_pairs", "removed_mapping_pairs"))
    if any(type(v) is not int or v < 0 for v in (total, retained, removed)) or total != retained + removed or total == 0:
        raise ValueError("Invalid submitted-pair accounting")
    scores, details = {}, {}
    for endpoint in ENDPOINTS:
        result = assessment["endpoints"][endpoint]
        native = result["native_participant"]
        if native["participant_id"] != participant:
            raise ValueError("Mixed endpoint participant")
        x, y = native["metric_x"], native["metric_y"]
        if any(type(v) not in (int, float) or not math.isfinite(v) or v < 0 for v in (x, y)) or y > 1:
            raise ValueError("Invalid native endpoint value")
        if endpoint in ("VGNC", "SwissTrees", "TreeFam-A"):
            if x > 1:
                raise ValueError("Recall outside unit interval")
            score = 2 * x * y / (x + y) if x + y else 0.
            details[endpoint] = {"recall": x, "precision": y, "statistic": "F1"}
        else:
            score = y
            details[endpoint] = {"assessed_relations": x, "statistic": result["score_semantics"]}
        if not math.isclose(result["score"], score, rel_tol=0, abs_tol=1e-12):
            raise ValueError("Score differs from native statistic")
        scores[endpoint] = score
    mean = sum(scores.values()) / len(ENDPOINTS)
    if not math.isclose(assessment["secondary_six_metric_mean"], mean, rel_tol=0, abs_tol=1e-12):
        raise ValueError("Secondary mean differs")
    return {"key": METHOD_KEYS[method], "status": "admitted", "scores": scores, "details": details,
            "secondary_mean": mean, "submitted_pairs": total, "retained_pairs": retained,
            "removed_mapping_pairs": removed, "prediction_semantics": conversion["semantics"]}


def export(sources, output):
    if output.exists():
        raise FileExistsError(output)
    rows, checked = {}, []
    for path, sha in sources:
        path = Path(path).resolve()
        report = read_frozen(path, sha)
        source_record = record(path)
        pair_record = report["pairs_manifest"]
        check(pair_record)
        conversion = json.loads(Path(pair_record["path"]).read_text())
        row = extract(report, conversion)
        if row["key"] in rows:
            raise ValueError("Duplicate corrected method admission")
        row["admission"] = source_record
        row["conversion"] = pair_record
        rows[row["key"]] = row
        checked.extend([source_record, pair_record])
    table = []
    for key, label, _, semantics in METHODS:
        row = rows.get(key, {"key": key, "status": "not_admitted", "scores": {e: None for e in ENDPOINTS},
                            "secondary_mean": None, "submitted_pairs": None, "retained_pairs": None,
                            "removed_mapping_pairs": None, "prediction_semantics": semantics})
        table.append({**row, "label": label})
    for item in checked:
        check(item)
    output.mkdir(parents=True, exist_ok=False)
    headers = ["Method", "Status", "GO similarity", "EC similarity", "VGNC F1", "SwissTrees F1",
               "TreeFam-A F1", "FAS", "Secondary mean", "Submitted pairs", "Retained pairs"]
    values = [[r["label"], r["status"], *[r["scores"][e] for e in ENDPOINTS], r["secondary_mean"],
               r["submitted_pairs"], r["retained_pairs"]] for r in table]
    with (output / "scores.tsv").open("x") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(headers)
        writer.writerows(values)
    lines = ["# Corrected-Release QfO Scores", "", "| " + " | ".join(headers) + " |",
             "| " + " | ".join(["---"] * len(headers)) + " |"]
    for row in values:
        rendered = ["pending" if v is None else f"{v:.6f}" if isinstance(v, float) else str(v) for v in row]
        lines.append("| " + " | ".join(rendered) + " |")
    limitations = ["Only admitted corrected-release reports are included; pending is not zero and does not describe scheduler state.",
        "GO/EC similarity and FAS are not F1; the six-metric mean is a project-defined secondary summary.",
        "No paired uncertainty or corrected-release ranking is established by this point-estimate table.",
        "This exporter checks admitted report hashes and native score arithmetic, not the complete inference/scorer workflow again."]
    (output / "scores.md").write_text("\n".join(lines) + "\n\n" + "\n\n".join(limitations) + "\n")
    result = {"status": "corrected_qfo_partial_comparison", "methods": table, "checked_records": checked,
              "source": record(__file__), "method_registry_source": record(Path(__file__).with_name("publication_comparison.py")),
              "publication_ready": False, "limitations": limitations,
              "outputs": [record(output / name) for name in ("scores.tsv", "scores.md")]}
    (output / "manifest.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assessment", nargs=2, action="append", required=True, metavar=("PATH", "SHA256"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    export(args.assessment, args.output.resolve())
