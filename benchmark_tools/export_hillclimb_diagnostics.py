"""Export endpoint diagnostics and historical attempts without mixing metric definitions."""

import argparse
import csv
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

ENDPOINTS = {"OrthoBench": ("OrthoBench F1",),
    "QfO": ("VGNC F", "SwissTrees F", "TreeFam-A F", "EC", "GO", "FAS"),
    "ThreeKingdoms": ("Three Kingdoms F1",)}
ADMITTED = {"retained_registered_result", "independent_content_readback_verified"}
COLUMNS = ("method", "label", "dataset", "endpoint", "statistic", "score", "precision", "recall",
    "cell_index", "cell_status", "native_genes", "native_species", "native_groups", "native_pairs",
    "submitted_pairs", "retained_pairs", "assessed_relations", "benchmark_genes_in_orthogroups",
    "benchmark_reported_gene_coverage", "reference_genes", "reference_genes_in_prediction",
    "reference_gene_coverage", "true_positive_gene_pairs", "false_positive_gene_pairs", "false_negative_gene_pairs",
    "benchmark_predicted_orthogroups", "exact_refogs", "semantics", "detail_source_document", "detail_source_location",
    "score_source_document", "score_source_location", "coverage_source_location")


def unit(value):
    if value is not None and (type(value) not in (int, float) or not math.isfinite(value) or not 0 <= value <= 1):
        raise ValueError("Expected a finite 0-to-1 score, precision, recall or coverage")
    return value


def count(value):
    if value is not None and (type(value) not in (int, float) or not math.isfinite(value) or value < 0 or int(value) != value):
        raise ValueError("Expected a nonnegative integral count")
    return value


def cell_rows(method, method_index, cell, index, source, retained, retained_source):
    dataset = cell["dataset"]
    if dataset not in ENDPOINTS:
        raise ValueError("Unknown diagnostic dataset")
    admitted = cell["status"] in ADMITTED
    legacy = cell["status"] == "retained_registered_result"
    scores = method["scores"] if legacy else cell.get("scores", {})
    details = cell.get("retained_details", {}) if legacy else cell.get("details", {})
    coverage = cell.get("coverage", {})
    location = ["methods", method_index, "cells", index]
    result = []
    for endpoint in ENDPOINTS[dataset]:
        name = dataset + ": " + endpoint
        score = unit(scores.get(name)) if admitted else None
        if admitted and (score is None or score != method["scores"].get(name)):
            raise ValueError("Admitted diagnostic score does not reproduce the current comparison")
        row = dict(method=method["key"], label=method["label"], dataset=dataset, endpoint=endpoint,
            statistic="conditional_similarity" if endpoint in ("EC", "GO", "FAS") else "F1",
            score=score, cell_index=index, cell_status=cell["status"], semantics=cell.get("semantics"),
            score_source_document=str(source),
            score_source_location=["methods", method_index, "scores", name] if legacy else [*location, "scores", name],
            detail_source_document=str(source), detail_source_location=[*location, "retained_details" if legacy else "details"],
            coverage_source_location=[*location, "coverage"] if coverage else None)
        for key in ("genes", "species", "groups", "pairs"):
            row["native_" + key] = count(coverage.get(key))
        if not admitted:
            result.append(row)
            continue
        if dataset == "OrthoBench":
            if legacy:
                matches = [(i, r) for i, r in enumerate(retained["rows"]) if r["key"] == method["key"]]
                if len(matches) != 1:
                    raise ValueError("Require one retained OrthoBench evidence row")
                retained_index, retained_row = matches[0]
                evidence = retained_row["orthobench_retained_evidence"]
                prediction = evidence.get("prediction_provenance") or retained_row.get("orthobench_supplemental_readback", {}).get("prediction")
                if (prediction is None or prediction != details["retained_ob_metadata"]["prediction"]
                        or abs(evidence["f_score_percent"] / 100 - score) > 1e-12):
                    raise ValueError("Retained OrthoBench diagnostics are not bound to the scored prediction")
                row.update(precision=unit(evidence["precision_percent"] / 100),
                    recall=unit(evidence["recall_percent"] / 100),
                    benchmark_genes_in_orthogroups=count(evidence.get("genes_in_orthogroups")),
                    benchmark_reported_gene_coverage=unit(evidence["coverage_percent"] / 100),
                    exact_refogs=count(evidence.get("exact_refogs")), detail_source_document=str(retained_source),
                    detail_source_location=["rows", retained_index, "orthobench_retained_evidence"])
                if method["key"] == "orthofinder_3_1_5_full":
                    row["semantics"] = "Official post-phylogenetic orthogroups, not the later debug N0.ids partition"
            else:
                row.update(precision=unit(details.get("precision")), recall=unit(details.get("recall")),
                    exact_refogs=count(details.get("family_details", {}).get("exact_refogs")))
        elif dataset == "QfO":
            key = endpoint.removesuffix(" F") if legacy else endpoint
            entry = details.get("score_details", {}).get(key, {}) if legacy else details.get(key, {})
            row["detail_source_location"] += ["score_details", key] if legacy else [key]
            if legacy:
                row.update(submitted_pairs=count(details.get("submitted_pairs")), retained_pairs=count(details.get("retained_pairs")))
            row.update(precision=unit(entry.get("precision")), recall=unit(entry.get("recall")),
                assessed_relations=count(entry.get("assessed_relations")))
            if row["statistic"] == "conditional_similarity" and (row["precision"] is not None or row["recall"] is not None):
                raise ValueError("Functional similarity must not acquire orthology precision or recall")
            if "similarity" in entry and entry["similarity"] != score:
                raise ValueError("Functional diagnostic does not reproduce the scored endpoint")
        else:
            entry = details["counts"] if legacy else details
            if entry.get("f_score") != score:
                raise ValueError("Three Kingdoms counts do not reproduce the scored endpoint")
            if legacy:
                row["detail_source_location"].append("counts")
            row.update(precision=unit(entry.get("precision")), recall=unit(entry.get("recall")),
                reference_gene_coverage=unit(entry.get("reference_gene_coverage")),
                benchmark_predicted_orthogroups=count(entry.get("predicted_orthogroups")))
            for key in ("reference_genes", "reference_genes_in_prediction", "true_positive_gene_pairs",
                        "false_positive_gene_pairs", "false_negative_gene_pairs"):
                row[key] = count(entry.get(key))
        result.append(row)
    return result


def build(summary, source, retained, retained_source):
    history, current = [], []
    for method_index, method in enumerate(summary["methods"]):
        observations = []
        for index, cell in enumerate(method["cells"]):
            observations.extend(cell_rows(method, method_index, cell, index, source, retained, retained_source))
        history.extend(observations)
        for dataset, endpoints in ENDPOINTS.items():
            for endpoint in endpoints:
                matches = [row for row in observations if row["dataset"] == dataset and row["endpoint"] == endpoint
                           and row["cell_status"] in ADMITTED]
                if matches:
                    # Latest admitted observation, never the maximum score or an average over repeated controls.
                    current.append(dict(matches[-1]))
                else:
                    if dataset + ": " + endpoint in method["scores"]:
                        raise ValueError("Scored endpoint lacks admitted diagnostic ownership")
                    current.append(dict(method=method["key"], label=method["label"], dataset=dataset,
                        endpoint=endpoint, statistic="conditional_similarity" if endpoint in ("EC", "GO", "FAS") else "F1",
                        cell_status="no_admitted_endpoint", score=None, score_source_document=str(source),
                        score_source_location=["methods", method_index, "scores"]))
    return current, history


def write_table(path, rows):
    with path.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "NA" if row.get(key) is None else json.dumps(row[key]) if isinstance(row[key], list)
                else row[key] for key in COLUMNS})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", type=Path, required=True)
    parser.add_argument("--retained-scores", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("Preserve prior diagnostic exports")
    inputs = [record(args.comparison), record(args.retained_scores), record(__file__),
        record(Path(__file__).with_name("prepare_ob_candidate_neighborhood.py"))]
    summary = json.loads(args.comparison.read_text())
    retained = json.loads(args.retained_scores.read_text())
    current, history = build(summary, inputs[0]["path"], retained, inputs[1]["path"])
    for item in inputs:
        check(item)
    args.output.mkdir(parents=True)
    tables = []
    for name, rows in (("endpoint_diagnostics.tsv", current), ("historical_cell_diagnostics.tsv", history)):
        path = args.output / name
        write_table(path, rows)
        tables.append(record(path))
    report = dict(status="hillclimb_endpoint_and_history_diagnostics_exported", inputs=inputs, tables=tables,
        current=current, history=history, defaults_promoted=False, goal_completed=False,
        limitations=["Scores and precision/recall are 0-to-1; OrthoBench percent diagnostics are explicitly divided by 100.",
            "Native universe counts, reported benchmark assignment coverage, reference-gene coverage, submitted/retained pairs and assessed functional relations are different quantities.",
            "Missing fields remain NA. Historical attempts stay separate; the current table selects the latest admitted cell, not the highest score or an average.",
            "OrthoBench is weighted group F1; QfO reference endpoints are pair F1; Three Kingdoms counts reference co-membership including within-species pairs.",
            "GO/EC/FAS conditional similarities have no orthology precision/recall. Assessed populations differ and FAS sampling is unseeded.",
            "Directly retained historical diagnostics are not new raw scoring or proof of historical input consumption. All panels are development-exposed."])
    with (args.output / "report.json").open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"current_endpoints": len(current), "historical_endpoint_observations": len(history),
        "scored_current_endpoints": sum(row.get("score") is not None for row in current)}))


if __name__ == "__main__":
    main()
