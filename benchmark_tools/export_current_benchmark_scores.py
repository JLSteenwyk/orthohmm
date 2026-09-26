"""Consolidate retained scores without mixing original and corrected QfO releases."""

import argparse
import csv
import json
import math
from pathlib import Path

from benchmark_tools.audit_failed_recovery_refinement import record

SOURCES = {
    "publication_comparison_orthomcl_complete_20260916.json": "094f842ad8211450519d4238c0b3a5020a7969049d4b7306f5465d5acbf914a3",
    "qfo_corrected_comparison_20260926_v7/manifest.json": "042aa221d01554114a1b0b413ca3e2ff56dad5f8c06362c69a1bf49f2de642fc",
    "three_kingdoms_comparison_matched_20260918/comparison.json": "ba1ed664ad8ee89ba65b72e7d4d4341eed864957b0720ab1ff4fbe91a398945a",
}
METRICS = ("GO", "EC", "VGNC", "SwissTrees", "TreeFam-A", "FAS")


def index(rows):
    result = {row["key"]: row for row in rows}
    if len(result) != len(rows):
        raise ValueError("Duplicate method key")
    return result


def assemble(ob, qfo, kingdoms):
    if (qfo["status"] != "corrected_qfo_publication_comparison" or qfo["admitted_methods"] != 8
            or kingdoms["status"] != "three_kingdoms_comparison_recounted"
            or any(report["publication_ready"] is not False for report in (ob, qfo, kingdoms))):
        raise ValueError("Unexpected source scope")
    old, corrected = index(ob["methods"]), index(qfo["methods"])
    selected = index([row for row in kingdoms["rows"] if row["use"] == "comparison"])
    if len(corrected) != 8 or set(old) != set(corrected) or set(selected) != set(corrected):
        raise ValueError("Method inventory differs across datasets")
    rows = []
    for key, current in corrected.items():
        if current["status"] != "admitted":
            raise ValueError("Unadmitted QfO method")
        scores = {metric: current["scores"][metric] for metric in METRICS}
        scores.update(OrthoBench=old[key]["orthobench"]["f_score_percent"] / 100,
                      ThreeKingdoms=selected[key]["counts"]["f_score"],
                      QfO_secondary_mean=current["secondary_mean"])
        if any(not isinstance(value, (int, float)) or not math.isfinite(value) or not 0 <= value <= 1
               for value in scores.values()):
            raise ValueError("Invalid score scale")
        if not math.isclose(scores["QfO_secondary_mean"], sum(scores[m] for m in METRICS) / 6, abs_tol=1e-12):
            raise ValueError("Secondary mean disagrees with endpoints")
        rows.append(dict(key=key, label=current["label"], scores=scores,
                         qfo_prediction_semantics=current["prediction_semantics"],
                         orthobench_retained_evidence=old[key]["orthobench"],
                         three_kingdoms={k: selected[key][k] for k in ("run", "input_status", "semantics", "groups")}))
    return rows


def export(results, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    inputs, reports = [], []
    for name, expected in SOURCES.items():
        path = results / name
        item = record(path.resolve())
        if item["sha256"] != expected:
            raise ValueError("Changed retained score source: " + name)
        inputs.append(item)
        reports.append(json.loads(path.read_text()))
    rows = assemble(*reports)
    columns = ("OrthoBench", *METRICS, "QfO_secondary_mean", "ThreeKingdoms")
    limitations = [
        "All displayed values use 0-to-1 units; OrthoBench was converted from percent.",
        "OrthoBench is weighted orthogroup F1; Three Kingdoms is BUSCO-reference pair F1.",
        "QfO GO/EC/FAS are similarities, not F1. Its six-endpoint mean is project-defined and secondary.",
        "No cross-dataset mean or universal ranking is defined; prediction semantics differ.",
        "Three Kingdoms input-consumption gaps remain; only SonicParanoid uses the contemporary matched-input row.",
        "FastOMA uses supplied-tree configurations; the OrthoFinder sequence checkpoint is diagnostic.",
        "Five OrthoBench rows lack direct prediction-file hashes in this source report; full transitive provenance is not consolidated here.",
        "Development-exposed evidence; no inference, raw scoring, uncertainty or resource comparison rerun."]
    output.mkdir(parents=True)
    table = output / "scores.tsv"
    with table.open("x", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["Method", *columns])
        writer.writerows([row["label"], *[row["scores"][c] for c in columns]] for row in rows)
    markdown = output / "scores.md"
    with markdown.open("x") as stream:
        stream.write("# Current Retained Benchmark Scores\n\n")
        stream.write("| Method | " + " | ".join(columns) + " |\n")
        stream.write("| --- | " + " | ".join(["---:"] * len(columns)) + " |\n")
        for row in rows:
            stream.write("| " + row["label"] + " | " + " | ".join(f'{row["scores"][c]:.6f}' for c in columns) + " |\n")
        stream.write("\n" + "\n".join("- " + note for note in limitations) + "\n")
    if [record(item["path"]) for item in inputs] != inputs:
        raise ValueError("Inputs changed during export")
    report = dict(status="current_retained_scores_consolidated", publication_ready=False,
                  inputs=inputs, rows=rows, limitations=limitations, source=record(Path(__file__).resolve()),
                  outputs=[record(table.resolve()), record(markdown.resolve())])
    with (output / "manifest.json").open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    export(args.results.resolve(), args.output.absolute())
