"""Report every prespecified recovered QfO contrast without selecting winners."""

import argparse
import hashlib
import json
import math
from pathlib import Path

ADMISSION_SHA = "89b683d0fc7fe9964ce5b6182832bb8eb28758fad3f6bbb4d706926143805956"
STAGES = ("multipass", "multipass_refined", "strict_profiles", "strict_profiles_refined")
METRICS = ("GO", "EC", "VGNC", "SwissTrees", "TreeFam-A", "FAS")
CONTRASTS = ((2, 0), (3, 1), (1, 0), (3, 2))


def summarize(admission):
    if admission["status"] != "four_stage_assessments_checked":
        raise ValueError("Require admitted assessments")
    records = admission["records"]
    if [r["stage"] for r in records] != list(STAGES):
        raise ValueError("Changed stage inventory/order")
    stages = []
    for i, row in enumerate(records):
        if row["status"] != "admitted" or row["index"] != i:
            raise ValueError("Unadmitted stage")
        endpoints = row["assessment"]["endpoints"]
        if set(endpoints) != set(METRICS):
            raise ValueError("Changed metric inventory")
        scores = {m: endpoints[m]["score"] for m in METRICS}
        if any(type(v) not in (int, float) or not math.isfinite(v) or not 0 <= v <= 1
               for v in scores.values()):
            raise ValueError("Invalid score")
        mean = sum(scores.values()) / len(scores)
        if not math.isclose(mean, row["assessment"]["secondary_six_metric_mean"], abs_tol=1e-12):
            raise ValueError("Secondary mean differs")
        conversion = row["conversion"]
        counts = {k: conversion[k] for k in ("total_pairs", "retained_pairs", "removed_mapping_pairs")}
        if any(type(v) is not int or v < 0 for v in counts.values()):
            raise ValueError("Invalid conversion counts")
        if counts["total_pairs"] != counts["retained_pairs"] + counts["removed_mapping_pairs"]:
            raise ValueError("Conversion counts do not reconcile")
        stages.append({"stage": row["stage"], "scores": scores, "native_endpoints": endpoints,
                       "secondary_six_metric_mean": mean, **counts})
    contrasts = []
    for a, b in CONTRASTS:
        contrasts.append({"candidate": STAGES[a], "reference": STAGES[b],
                          "differences": {m: stages[a]["scores"][m] - stages[b]["scores"][m]
                                          for m in METRICS},
                          "secondary_mean_difference": stages[a]["secondary_six_metric_mean"] -
                          stages[b]["secondary_six_metric_mean"],
                          "paired_uncertainty": "not_yet_established"})
    return {"status": "descriptive_prespecified_stage_report", "publication_ready": False,
            "stages": stages, "contrasts": contrasts,
            "limitations": ["Development-exposed cluster-derived predictions, not native reconciled ortholog pairs.",
                            "Profile contrasts include downstream graph and singleton-assignment responses.",
                            "Refined means sequence-based post-clustering refinement, not phylogeny.",
                            "Six-metric mean is a project-defined secondary summary, not official F1.",
                            "Native standard errors do not supply paired difference confidence intervals.",
                            "No significance, independent generalization, superiority or controlled timing claim."]}


def markdown(report):
    lines = ["# Recovered QfO Stage Results", "", "Descriptive results; paired uncertainty remains outstanding.", "",
             "| Stage | " + " | ".join(METRICS) + " | Secondary mean | Retained pairs |",
             "| --- | " + " | ".join(["---:"] * 8) + " |"]
    for row in report["stages"]:
        values = [row["scores"][m] for m in METRICS] + [row["secondary_six_metric_mean"]]
        lines.append("| " + row["stage"] + " | " + " | ".join(f"{v:.6f}" for v in values)
                     + f" | {row['retained_pairs']} |")
    lines += ["", "## Prespecified Differences", "", "Candidate minus reference; raw metric units, not percentage changes.", "",
              "| Contrast | " + " | ".join(METRICS) + " | Secondary mean |",
              "| --- | " + " | ".join(["---:"] * 7) + " |"]
    for row in report["contrasts"]:
        values = [row["differences"][m] for m in METRICS] + [row["secondary_mean_difference"]]
        lines.append(f"| {row['candidate']} - {row['reference']} | " +
                     " | ".join(f"{v:+.6f}" for v in values) + " |")
    lines += ["", "## Native Coordinates", "",
              "Native challenge-specific axes and counts are retained, not pooled across challenges.", "",
              "| Stage | Challenge | X axis | X | Y axis | Y |",
              "| --- | --- | --- | ---: | --- | ---: |"]
    for row in report["stages"]:
        for metric in METRICS:
            endpoint = row["native_endpoints"][metric]
            axes, native = endpoint["axes"], endpoint["native_participant"]
            lines.append(f"| {row['stage']} | {metric} | {axes['x_axis']} | {native['metric_x']} | "
                         f"{axes['y_axis']} | {native['metric_y']} |")
    lines += ["", "## Limitations", "", *["- " + item for item in report["limitations"]], ""]
    return "\n".join(lines)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--admission", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.markdown.exists():
        raise FileExistsError("Require fresh reports")
    data = args.admission.read_bytes()
    if hashlib.sha256(data).hexdigest() != ADMISSION_SHA:
        raise ValueError("Changed frozen admission")
    report = summarize(json.loads(data))
    report["admission_sha256"] = ADMISSION_SHA
    report["source_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open("x") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    with args.markdown.open("x") as handle:
        handle.write(markdown(report))
