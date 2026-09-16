"""Render seed-level simulation results without imputing failed scores."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.summarize_simulation_panel import CONDITIONS, METHODS, PANELS

LABELS = {"orthohmm_high_sensitivity": "OrthoHMM high sensitivity",
          "orthohmm_satellite_v2": "OrthoHMM satellite_v2",
          "orthofinder_full": "OrthoFinder full",
          "orthofinder_sequence_only": "OrthoFinder sequence checkpoint"}


def percent(value):
    return "NA" if value is None else f"{100 * value:.2f}"


def interval(value):
    return "NA" if value is None else f"[{value[0]:.2f}, {value[1]:.2f}]"


def render(report):
    panel = report["panel"]
    if panel not in PANELS:
        raise ValueError("Unknown simulation panel")
    lines = [f"# Simulation Results: {panel}", "",
             "Scores are percentages and means of completed seed-level metrics, not pooled gene pairs.",
             "Failed or inapplicable seeds have no imputed score. Available-case means are conditional on success.", "",
             "| Condition | Method | Completed | Failed | Inapplicable | Precision | Recall | F1 |",
             "|---|---|---:|---:|---:|---:|---:|---:|"]
    for condition in CONDITIONS:
        for method in METHODS:
            row = report["conditions"][condition]["methods"][method]
            means = row["available_case_means"]
            lines.append(f"| {condition} | {LABELS[method]} | {len(row['complete_seeds'])} | {len(row['failed_seeds'])} | "
                         f"{len(row['inapplicable_seeds'])} | {percent(means['precision'])} | {percent(means['recall'])} | {percent(means['f1'])} |")
    lines.extend(["", "## Paired F1 Differences", "",
                  "OrthoHMM minus full OrthoFinder, in percentage points, using only successful paired seeds.", "",
                  "| Condition | Method | Paired Seeds | Difference | Nominal 95% CI | Bonferroni-14 CI |",
                  "|---|---|---:|---:|---|---|"])
    for condition in CONDITIONS:
        for method in METHODS[:2]:
            row = report["conditions"][condition]["contrasts"][method]
            metric = row.get("metrics", {}).get("f1", {})
            difference = metric.get("difference_percentage_points")
            value = "NA" if difference is None else f"{difference:.2f}"
            lines.append(f"| {condition} | {LABELS[method]} | {len(row['included_seeds'])} | {value} | "
                         f"{interval(metric.get('paired_95_percent_ci'))} | {interval(metric.get('bonferroni_14_ci'))} |")
    lines.extend(["", "The sequence checkpoint requires a valid full parent run and has no independent runtime.",
                  "Intervals use 20,000 whole-seed bootstrap replicates; ten seeds limit tail resolution.",
                  "Resource logs are descriptive shared-machine measurements, not controlled scaling results.",
                  "This panel is not pooled with the other simulation panel and does not establish publication readiness."])
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = args.results.read_bytes()
    text = render(json.loads(data))
    text += f"\nSource results: `{args.results.name}`, SHA-256 `{hashlib.sha256(data).hexdigest()}`.\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as handle:
        handle.write(text)


if __name__ == "__main__":
    main()
