#!/usr/bin/env python3
"""Summarize the six QfO metrics retained in the OrthoHMM comparison."""

import argparse
import json
from pathlib import Path


METRICS = {
    "VGNC F": ("VGNC/VGNC.json", "f_score"),
    "SwissTrees F": ("SwissTrees/SwissTrees.json", "f_score"),
    "TreeFam-A F": ("TreeFam-A/TreeFam-A.json", "f_score"),
    "EC": ("EC/EC.json", "metric_y"),
    "GO": ("GO/GO.json", "metric_y"),
    "FAS": ("FAS/FAS.json", "metric_y"),
}


def harmonic_mean(x, y):
    return 0.0 if x + y == 0 else 2.0 * x * y / (x + y)


def read_metric(path, mode):
    data = json.loads(path.read_text())
    participants = data["datalink"]["inline_data"]["challenge_participants"]
    if len(participants) != 1:
        raise ValueError(f"Expected one participant in {path}, found {len(participants)}")
    participant = participants[0]
    if mode == "f_score":
        return harmonic_mean(participant["metric_x"], participant["metric_y"])
    return participant["metric_y"]


def summarize(scoring_dir):
    results_dir = scoring_dir / "results"
    scores = {
        name: read_metric(results_dir / relative_path, mode)
        for name, (relative_path, mode) in METRICS.items()
    }
    return {
        "prediction": scoring_dir.name,
        "scores": scores,
        "mean": sum(scores.values()) / len(scores),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scoring_dirs", nargs="+", type=Path)
    parser.add_argument("--json", type=Path, dest="json_path")
    args = parser.parse_args()

    summaries = [summarize(path.resolve()) for path in args.scoring_dirs]
    if args.json_path:
        args.json_path.parent.mkdir(parents=True, exist_ok=True)
        args.json_path.write_text(json.dumps({"records": summaries}, indent=2) + "\n")

    columns = ["prediction", *METRICS, "mean"]
    print("| " + " | ".join(columns) + " |")
    print("| " + " | ".join(["---"] + ["---:"] * (len(columns) - 1)) + " |")
    for summary in summaries:
        values = [summary["scores"][name] for name in METRICS]
        row = [summary["prediction"], *(f"{value:.3f}" for value in values)]
        row.append(f"{summary['mean']:.3f}")
        print("| " + " | ".join(row) + " |")


if __name__ == "__main__":
    main()
