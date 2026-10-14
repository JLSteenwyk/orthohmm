"""Post hoc paired-family uncertainty for independently admitted hill-climb cells."""

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.bootstrap_orthobench import METRICS, paired_bootstrap, statistics, weighted_records
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_full_search_hillclimb import runtime
from benchmark_tools.run_membership_hillclimb import write_json

METHODS = ("shared_hit:shared_hit1", "singleton_tie:defer_ambiguous")
PROTOCOL = "hillclimb_recent_ob_uncertainty_20261014.json"


def selected(summary, key):
    methods = [method for method in summary["methods"] if method["key"] == key]
    if len(methods) != 1:
        raise ValueError("Require a unique admitted method")
    cells = [cell for cell in methods[0]["cells"] if cell["dataset"] == "OrthoBench"
             and cell["status"] == "independent_content_readback_verified"]
    if len(cells) != 1 or methods[0]["scores"].get("OrthoBench: OrthoBench F1") != cells[0]["scores"]["OrthoBench: OrthoBench F1"]:
        raise ValueError("Require unique independently verified OrthoBench scores")
    return cells[0]


def sufficient_statistics(cell):
    if (cell["status"] != "independent_content_readback_verified"
            or (cell["coverage"]["genes"], cell["coverage"]["species"]) != (251378, 12)):
        raise ValueError("Require a verified complete OrthoBench universe")
    score = cell["details"]["family_details"]
    names, _, counts = weighted_records(score["refog_records"])
    if len(names) != 70:
        raise ValueError("Require exactly 70 RefOGs")
    actual = statistics(counts.sum(axis=0))
    expected = np.array([score[key] for key in METRICS])
    if (not np.isfinite(expected).all() or not np.allclose(actual, expected, atol=1e-8, rtol=0)
            or abs(actual[0] / 100 - cell["scores"]["OrthoBench: OrthoBench F1"]) > 1e-12):
        raise ValueError("Family sufficient statistics do not reproduce admitted scores")
    return score["refog_records"], expected


def collect(summary):
    cells = {key: selected(summary, key) for key in METHODS}
    controls = {METHODS[0]: cells[METHODS[0]]["matched_control"]["cell"],
                METHODS[1]: cells[METHODS[1]]["matched_control"]}
    baseline_records, baseline_score = sufficient_statistics(controls[METHODS[0]])
    other_records, _ = sufficient_statistics(controls[METHODS[1]])
    if {row["refog"]: row for row in baseline_records} != {row["refog"]: row for row in other_records}:
        raise ValueError("Matched controls have different RefOG records despite aggregate agreement")
    records, points = {"matched_control": baseline_records}, {"matched_control": baseline_score}
    pins = []
    for cell in (*cells.values(), *controls.values()):
        pins.extend(cell["checked"])
    for key, cell in cells.items():
        records[key], points[key] = sufficient_statistics(cell)
    return records, points, pins


def run(readback, output):
    if output.exists():
        raise FileExistsError("Preserve every uncertainty readout")
    protocol_path = Path(__file__).with_name(PROTOCOL)
    protocol = json.loads(protocol_path.read_text())
    if (tuple(protocol["candidates"]) != METHODS or protocol["replicates"] != 100000
            or protocol["seed"] != 20261014 or protocol["alpha"] != .05
            or protocol["multiplicity_endpoints"] != 6):
        raise ValueError("Changed declared uncertainty readout")
    summary = json.loads(readback.read_text())
    records, points, pins = collect(summary)
    scientific_runtime = runtime()
    pins.extend([record(readback), summary["score_table"], record(protocol_path), record(__file__),
                 record(Path(__file__).with_name("bootstrap_orthobench.py")),
                 record(Path(__file__).with_name("score_orthobench_partition.py")),
                 *summary["reporter_files"], *scientific_runtime["files"]])
    pins = list({pin["path"]: pin for pin in pins}.values())
    for pin in pins:
        check(pin)
    arguments = {key: protocol[key] for key in ("replicates", "seed", "alpha", "multiplicity_endpoints")}
    result = paired_bootstrap(records, "matched_control", **arguments)
    repeated = paired_bootstrap(records, "matched_control", **arguments)
    if result != repeated:
        raise ValueError("Paired-family bootstrap does not reproduce deterministically")
    for method, point in result["point_estimates_percent"].items():
        if not np.allclose([point[key] for key in METRICS], points[method], atol=1e-8, rtol=0):
            raise ValueError("Bootstrap point estimate differs from retained scores")
    for pin in pins:
        check(pin)
    if runtime() != scientific_runtime:
        raise ValueError("Scientific runtime changed during readout")
    result.update(status="recent_orthobench_post_hoc_paired_uncertainty_verified", protocol=record(protocol_path),
                  readback=record(readback), runtime=scientific_runtime, checked=pins, deterministic_replay=True,
                  inference_rerun=False, independent_confirmation=False, defaults_promoted=False, goal_completed=False)
    result["limitations"].extend([
        "This readout was declared after observing point estimates; its six-endpoint adjustment does not cover historical tuning or method selection.",
        "Intervals containing zero do not demonstrate equivalence. The 0.005 decline heuristic remains a development screen.",
        "This analysis does not quantify run-to-run, species-sampling, or unseen-dataset uncertainty; QfO and Three Kingdoms are not bootstrapped here."])
    output.mkdir(parents=True)
    with (output / "intervals.tsv").open("x", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["method", "metric", "difference_pp", "paired_lower_pp", "paired_upper_pp", "adjusted_lower_pp", "adjusted_upper_pp"])
        for method, comparison in result["comparisons"].items():
            for metric, values in comparison["metrics"].items():
                writer.writerow([method, metric, values["difference_percentage_points"],
                                 *values["paired_percentile_ci"], *values["bonferroni_percentile_ci"]])
    result["interval_table"] = record(output / "intervals.tsv")
    write_json(output / "report.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--readback", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.readback.resolve(), args.output.resolve())
    print(result["status"])
