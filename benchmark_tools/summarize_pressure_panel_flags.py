"""Describe previously audited interval flags without changing timing eligibility."""

import argparse
from collections import Counter
import gzip
import json
import math
from pathlib import Path
import statistics

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


def distribution(values):
    if not values or any(type(v) not in (int, float) or not math.isfinite(v) for v in values):
        raise ValueError("Require nonempty finite observations")
    return dict(minimum=min(values), median=statistics.median(values), maximum=max(values))


def unique_report(evidence):
    candidates = {}
    for item in evidence:
        if Path(item["path"]).name != "frontier_report.json":
            continue
        previous = candidates.setdefault(item["path"], item)
        if previous != item:
            raise ValueError("Conflicting audited report records")
    if len(candidates) != 1:
        raise ValueError("Require unique audited report")
    return next(iter(candidates.values()))


def summarize(screen):
    original = screen["original_threshold_screen"]
    intervals = original["intervals"]
    flags = [i for i, row in enumerate(intervals) if not row["screen_passed"]]
    if flags != original["flagged_intervals"]:
        raise ValueError("Stored flags differ from interval screens")
    frontier = screen["frontier_intervals"]
    if len(frontier) != len(intervals):
        raise ValueError("Interval inventories differ")
    reasons = Counter(reason for row in intervals for reason in row["reasons"])
    return dict(intervals=len(intervals), flagged_intervals=len(flags), reasons=dict(reasons),
                whole_command=original["whole_command_screen"],
                signed_residual_cores=distribution([r["signed_unassigned_average_cores"] for r in intervals]),
                outer_read_overhang_s=distribution([r["outer_read_overhang_s"] for r in intervals]),
                outside_frontier_cpu_s=distribution([r["outside_target_frontier_cpu_s"] for r in frontier]),
                root_minus_frontier_cpu_s=distribution([r["root_minus_frontier_cpu_s"] for r in frontier]),
                native_pressure_whole_command=screen["native_pressure_whole_command"])


def report(audit_path, expected_sha):
    source = record(audit_path)
    if source["sha256"] != expected_sha:
        raise ValueError("Audit checksum differs")
    audit = json.loads(gzip.decompress(audit_path.read_bytes()))
    if [r["index"] for r in audit["runs"]] != list(range(18)):
        raise ValueError("Require all 18 audited tasks")
    rows, evidence = [], [source]
    for run in audit["runs"]:
        row = {k: run[k] for k in ("index", "method", "mode", "status")}
        if run["status"] == "validated" and run["mode"] == "periodic":
            item = unique_report(run["evidence"])
            check(item)
            row["diagnostic"] = summarize(json.loads(Path(item["path"]).read_text())["screening"])
            evidence.append(item)
        rows.append(row)
    for item in evidence:
        check(item)
    return dict(status="audited_pressure_flags_described", runs=rows, evidence=evidence,
                source=record(__file__), scientific_timings_admitted=False,
                limitations=["Descriptive summaries of audited reports, not an independent raw-counter replay.",
                             "Different counter windows prevent causal subtraction or timing correction.",
                             "No thresholds changed; failed tasks and unavailable intervals remain visible.",
                             "Pressure includes native work and does not identify external interference."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", required=True, type=Path)
    parser.add_argument("--audit-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = report(args.audit.resolve(), args.audit_sha256)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
