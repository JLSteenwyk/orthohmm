"""Exploratory endpoint reaggregation of audited CPU counters, never eligibility."""

import argparse
import gzip
import json
from pathlib import Path

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dual_cpu_brackets import compare

AUDIT_SHA = "f6451d1fd5a96ff5f1a9b2155c6e4d9f17fabeb4de8a86cd53fbd1e70f37bff3"
WIDTHS = (5, 10, 30, 60)


def blocks(interval_count, width):
    if (type(interval_count) is not int or interval_count < 1
            or type(width) is not int or width < 1):
        raise ValueError("Require positive integer interval count and width")
    return [(left, min(left + width, interval_count)) for left in range(0, interval_count, width)]


def aggregate(run, measured):
    if run["screening"] != measured["screening"] or run["job_id"] != measured["job_id"]:
        raise ValueError("Audited and retained measurement differ")
    points, screen = measured["points"], measured["screening"]
    if len(points) < 2 or len(screen["narrow_intervals"]) != len(points) - 1:
        raise ValueError("Incomplete single-interval evidence")
    singles = [compare(a, b, run["job_id"]) for a, b in zip(points, points[1:])]
    if ([r["narrow"] for r in singles] != screen["narrow_intervals"]
            or [r["outer"] for r in singles] != screen["original_screening"]["original_threshold_screen"]["intervals"]):
        raise ValueError("Single-interval counter replay differs")
    flags = [i for i, r in enumerate(singles) if not r["narrow"]["screen_passed"]]
    if flags != run["narrow_flagged_intervals"] or flags != screen["narrow_flagged_intervals"]:
        raise ValueError("Original narrow flags differ")
    scales = []
    for width in WIDTHS:
        rows = []
        for start, end in blocks(len(singles), width):
            # Re-read endpoints arithmetically; do not sum overlapping host brackets.
            value = compare(points[start], points[end], run["job_id"], enforce_gap=False)
            rows.append(dict(start_point=start, end_point=end, observation_steps=end-start,
                remainder=end-start < width, result=value,
                original_narrow_flagged_intervals=[i for i in flags if start <= i < end]))
        scales.append(dict(observation_steps_per_block=width, blocks=rows,
            narrow_flagged_blocks=sum(not r["result"]["narrow"]["screen_passed"] for r in rows),
            outer_flagged_blocks=sum(not r["result"]["outer"]["screen_passed"] for r in rows),
            flagged_original_intervals_in_narrow_passing_blocks=sum(
                len(r["original_narrow_flagged_intervals"]) for r in rows if r["result"]["narrow"]["screen_passed"])))
    whole = compare(points[0], points[-1], run["job_id"], enforce_gap=False)
    if whole != screen["observation_window"]:
        raise ValueError("Whole observation-window replay differs")
    return dict(index=run["index"], method=run["method"], job_id=run["job_id"],
        original_intervals=singles, original_narrow_flagged_intervals=flags,
        scales=scales, whole_observation_window=whole,
        scientific_timings_admitted=False, controlled_workload_verified=False)


def report(audit_path):
    source = record(__file__)
    helpers = [record(Path(__file__).with_name(name)) for name in (
        "probe_dual_cpu_brackets.py", "probe_interval_cpu.py", "probe_cgroup_frontier.py",
        "measure_native_frontier_step.py", "measure_native_hierarchy_step.py",
        "probe_native_pressure.py", "probe_host_counters.py")]
    audit_record = record(audit_path)
    if audit_record["sha256"] != AUDIT_SHA:
        raise ValueError("Require the pinned complete native audit")
    audit = json.loads(gzip.decompress(audit_path.read_bytes()))
    if audit["validated_tasks"] != 3 or [r["index"] for r in audit["runs"]] != [0, 1, 2]:
        raise ValueError("Incomplete native diagnostic panel")
    evidence, results = [audit_record], []
    for run in audit["runs"]:
        if run["status"] != "validated":
            raise ValueError("Cannot omit failed native tasks")
        items = [r for r in run["inventory"] if Path(r["path"]).name == "dual_bracket_report.json"]
        if len(items) != 1:
            raise ValueError("Ambiguous measurement report")
        check(items[0])
        evidence.extend(items)
        results.append(aggregate(run, json.loads(Path(items[0]["path"]).read_text())))
    for item in [source, *evidence, *helpers]:
        check(item)
    return dict(status="dual_cpu_aggregation_described_not_admitted", runs=results, evidence=evidence,
        source=source, helpers=helpers, scientific_timings_admitted=False,
        controlled_workload_verified=False, publication_ready=False,
        limitations=["Post-outcome exploratory scales, not a prospective scientific inclusion policy.",
            "Widths count observation steps, not exact seconds; actual durations remain in each result.",
            "Every original interval and final partial block is retained; coarser passing blocks do not erase flags.",
            "Shared host read brackets overlap even across adjacent endpoint blocks; do not sum their residuals.",
            "Averaging can dilute bursts; passing longer windows does not establish absence of interference.",
            "Observation windows include wrapper work and are not exact native-command boundaries.",
            "No historical timing admission, causal attribution, overhead subtraction or selective rerun."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = report(args.audit.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
