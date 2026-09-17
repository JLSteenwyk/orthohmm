"""Explain inconclusive host intervals without upgrading timing admission."""

import argparse
from collections import Counter
import io
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.command_host_monitor import HostMonitor
from benchmark_tools.snapshot_orthohmm_input_order import record


def unmatched_inventory(rows):
    previous, events, statuses, errors, uncertain = {}, Counter(), Counter(), Counter(), Counter()
    for row in rows:
        if "observation_error" in row:
            raise ValueError("Observation errors require separate review")
        sample = row["snapshot"]
        current = {(p["pid"], p["created"]): p for p in sample["processes"]}
        if len(current) != len(sample["processes"]):
            raise ValueError("Duplicate process identity")
        errors.update(r["type"] for r in sample["errors"])
        interval = row["interval"]
        if interval is not None:
            statuses[interval["status"]] += 1
            uncertain.update(r["reason"] for r in interval["uncertain_processes"])
            for item in interval["unmatched_foreign_processes"]:
                key = (item["pid"], item["created"])
                if key in current and key not in previous:
                    process, direction = current[key], "appeared"
                elif key in previous and key not in current:
                    process, direction = previous[key], "disappeared"
                else:
                    raise ValueError("Unmatched identity not exclusive to one snapshot")
                events[(process["name"], process["cgroup"], direction)] += 1
        previous = current
    return {"interval_counts": dict(statuses), "snapshot_error_types": dict(errors),
            "uncertain_reasons": dict(uncertain),
            "unmatched_identity_events": sum(events.values()),
            "kworker_named_identity_events": sum(n for (name, _, _), n in events.items() if name.startswith("kworker/")),
            "unmatched_events": [{"name": name, "cgroup": group, "direction": direction, "count": count}
                                 for (name, group, direction), count in sorted(events.items())]}


def review(directory):
    measurement_path, raw_path = directory / "results.json", directory / "host_samples.jsonl"
    identities = [record(measurement_path), record(raw_path)]
    measured = json.loads(measurement_path.read_text())
    if identities[1]["sha256"] != measured["host_samples.jsonl"]["sha256"]:
        raise ValueError("Raw host evidence differs from measurement")
    rows = [json.loads(line) for line in raw_path.read_text().splitlines()]
    diagnostics = unmatched_inventory(rows)
    output = io.StringIO()
    monitor = HostMonitor(output, measured["baseline_scope"])
    monitor.observer_pid = measured["wrapper_pid"]
    for index, row in enumerate(rows):
        if row["index"] != index or row["observer_pid"] != measured["wrapper_pid"]:
            raise ValueError("Host observer/index differs")
        monitor.sample_fn = lambda sample=row["snapshot"]: sample
        monitor.observe()
    if [json.loads(line) for line in output.getvalue().splitlines()] != rows:
        raise ValueError("Retained host intervals do not replay")
    summary = monitor.summary(measured["command_launch_started_monotonic_s"], measured["command_wait_finished_monotonic_s"])
    retained = measured["host_workload"]
    for key in ("source", "observer_source"):
        if summary[key]["sha256"] != retained[key]["sha256"]:
            raise ValueError("Host observer source changed")
    if {k: v for k, v in summary.items() if k not in ("source", "observer_source")} != {
            k: v for k, v in retained.items() if k not in ("source", "observer_source")}:
        raise ValueError("Retained host summary does not replay")
    if [record(measurement_path), record(raw_path)] != identities:
        raise ValueError("Evidence changed during review")
    return {"status": "host_interval_causes_reviewed_not_admitted", "source": record(__file__),
            "evidence": identities, "job_id": measured["job_id"], "clock_domain": measured["clock_domain"],
            "native_host_evidence_path": measured["host_samples.jsonl"]["path"],
            "retained_host_summary": retained, "diagnostics": diagnostics,
            "scientific_timing_admitted": False, "controlled_workload_verified": False,
            "limitations": ["Process names are descriptive, not independently authenticated kernel-thread classifications.",
                            "Unmatched process CPU use and short-lived unsampled work remain unmeasured.",
                            "No monitor threshold, classification or deployed recipe was changed.",
                            "This host-only replay does not validate native output, resource samples or complete timing admission."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = review(args.directory)
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
