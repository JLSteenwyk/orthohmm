"""Replay a completed collector record; this does not certify benchmark admission."""

import argparse
import io
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from command_host_monitor import HostMonitor
from measure_slurm_command import check_allocation
from monitor_slurm_resources import source_record, summarize


def require(condition, message):
    if not condition:
        raise ValueError(message)


def audit(directory, expected_sha):
    directory = Path(directory).resolve()
    report_path = directory / "results.json"
    report_record = source_record(report_path)
    require(report_record["sha256"] == expected_sha, "Collector report hash differs")
    report = json.loads(report_path.read_text())
    require(report["status"] == "command_exited_zero" and report["exit_code"] == 0
            and report["timed_out"] is False and "measurement_error" not in report,
            "Require completed error-free command measurement")
    raw_records = []
    for name in ("samples.jsonl", "host_samples.jsonl", "command.log"):
        item = source_record(directory / name)
        require(item["sha256"] == report[name]["sha256"], "Raw evidence hash differs: " + name)
        raw_records.append(item)
    for key, filename in (("source", "measure_slurm_command.py"),
                          ("snapshot_source", "slurm_resource_snapshot.py")):
        require(report[key]["sha256"] == source_record(Path(__file__).with_name(filename))["sha256"],
                "Unsupported collector source revision")
    times = [report[k] for k in ("wrapper_started_monotonic_s", "command_launch_started_monotonic_s",
                                "command_wait_finished_monotonic_s", "wrapper_finished_monotonic_s")]
    require(all(isinstance(t, (int, float)) and math.isfinite(t) for t in times), "Invalid clock bounds")
    start, launch, end, finish = times
    require(start <= launch < end <= finish, "Invalid command boundary order")
    require(report["command_wall_s"] == end - launch and report["wrapper_wall_s"] == finish - start,
            "Wall duration differs from recorded bounds")
    domain = report["clock_domain"]
    require(domain["hostname"] and domain["boot_id"] and domain["clock"] == "time.monotonic"
            and domain["unit"] == "seconds", "Missing clock-domain provenance")
    samples = [json.loads(line) for line in (directory / "samples.jsonl").read_text().splitlines()]
    require(len(samples) >= 2, "Insufficient resource observations")
    previous_end = start
    baseline_cpus = samples[0]["snapshot"]["metrics"]["effective_cpus"]
    for index, row in enumerate(samples):
        a, b = row["started_monotonic_s"], row["finished_monotonic_s"]
        require(row["index"] == index and previous_end <= a <= b <= finish,
                "Resource observation bounds or indices differ")
        require(row["elapsed_s"] == b - start, "Relative resource time differs")
        previous_end = b
        snapshot = row["snapshot"]
        require(snapshot["job_id"] == report["job_id"] and snapshot["pid"] == report["wrapper_pid"]
                and snapshot["scope"] == report["baseline_scope"], "Resource scope differs")
        require(not snapshot["sampling_errors"], "Resource observation has sampling errors")
        require(snapshot["metrics"]["effective_cpus"] == baseline_cpus, "Effective cpuset changed")
        limits = [int(v["memory_max"]) for v in snapshot["ancestor_limits"] if v["memory_max"] != "max"]
        require(limits and min(limits) == report["requested_memory_bytes"], "Memory cap changed")
    check_allocation(samples[0]["snapshot"], report["wrapper_pid"],
                     report["requested_cpus"], report["requested_memory_bytes"])
    check_allocation(samples[-1]["snapshot"], report["wrapper_pid"],
                     report["requested_cpus"], report["requested_memory_bytes"])
    require(samples[0]["finished_monotonic_s"] <= launch and samples[-1]["started_monotonic_s"] >= end,
            "Resource observations do not bracket command")
    require(summarize(samples) == report["summary"], "Resource summary does not replay")
    rows = [json.loads(line) for line in (directory / "host_samples.jsonl").read_text().splitlines()]
    require(len(rows) >= 2, "Insufficient host observations")
    output = io.StringIO()
    monitor = HostMonitor(output, report["baseline_scope"])
    monitor.observer_pid = report["wrapper_pid"]
    previous_end = start
    for index, row in enumerate(rows):
        require("observation_error" not in row, "Host observation error retained")
        require(row["index"] == index and row["observer_pid"] == report["wrapper_pid"],
                "Host observer identity or indices differ")
        snapshot = row["snapshot"]
        require(previous_end <= snapshot["started_monotonic_s"] <= snapshot["finished_monotonic_s"] <= finish,
                "Host observation bounds differ")
        previous_end = snapshot["finished_monotonic_s"]
        monitor.sample_fn = lambda: snapshot
        monitor.observe()
    replayed = [json.loads(line) for line in output.getvalue().splitlines()]
    require(replayed == rows, "Host intervals do not replay")
    rebuilt, retained = monitor.summary(launch, end), report["host_workload"]
    for key in ("source", "observer_source"):
        require(rebuilt[key]["sha256"] == retained[key]["sha256"], "Unsupported host observer source")
    require({k: v for k, v in rebuilt.items() if k not in {"source", "observer_source"}} ==
            {k: v for k, v in retained.items() if k not in {"source", "observer_source"}},
            "Host summary does not replay")
    require(rebuilt["command_bracketed_by_samples"], "Host observations do not bracket command")
    require(source_record(report_path) == report_record, "Report changed during audit")
    for item in raw_records:
        require(source_record(item["path"]) == item, "Evidence changed during audit")
    return {"status": "collector_evidence_replayed", "report": report_record,
            "source": source_record(__file__), "raw_evidence": raw_records,
            "resource_observations": len(samples), "host_observations": len(rows),
            "host_status": rebuilt["status"], "controlled_workload_verified": False,
            "benchmark_admitted": False,
            "limitations": ["Does not certify source/input/native output correctness, scheduler accounting or host exclusivity.",
                            "Observed resource span includes collector overhead; sampling gaps remain.",
                            "Clock domain is retained provenance, not independent authentication of the originating host."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.directory, args.report_sha256)
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
