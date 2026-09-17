"""Retain a bounded read-only resource time series for one live Slurm subtree."""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time

import psutil

sys.path.insert(0, str(Path(__file__).resolve().parent))
import slurm_resource_snapshot


def summarize(samples):
    if not samples:
        raise ValueError("No successful resource observations")
    first, last = samples[0], samples[-1]
    identity = (first["snapshot"]["job_id"], first["snapshot"]["pid"], first["snapshot"]["scope"], first["anchor_created"])
    for index, sample in enumerate(samples):
        current = sample["snapshot"]
        if (current["job_id"], current["pid"], current["scope"], sample["anchor_created"]) != identity:
            raise ValueError("Resource observations changed job, anchor identity or scope")
        if index:
            previous = samples[index - 1]
            if sample["elapsed_s"] <= previous["elapsed_s"]:
                raise ValueError("Non-increasing sample time")
            for key in ("usage_usec", "user_usec", "system_usec"):
                if current["metrics"]["cpu"][key] < previous["snapshot"]["metrics"]["cpu"][key]:
                    raise ValueError("Cgroup CPU counter decreased")
    delta = {key: last["snapshot"]["metrics"]["cpu"][key] - first["snapshot"]["metrics"]["cpu"][key]
             for key in ("usage_usec", "user_usec", "system_usec")}
    interval = last["elapsed_s"] - first["elapsed_s"]
    return {"observations": len(samples), "observation_span_s": interval, "cpu_delta_usec": delta,
            "mean_cpu_cores_over_observed_span": delta["usage_usec"] / 1e6 / interval if interval else None,
            "maximum_sampled_sum_rss_bytes": max(s["snapshot"]["sampled_sum_process_rss_bytes"] for s in samples),
            "maximum_reported_cgroup_peak_bytes": max(s["snapshot"]["metrics"]["memory_peak_since_creation_or_reset_bytes"] for s in samples),
            "samples_with_process_errors": sum(bool(s["snapshot"]["sampling_errors"]) for s in samples)}


def source_record(path):
    path = Path(path).resolve()
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def collect(pid, job_id, output, count, interval, sample_fn=None, anchor_fn=None, host_fn=None,
            clock=time.monotonic, sleep=time.sleep):
    if count < 2 or not math.isfinite(interval) or interval <= 0:
        raise ValueError("Require at least two samples and a finite positive interval")
    if output.exists():
        raise FileExistsError(output)
    sample_fn = sample_fn or slurm_resource_snapshot.snapshot
    anchor_fn = anchor_fn or (lambda value: psutil.Process(value).create_time())
    host_fn = host_fn or (lambda: {"load_average": list(os.getloadavg()),
                                  "per_cpu_times": [value._asdict() for value in psutil.cpu_times(percpu=True)]})
    output.mkdir(parents=True)
    samples = []
    started = clock()
    report = {"status": "collecting", "job_id": job_id, "pid": pid, "requested_samples": count,
              "interval_s": interval, "source": source_record(__file__), "snapshot_source": source_record(slurm_resource_snapshot.__file__),
              "psutil_version": psutil.__version__, "accuracy_evaluated": False,
              "limitations": ["Bounded observation window, not evidence that the job or inference completed.",
                  "A sampling failure does not establish job termination; recheck authoritative scheduler state.",
                  "CPU deltas cover only the recorded span and cgroup subtree, including any launcher or validation activity there.",
                  "Reported cgroup peak may precede this window; sampled RSS misses between-sample peaks and can double-count shared pages.",
                  "Host load/CPU counters include the target job and do not isolate unrelated workloads or certify quiet-host conditions."]}
    try:
        with (output / "samples.jsonl").open("x") as handle:
            for index in range(count):
                if index:
                    sleep(max(0., started + index * interval - clock()))
                anchor = anchor_fn(pid)
                observed = sample_fn(pid, job_id)
                if anchor_fn(pid) != anchor:
                    raise ValueError("Anchor PID was reused during observation")
                row = {"index": index, "elapsed_s": clock() - started, "anchor_created": anchor,
                       "snapshot": observed, "host": host_fn()}
                samples.append(row)
                handle.write(json.dumps(row, sort_keys=True) + "\n")
                handle.flush()
                report["summary"] = summarize(samples)
        report["status"] = "bounded_observations_complete"
    except BaseException as error:
        report.update(status="observation_failed", error_type=type(error).__name__, error=str(error), observations_retained=len(samples))
        raise
    finally:
        report["collector_wall_s"] = clock() - started
        report["samples"] = source_record(output / "samples.jsonl")
        (output / "results.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pid", required=True, type=int)
    parser.add_argument("--job-id", required=True, type=int)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--samples", type=int, default=5)
    parser.add_argument("--interval", type=float, default=5.)
    args = parser.parse_args()
    collect(args.pid, args.job_id, args.output.resolve(), args.samples, args.interval)
