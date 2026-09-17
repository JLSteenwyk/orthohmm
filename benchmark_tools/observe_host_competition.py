"""Read-only CPU competition evidence; absence of detection is not host exclusivity."""

import argparse
import json
import math
import os
from pathlib import Path, PurePosixPath
import sys
import time

import psutil

sys.path.insert(0, str(Path(__file__).resolve().parent))
from slurm_resource_snapshot import scoped_path
from monitor_slurm_resources import source_record


def membership(pid):
    lines = Path(f"/proc/{pid}/cgroup").read_text().splitlines()
    rows = [line[3:] for line in lines if line.startswith("0::")]
    if len(rows) != 1:
        raise ValueError("Require one unified cgroup")
    path = PurePosixPath(rows[0])
    if not path.is_absolute() or ".." in path.parts:
        raise ValueError("Invalid cgroup membership")
    return str(path)


def snapshot():
    started = time.monotonic()
    rows, errors = [], []
    for pid in psutil.pids():
        try:
            process = psutil.Process(pid)
            created = process.create_time()
            group = membership(pid)
            cpu = process.cpu_times()
            observed = time.monotonic()
            name = process.name()
            if psutil.Process(pid).create_time() != created or membership(pid) != group:
                raise ValueError("Process identity or cgroup changed during collection")
            rows.append({"pid": pid, "created": created, "cgroup": group, "name": name,
                         "user_s": cpu.user, "system_s": cpu.system, "observed_monotonic_s": observed})
        except (psutil.Error, OSError, ValueError) as error:
            errors.append({"pid": pid, "type": type(error).__name__})
    return {"processes": rows, "errors": errors, "started_monotonic_s": started,
            "finished_monotonic_s": time.monotonic(), "psutil_version": psutil.__version__}


def analyze(before, after, scope, monitor_pid, threshold=.25):
    scope = PurePosixPath(scope)
    if (not scope.is_absolute() or scope == PurePosixPath("/") or ".." in scope.parts
            or not math.isfinite(threshold) or threshold <= 0):
        raise ValueError("Require a scoped subtree and positive CPU threshold")
    def index(sample):
        rows = sample["processes"]
        keyed = {(r["pid"], r["created"]): r for r in rows}
        if len(keyed) != len(rows) or len({r["pid"] for r in rows}) != len(rows):
            raise ValueError("Duplicate process identity in snapshot")
        return keyed
    a, b = index(before), index(after)
    foreign = lambda row: row["pid"] != monitor_pid and not PurePosixPath(row["cgroup"]).is_relative_to(scope)
    unmatched = [{"pid": pid, "created": created} for pid, created in sorted(set(a) ^ set(b))
                 if foreign((a if (pid, created) in a else b)[pid, created])]
    uncertain, measured = [], []
    for key in sorted(set(a) & set(b)):
        first, last = a[key], b[key]
        if not foreign(first) and not foreign(last):
            continue
        if first["cgroup"] != last["cgroup"]:
            uncertain.append({"pid": key[0], "reason": "cgroup_changed"})
            continue
        span = last["observed_monotonic_s"] - first["observed_monotonic_s"]
        deltas = [last[k] - first[k] for k in ("user_s", "system_s")]
        if not math.isfinite(span) or span <= 0 or any(not math.isfinite(x) or x < 0 for x in deltas):
            uncertain.append({"pid": key[0], "reason": "invalid_cpu_counter_or_time"})
            continue
        measured.append({"pid": key[0], "created": key[1], "name": last["name"], "cgroup": last["cgroup"],
                         "cpu_s": sum(deltas), "span_s": span, "average_cores": sum(deltas) / span})
    total = sum(r["average_cores"] for r in measured)
    error_count = len(before["errors"]) + len(after["errors"])
    state = ("competing_cpu_observed" if total >= threshold else "inconclusive"
             if unmatched or uncertain or error_count else "no_large_persistent_competitor_observed")
    return {"status": state, "controlled_workload_verified": False, "scope": str(scope),
            "threshold_average_cores": threshold, "sum_observed_foreign_average_cores": total,
            "persistent_foreign_processes": sorted(measured, key=lambda r: (-r["average_cores"], r["pid"])),
            "unmatched_foreign_processes": unmatched, "uncertain_processes": uncertain, "sampling_error_count": error_count}


def observe(pid, job_id, interval, output):
    if output.exists():
        raise FileExistsError(output)
    if not math.isfinite(interval) or interval <= 0:
        raise ValueError("Require a positive observation interval")
    anchor = psutil.Process(pid).create_time()
    scope = scoped_path(Path(f"/proc/{pid}/cgroup").read_text(), job_id)
    first = snapshot()
    time.sleep(interval)
    second = snapshot()
    if psutil.Process(pid).create_time() != anchor or membership(pid) != str(scope):
        raise ValueError("Observation anchor changed identity or cgroup")
    report = {"source": source_record(__file__), "anchor_pid": pid, "anchor_created": anchor, "job_id": job_id,
              "monitor_pid": os.getpid(), "requested_interval_s": interval, "snapshots": [first, second],
              "summary": analyze(first, second, str(scope), os.getpid()), "accuracy_evaluated": False,
              "limitations": ["Visible persistent processes only; short-lived work between samples can be missed.",
                  "CPU deltas use user+system only; child accounting is excluded to avoid counting descendants twice.",
                  "Per-process intervals differ slightly because collection is not atomic.",
                  "Does not establish CPU exclusivity, quiet I/O, memory bandwidth, GPU contention, or a complete run's workload.",
                  "No command lines, environments or signals to other processes are collected or issued."]}
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pid", type=int, required=True)
    parser.add_argument("--job-id", type=int, required=True)
    parser.add_argument("--interval", type=float, default=3.)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    observe(args.pid, args.job_id, args.interval, args.output.resolve())
