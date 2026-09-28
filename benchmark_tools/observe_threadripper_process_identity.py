"""Capture identity-checked Linux Kthread fields without approving background work."""

import argparse
import json
import math
import os
from pathlib import Path
import time

import psutil

from benchmark_tools.observe_host_competition import membership, snapshot
from benchmark_tools.monitor_slurm_resources import source_record


def kernel_fields(text, pid):
    fields = {}
    for line in text.splitlines():
        key, sep, value = line.partition(":")
        if sep and key in {"Pid", "Tgid", "Kthread"}:
            if key in fields:
                raise ValueError("Duplicate proc status identity field")
            fields[key] = value.strip()
    if (fields.get("Pid") != str(pid) or fields.get("Tgid") != str(pid)
            or fields.get("Kthread") not in {"0", "1"}):
        raise ValueError("Missing or inconsistent proc status identity fields")
    return dict(pid=pid, tgid=pid, kthread=int(fields["Kthread"]))


def details(row):
    pid = row["pid"]
    started = time.monotonic()
    process = psutil.Process(pid)
    if process.create_time() != row["created"] or membership(pid) != row["cgroup"]:
        raise ValueError("Process identity changed before status observation")
    path = Path(f"/proc/{pid}/status")
    first = kernel_fields(path.read_text(), pid)
    last = kernel_fields(path.read_text(), pid)
    # Fresh Process objects avoid accepting the first object's cached create time.
    final = psutil.Process(pid)
    if (first != last or final.create_time() != row["created"]
            or membership(pid) != row["cgroup"] or not process.is_running()
            or (first["kthread"] == 0 and final.name() != row["name"])):
        raise ValueError("Process identity or kernel type changed during status observation")
    return dict(first, started_monotonic_s=started, finished_monotonic_s=time.monotonic())


def enriched_snapshot():
    boot_path = Path("/proc/sys/kernel/random/boot_id")
    boot = boot_path.read_text().strip()
    sample = snapshot()
    for row in sample["processes"]:
        try:
            row["kernel_identity"] = details(row)
        except (psutil.Error, OSError, ValueError) as error:
            # Preserve the original CPU row; missing type evidence remains an error.
            row["kernel_identity_error"] = type(error).__name__
            sample["errors"].append(dict(pid=row["pid"], type=type(error).__name__,
                                         stage="kernel_identity"))
    if boot_path.read_text().strip() != boot:
        raise ValueError("Boot changed during observation")
    sample.update(boot_id=boot, schema="threadripper_typed_process_snapshot_v1",
                  finished_monotonic_s=time.monotonic())
    return sample


def capture(output, interval=3.):
    if not math.isfinite(interval) or interval <= 0:
        raise ValueError("Require a positive observation interval")
    # Reserve the output before doing work, so failed captures cannot be overwritten.
    with Path(output).open("x") as handle:
        report = dict(schema="threadripper_typed_process_observation_v1", status="started",
                      observer_pid=os.getpid(), requested_interval_s=interval,
                      source=source_record(__file__), snapshots=[],
                      observer_source=source_record(Path(__file__).with_name("observe_host_competition.py")),
                      host=os.uname().nodename, kernel_release=os.uname().release,
                      psutil_version=psutil.__version__,
                      controlled_workload_verified=False, scientific_timings_admitted=False)
        try:
            report["snapshots"].append(enriched_snapshot())
            time.sleep(interval)
            report["snapshots"].append(enriched_snapshot())
            if report["snapshots"][0]["boot_id"] != report["snapshots"][1]["boot_id"]:
                raise ValueError("Snapshots belong to different boots")
            report["status"] = "observations_retained"
        except BaseException as error:
            report.update(status="capture_failed", error_type=type(error).__name__)
            raise
        finally:
            json.dump(report, handle, indent=2, sort_keys=True)
            handle.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--interval", type=float, default=3.)
    args = parser.parse_args()
    capture(args.output, args.interval)
