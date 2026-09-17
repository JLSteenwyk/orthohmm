"""Measure one fresh command inside its existing dedicated Slurm task subtree."""

import argparse
from contextlib import ExitStack
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import psutil

sys.path.insert(0, str(Path(__file__).resolve().parent))
import slurm_resource_snapshot
from monitor_slurm_resources import source_record, summarize
from command_host_monitor import HostMonitor


def check_allocation(sample, pid, cpu_count, memory_bytes):
    if {p["pid"] for p in sample["processes"]} != {pid} or sample["sampling_errors"]:
        raise ValueError("Require a dedicated task subtree containing only the measurement wrapper before launch")
    if len(sample["metrics"]["effective_cpus"]) != cpu_count:
        raise ValueError("Effective cpuset differs from requested CPU allocation")
    limits = [int(item["memory_max"]) for item in sample["ancestor_limits"] if item["memory_max"] != "max"]
    if not limits or min(limits) != memory_bytes:
        raise ValueError("Inherited task/job memory cap differs from requested allocation")


def stop_owned_group(process, grace=5.):
    # Only the new session created by this wrapper is targeted, never the Slurm job.
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return process.wait()
    deadline = time.monotonic() + grace
    while time.monotonic() < deadline:
        process.poll()
        try:
            os.killpg(process.pid, 0)
        except ProcessLookupError:
            return process.wait()
        time.sleep(.02)
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    return process.wait()


def measure(command, output, job_id, cpu_count, memory_bytes, timeout_s, interval_s=1., snapshot_fn=None,
            monitor_host=False, host_snapshot_fn=None):
    if (not command or cpu_count < 1 or memory_bytes < 1 or not math.isfinite(timeout_s)
            or timeout_s <= 0 or not math.isfinite(interval_s) or interval_s <= 0):
        raise ValueError("Invalid command or resource measurement plan")
    if output.exists():
        raise FileExistsError(output)
    snapshot_fn = snapshot_fn or slurm_resource_snapshot.snapshot
    output.mkdir(parents=True)
    pid = os.getpid()
    anchor = psutil.Process(pid).create_time()
    started = time.monotonic()
    report = {"status": "preflight", "command": command, "cwd": str(Path.cwd()), "job_id": job_id,
              "requested_cpus": cpu_count, "requested_memory_bytes": memory_bytes, "timeout_s": timeout_s,
              "interval_s": interval_s, "source": source_record(__file__),
              "snapshot_source": source_record(slurm_resource_snapshot.__file__), "accuracy_evaluated": False,
              "controlled_workload_verified": False, "limitations": [
                  "Caller must independently verify source/runtime, inputs, host workload and native output correctness.",
                  "Cgroup CPU and memory include this wrapper; sampling and accounting overhead are retained.",
                  "Cgroup peak may include preflight and predate command launch; sampled RSS misses transient peaks and may double-count shared pages.",
                  "Command wall time includes spawn/wait scheduling overhead; separate conversion/scoring must use separate measurements.",
                  "A command exit does not certify biological output success. Measurement failure does not silently terminate a still-running command."]}
    process, samples, timed_out = None, [], False
    try:
        with ExitStack() as stack:
            series = stack.enter_context((output / "samples.jsonl").open("x"))
            log = stack.enter_context((output / "command.log").open("x"))
            def observe():
                snapshot = snapshot_fn(pid, job_id)
                row = {"index": len(samples), "elapsed_s": time.monotonic() - started, "anchor_created": anchor,
                       "snapshot": snapshot, "host": {"load_average": list(os.getloadavg())}}
                samples.append(row)
                series.write(json.dumps(row, sort_keys=True) + "\n")
                series.flush()
                report["summary"] = summarize(samples)
                return snapshot
            baseline = observe()
            check_allocation(baseline, pid, cpu_count, memory_bytes)
            report["baseline_scope"] = baseline["scope"]
            host = None
            if monitor_host:
                handle = stack.enter_context((output / "host_samples.jsonl").open("x"))
                host = HostMonitor(handle, baseline["scope"], host_snapshot_fn)
                host.observe()
            launch = time.monotonic()
            process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            report.update(status="running", child_pid=process.pid)
            while True:
                remaining = timeout_s - (time.monotonic() - launch)
                if remaining <= 0:
                    timed_out = True
                    stop_owned_group(process)
                    break
                try:
                    process.wait(timeout=min(interval_s, remaining))
                    break
                except subprocess.TimeoutExpired:
                    if host is not None:
                        host.observe()
                    if "measurement_error" not in report:
                        try:
                            observe()
                        except Exception as error:
                            report["measurement_error"] = {"type": type(error).__name__, "message": str(error)}
            report.update(exit_code=process.returncode, command_wall_s=time.monotonic() - launch, timed_out=timed_out)
            end = time.monotonic()
            if host is not None:
                host.observe()
                report["host_workload"] = host.summary(launch, end)
            if "measurement_error" not in report:
                try:
                    final = observe()
                    leftovers = [p for p in final["processes"] if p["pid"] != pid]
                    if leftovers:
                        report["remaining_processes"] = leftovers
                        report["measurement_error"] = {"type": "RemainingProcesses", "message": "Command left processes in the dedicated task subtree"}
                        stop_owned_group(process)
                except Exception as error:
                    report["measurement_error"] = {"type": type(error).__name__, "message": str(error)}
            report["status"] = ("measurement_failed" if "measurement_error" in report else
                                "command_timed_out" if timed_out else "command_failed" if process.returncode else "command_exited_zero")
    except BaseException as error:
        report.update(status="wrapper_failed", error_type=type(error).__name__, error=str(error))
        if process is not None and process.poll() is None:
            stop_owned_group(process)
            report["exit_code"] = process.returncode
        raise
    finally:
        report["wrapper_wall_s"] = time.monotonic() - started
        for name in ("samples.jsonl", "command.log", "host_samples.jsonl"):
            if (output / name).exists():
                report[name] = source_record(output / name)
        (output / "results.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--job-id", type=int, required=True)
    parser.add_argument("--cpus", type=int, required=True)
    parser.add_argument("--memory-gib", type=int, required=True)
    parser.add_argument("--timeout", type=float, required=True)
    parser.add_argument("--interval", type=float, default=1.)
    parser.add_argument("--monitor-host", action="store_true")
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    result = measure(command, args.output.resolve(), args.job_id, args.cpus, args.memory_gib * 1024 ** 3, args.timeout,
                     args.interval, monitor_host=args.monitor_host)
    raise SystemExit(0 if result["status"] == "command_exited_zero" else 1)
