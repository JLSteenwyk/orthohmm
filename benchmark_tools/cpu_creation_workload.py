"""Bounded Linux CPU/process-creation diagnostic; no timing admission."""

import argparse
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import time


def validate(mode, seconds, child_cpu_s, max_children):
    if mode not in ("steady", "creation"):
        raise ValueError("Unknown workload mode")
    if (type(seconds) not in (int, float) or not math.isfinite(seconds) or not 0 < seconds <= 30
            or type(child_cpu_s) not in (int, float) or not math.isfinite(child_cpu_s)
            or not .0001 <= child_cpu_s <= .1
            or type(max_children) is not int or not 1 <= max_children <= 10000):
        raise ValueError("Workload exceeds fixed safety bounds")


def burn_until(deadline, cpu_seconds=None):
    started = time.process_time()
    value = 1
    while time.monotonic() < deadline:
        if cpu_seconds is not None and time.process_time() - started >= cpu_seconds:
            return True
        for _ in range(256):
            value = (value * 1664525 + 1013904223) & 0xffffffff
    return cpu_seconds is None


def worker(mode, cpu, seconds, child_cpu_s, max_children):
    validate(mode, seconds, child_cpu_s, max_children)
    allowed = sorted(os.sched_getaffinity(0))
    if type(cpu) is not int or cpu not in allowed:
        raise ValueError("Requested CPU outside inherited affinity")
    os.sched_setaffinity(0, {cpu})
    membership = Path("/proc/self/cgroup").read_text()
    started_ns, started_cpu = time.monotonic_ns(), time.process_time()
    deadline = time.monotonic() + seconds
    child_count, child_cpu, child_max_wall, limited_children = 0, 0., 0., 0
    capped = False
    if mode == "steady":
        burn_until(deadline)
    else:
        while time.monotonic() < deadline:
            if child_count == max_children:
                capped = True
                break
            start = time.monotonic()
            pid = os.fork()
            if pid == 0:
                try:
                    met = burn_until(min(deadline, start + 1.), child_cpu_s)
                    os._exit(0 if met else 3)
                except BaseException:
                    os._exit(4)
            _, status, usage = os.wait4(pid, 0)
            code = os.waitstatus_to_exitcode(status)
            if code not in (0, 3):
                raise RuntimeError(f"Creation child failed: {code}")
            child_count += 1
            limited_children += code == 3
            child_cpu += usage.ru_utime + usage.ru_stime
            child_max_wall = max(child_max_wall, time.monotonic() - start)
    result = dict(mode=mode, pid=os.getpid(), cpu=cpu, inherited_affinity=allowed,
        affinity=sorted(os.sched_getaffinity(0)), membership=membership,
        started_ns=started_ns, finished_ns=time.monotonic_ns(), self_cpu_s=time.process_time()-started_cpu,
        children_cpu_s=child_cpu, reaped_children=child_count, wall_limited_children=limited_children,
        maximum_child_wall_s=child_max_wall, child_count_cap_reached=capped,
        requested_seconds=seconds, child_cpu_target_s=child_cpu_s, max_children=max_children)
    if result["affinity"] != [cpu] or Path("/proc/self/cgroup").read_text() != membership:
        raise ValueError("Workload affinity or membership changed")
    return result


def run(mode, cpus, seconds, child_cpu_s=.005, max_children=10000):
    validate(mode, seconds, child_cpu_s, max_children)
    allowed = set(os.sched_getaffinity(0))
    if (not cpus or len(cpus) > 20 or len(set(cpus)) != len(cpus)
            or any(type(cpu) is not int for cpu in cpus) or not set(cpus) <= allowed):
        raise ValueError("Require 1-20 distinct allowed CPUs")
    processes, rows = [], []
    source = str(Path(__file__).resolve())
    started = time.monotonic_ns()
    deadline = time.monotonic() + seconds + 15
    try:
        for cpu in cpus:
            command = [sys.executable, "-I", "-B", source, "--worker", "--mode", mode,
                "--cpus", str(cpu), "--seconds", str(seconds), "--child-cpu", str(child_cpu_s),
                "--max-children", str(max_children)]
            processes.append(subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                              text=True, start_new_session=True))
        for cpu, process in zip(cpus, processes):
            stdout, stderr = process.communicate(timeout=max(.001, deadline - time.monotonic()))
            if process.returncode:
                raise RuntimeError(f"Worker failed on CPU {cpu}: {process.returncode}: {stderr}")
            row = json.loads(stdout)
            if row["pid"] != process.pid or row["cpu"] != cpu or row["mode"] != mode:
                raise ValueError("Worker identity mismatch")
            rows.append(row)
    finally:
        for process in processes:
            # Each worker owns its session; kill its group to include any live fork child.
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait()
            for stream in (process.stdout, process.stderr):
                if stream is not None:
                    stream.close()
    return dict(status="bounded_workload_completed", mode=mode, cpus=cpus,
        started_ns=started, finished_ns=time.monotonic_ns(), workers=rows,
        reaped_children=sum(r["reaped_children"] for r in rows),
        worker_cpu_s=sum(r["self_cpu_s"] + r["children_cpu_s"] for r in rows),
        child_count_cap_reached=any(r["child_count_cap_reached"] for r in rows),
        scientific_timings_admitted=False, controlled_workload_verified=False, publication_ready=False,
        limitations=["CPU includes workers and reaped fork children, not launcher/startup or kernel-wide activity.",
                     "Workers start asynchronously; utilization is observed, not assumed equal across modes.",
                     "Fork workload is diagnostic and does not reproduce native tool process execution."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--mode", choices=("steady", "creation"), required=True)
    parser.add_argument("--cpus", type=int, nargs="+", required=True)
    parser.add_argument("--seconds", type=float, required=True)
    parser.add_argument("--child-cpu", type=float, default=.005)
    parser.add_argument("--max-children", type=int, default=10000)
    args = parser.parse_args()
    if args.worker:
        if len(args.cpus) != 1:
            parser.error("Worker requires exactly one CPU")
        result = worker(args.mode, args.cpus[0], args.seconds, args.child_cpu, args.max_children)
    else:
        result = run(args.mode, args.cpus, args.seconds, args.child_cpu, args.max_children)
    print(json.dumps(result, sort_keys=True, allow_nan=False))
