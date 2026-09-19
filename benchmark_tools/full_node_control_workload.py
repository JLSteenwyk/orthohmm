"""Bounded steady/process-creation workloads with per-process witnesses."""

import argparse
import math
import os
from pathlib import Path
import resource
import signal
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.probe_dgx_step_separation import save, wait_file

DURATION = 20.
MAX_CREATIONS = 200000


def require_start(value):
    if type(value) is not dict or set(value) != {"go"} or value["go"] is not True:
        raise ValueError("Invalid shared start signal")


def validate(mode, cpus, duration, cap):
    if mode not in {"steady", "churn"}:
        raise ValueError("Unknown bounded workload")
    if (not cpus or len(cpus) > 20 or len(set(cpus)) != len(cpus)
            or any(type(cpu) is not int or cpu < 0 for cpu in cpus)
            or not set(cpus) <= os.sched_getaffinity(0)):
        raise ValueError("Invalid bounded worker affinities")
    if type(duration) not in (int, float) or not math.isfinite(duration) or not 0 < duration <= DURATION:
        raise ValueError("Invalid bounded duration")
    if type(cap) is not int or not 0 < cap <= MAX_CREATIONS:
        raise ValueError("Invalid creation cap")


def work(mode, duration, cap, parent_pid):
    started = time.monotonic_ns()
    cpu_start = time.process_time()
    children_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    deadline = started + int(duration * 1e9)
    created, value = 0, 1
    while time.monotonic_ns() < deadline:
        if os.getppid() != parent_pid:
            raise RuntimeError("Workload parent disappeared")
        if mode == "steady":
            for _ in range(1000):
                value = (value * 1664525 + 1013904223) % 4294967296
        else:
            if created == cap:
                break
            child = os.fork()
            if child == 0:
                os._exit(0)
            waited, status = os.waitpid(child, 0)
            if waited != child or status != 0:
                raise RuntimeError(f"Short-lived child failed: pid={child}, waited={waited}, status={status}")
            created += 1
    finished = time.monotonic_ns()
    child_cpu = resource.getrusage(resource.RUSAGE_CHILDREN)
    return dict(started_ns=started, finished_ns=finished, self_cpu_s=time.process_time()-cpu_start,
        waited_child_user_s=child_cpu.ru_utime-children_before.ru_utime,
        waited_child_system_s=child_cpu.ru_stime-children_before.ru_stime,
        creations=created, creation_cap_reached=mode == "churn" and created == cap,
        checksum=value if mode == "steady" else None)


def child_worker(directory, cpu, mode, duration, cap, parent_pid):
    allowed = sorted(os.sched_getaffinity(0))
    os.sched_setaffinity(0, {cpu})
    identity = dict(pid=os.getpid(), parent_pid=parent_pid, cpu=cpu, allowed=allowed,
                    affinity=sorted(os.sched_getaffinity(0)), membership=Path("/proc/self/cgroup").read_text())
    save(directory / f"ready_{cpu}.json", identity)
    go = wait_file(directory / "workload_go.json", seconds=30)
    require_start(go)
    result = dict(identity, mode=mode, duration_s=duration, creation_cap=cap,
                  **work(mode, duration, cap, parent_pid))
    result["final_affinity"] = sorted(os.sched_getaffinity(0))
    result["final_membership"] = Path("/proc/self/cgroup").read_text()
    save(directory / f"done_{cpu}.json", result)


def workload(directory, mode, cpus, duration=DURATION, cap=MAX_CREATIONS):
    validate(mode, cpus, duration, cap)
    if (directory / "workload_ready.json").exists() or any(directory.glob("ready_*.json")):
        raise FileExistsError("Require fresh workload witnesses")
    parent = os.getpid()
    living, statuses = set(), {}
    started = time.monotonic_ns()
    try:
        for cpu in cpus:
            child = os.fork()
            if child == 0:
                code = 0
                try:
                    child_worker(directory, cpu, mode, duration, cap, parent)
                except BaseException as error:
                    code = 1
                    try:
                        save(directory / f"failed_{cpu}.json", dict(error_type=type(error).__name__, error=str(error)))
                    except OSError:
                        pass
                finally:
                    # A child must never execute the parent's sibling-cleanup path.
                    os._exit(code)
            living.add(child)
        ready = [wait_file(directory / f"ready_{cpu}.json", seconds=30) for cpu in cpus]
        save(directory / "workload_ready.json", dict(pid=parent, cpus=cpus, mode=mode,
             membership=Path("/proc/self/cgroup").read_text(), workers=ready))
        deadline = time.monotonic() + duration + 35
        while living:
            for pid in list(living):
                waited, status = os.waitpid(pid, os.WNOHANG)
                if waited:
                    living.remove(pid)
                    statuses[str(pid)] = status
                    if status != 0:
                        raise RuntimeError("Bounded worker failed")
            if time.monotonic() >= deadline:
                raise TimeoutError("Bounded worker deadline exceeded")
            if living:
                time.sleep(.01)
        workers = [wait_file(directory / f"done_{cpu}.json", seconds=1) for cpu in cpus]
        result = dict(mode=mode, cpus=cpus, duration_s=duration, creation_cap=cap, workers=workers,
                      statuses=statuses, started_ns=started, finished_ns=time.monotonic_ns(),
                      pid=parent, membership=Path("/proc/self/cgroup").read_text(),
                      scientific_timings_admitted=False)
        save(directory / "workload_done.json", result)
        if any(row["creation_cap_reached"] for row in workers):
            raise RuntimeError("Creation cap reached; workload is invalid")
        return result
    finally:
        for pid in living:
            try:
                os.kill(pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
        for pid in living:
            os.waitpid(pid, 0)


def competitor(directory, cpu, duration=DURATION):
    validate("steady", [cpu], duration, MAX_CREATIONS)
    os.sched_setaffinity(0, {cpu})
    parent = os.getppid()
    identity = dict(pid=os.getpid(), parent_pid=parent, affinity=sorted(os.sched_getaffinity(0)),
                    membership=Path("/proc/self/cgroup").read_text())
    save(directory / "competitor_ready.json", identity)
    require_start(wait_file(directory / "workload_go.json", seconds=30))
    result = dict(identity, **work("steady", duration, MAX_CREATIONS, parent))
    result.update(final_affinity=sorted(os.sched_getaffinity(0)), final_membership=Path("/proc/self/cgroup").read_text())
    save(directory / "competitor_done.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--worker", choices=("steady", "churn"))
    mode.add_argument("--competitor", type=int)
    args = parser.parse_args()
    if args.worker:
        cpus = sorted(os.sched_getaffinity(0))
        if len(cpus) != 20:
            raise ValueError("Native production control requires exactly20 allowed CPUs")
        workload(args.directory.resolve(), args.worker, cpus)
    else:
        competitor(args.directory.resolve(), args.competitor)
