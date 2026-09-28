"""Threadripper scaling collector; no timing admission or scientific authorization.

Separate from frozen historical DGX collectors. Native affinity is 32 CPUs;
the scheduler reserves 64 SMT slots, with no hard 32-CPU quota.
"""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.measure_native_interval_step import run_command, step_memory
from benchmark_tools.measure_native_hierarchy_step import interval_point
from benchmark_tools.measure_native_lineage_step import evaluate as evaluate_lineage
from benchmark_tools.measure_native_root_context import read_point as read_root_point, evaluate, lineage_identity
from benchmark_tools.probe_host_counters import snapshot
from benchmark_tools.probe_dgx_step_separation import save, wait_file
from benchmark_tools.command_host_monitor import HostMonitor
from benchmark_tools.slurm_resource_snapshot import scoped_path

from benchmark_tools.observe_thread_affinity import observe
from benchmark_tools.probe_threadripper_allocation import inspect, validate as validate_allocation
from benchmark_tools.disk_observation_sequence import DiskObservations

TIMEOUT = 85800


def read_job_memory(scope):
    result = step_memory({"native_cpu_scope": str(scope)})
    directory = Path("/sys/fs/cgroup") / str(scope).lstrip("/")
    for name in ("memory.max", "memory.swap.max", "memory.stat"):
        try:
            result["raw"][name] = (directory / name).read_text()
        except OSError as error:
            result["errors"].append(dict(field=name, type=type(error).__name__, errno=error.errno))
    result["finished_ns"] = time.monotonic_ns()
    return result


def read_point(pid, membership, job_id, failure_path):
    point = read_root_point(pid, membership, job_id, failure_path)
    scope = Path("/sys/fs/cgroup") / str(scoped_path(membership, job_id)).lstrip("/")
    # Use the full user step subtree, not merely the anchor's task cgroup.
    while scope.name != "user":
        if scope.name.startswith("step_") or scope == scope.parent:
            raise ValueError("Require a Slurm user subtree")
        scope = scope.parent
    point["thread_affinity"] = observe(scope, range(32))
    return point


def validate(command, cpus, timeout_s, interval_s):
    if (not isinstance(command, list) or not command
            or not all(isinstance(v, str) and v for v in command)
            or not Path(command[0]).is_absolute()):
        raise ValueError("Require nonempty command with absolute executable")
    if (type(cpus) is not int or cpus != 32 or type(timeout_s) is not int or timeout_s != TIMEOUT
            or type(interval_s) not in (int, float) or interval_s != 1.):
        raise ValueError("Require frozen 32CPU/85800s/1s scaling settings")


def worker(directory):
    plan = json.loads((directory / "command.json").read_text())
    validate(plan["command"], plan["cpus"], plan["timeout_s"], plan["interval_s"])
    placement = inspect()
    validate_allocation(placement, placement, 64)
    if placement["affinity"] != list(range(32)):
        raise ValueError("Require the frozen CPU IDs 0-31")
    save(directory / "ready.json", dict(pid=os.getpid(), cgroup=Path("/proc/self/cgroup").read_text(),
                                       placement=placement))
    gate = wait_file(directory / "go.json")
    if not isinstance(gate, dict) or set(gate) != {"go"} or gate["go"] is not True:
        save(directory / "aborted_before_native.json", dict(status="observer_did_not_release_native"))
        return
    before = snapshot()
    start = time.monotonic_ns()
    with (directory / "native.log").open("x") as log:
        code, timed_out = run_command(plan["command"], log, TIMEOUT)
    finish = time.monotonic_ns()
    after = snapshot()
    save(directory / "done.json", dict(exit_code=code, timed_out=timed_out,
        started_ns=start, finished_ns=finish, snapshots=[before, after]))
    wait_file(directory / "release.json")


def measure(command, directory, job_id, cpus, memory_bytes, timeout_s, interval_s,
            monitor_host=True, host_interval_s=30.):
    validate(command, cpus, timeout_s, interval_s)
    if (type(memory_bytes) is not int or memory_bytes != 128*1024**3
            or monitor_host is not True or type(host_interval_s) not in (int, float) or host_interval_s != 30.):
        raise ValueError("Require fixed memory and observer settings")
    if (type(job_id) is not int or job_id <= 0 or os.environ.get("SLURM_JOB_ID") != str(job_id)
            or os.environ.get("SLURM_CPUS_PER_TASK") != "64"
            or os.environ.get("SLURM_MEM_PER_NODE") != "131072" or os.uname().nodename != "bizon"):
        raise ValueError("Require matching Threadripper64-slot/128GiB allocation")
    directory = Path(directory).absolute()
    directory.mkdir(exist_ok=False)
    save(directory / "command.json", dict(command=command, cpus=cpus, timeout_s=timeout_s, interval_s=interval_s))
    launched = ["srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1", "--cpus-per-task=64", "--cpu-bind=mask_cpu:0xffffffff",
                sys.executable, "-B", str(Path(__file__).resolve()), "--worker", str(directory)]
    with (directory / "step.log").open("x") as log, (directory / "host_processes.jsonl").open("x") as host_log:
        process = subprocess.Popen(launched, stdout=log, stderr=subprocess.STDOUT)
        try:
            ready = wait_file(directory / "ready.json")
            scope = scoped_path(ready["cgroup"], job_id)
            # Exclude all of this job, including its observer and sibling steps.
            job_scope = next(parent for parent in scope.parents if parent.name == f"job_{job_id}")
            host = HostMonitor(host_log, str(job_scope))
            host.observe()
            next_host = time.monotonic() + host_interval_s
            job_memory_before = read_job_memory(job_scope)
            save(directory / "job_memory_before.json", job_memory_before)
            # Inventory the host before starting the one-second point cadence.
            points = DiskObservations(directory)
            points.append(read_point(ready["pid"], ready["cgroup"], job_id, directory / "failed_point.json"))
            save(directory / "go.json", {"go": True})
            start = time.monotonic()
            index = 0
            while True:
                index += 1
                time.sleep(max(0, start + index * interval_s - time.monotonic()))
                completed = (directory / "done.json").exists()
                if process.poll() is not None:
                    raise RuntimeError("Worker exited before final observation")
                points.append(read_point(ready["pid"], ready["cgroup"], job_id, directory / "failed_point.json"))
                if completed or time.monotonic() >= next_host:
                    host.observe()
                    next_host = time.monotonic() + host_interval_s
                if completed:
                    break
                if time.monotonic() - start > TIMEOUT + 30:
                    raise TimeoutError("Native command exceeded timeout and cleanup allowance")
            done = json.loads((directory / "done.json").read_text())
            host_summary = host.summary(done["started_ns"] / 1e9, done["finished_ns"] / 1e9)
            save(directory / "host_process_summary.json", host_summary)
            memory = step_memory(interval_point(points[-1], job_id))
            save(directory / "step_memory.json", memory)
            job_memory_after = read_job_memory(job_scope)
            save(directory / "job_memory_after.json", job_memory_after)
            save(directory / "release.json", {"release": True})
            if process.wait(timeout=45) != 0:
                raise RuntimeError("Native step wrapper failed")
            measured = dict(status="command_exited_zero" if done["exit_code"] == 0 else "command_failed",
                schema="threadripper_scaling_v3",
                native=done, native_wall_s=(done["finished_ns"]-done["started_ns"])/1e9,
                job_id=job_id, launched=launched, placement=ready["placement"], point_records=points.records(), step_memory=memory,
                host_process_observation=host_summary,
                job_memory=dict(before=job_memory_before, after=job_memory_after),
                screening=evaluate_lineage(points, done, job_id), scientific_timings_admitted=False,
                controlled_workload_verified=False, publication_ready=False,
                limitations=["Threadripper collector, not controlled comparative timing admission.",
                    "Thread affinity is periodically observed, not a hard CPU quota or continuous guarantee.",
                    "Job peak is since cgroup creation through the post-native read, including preparation and observer; not process RSS.",
                    "Native-step and job peaks overlap and must not be added or baseline-subtracted.",
                    "All CPU flags and non-atomic read windows retained; no overhead subtraction.",
                    "Periodic process CPU observations miss short-lived work and do not establish quiet GPU/device-I/O activity.",
                    "Decoded raw history is disk-backed; final interval reports still scale with observation count after inference.",
                    "Raw observation writes and collection still consume resources; long-run overhead remains unvalidated."])
            save(directory / "lineage_report.json", measured)
            context = dict(status="native_root_context_measured", job_id=job_id,
                native_wall_s=measured["native_wall_s"], context=evaluate(points, job_id),
                lineage_report=lineage_identity(directory), scientific_timings_admitted=False,
                environmental_validity_established=False)
            save(directory / "root_context_report.json", context)
            return measured
        finally:
            # A failed initial observation must not release an unobserved native command.
            if not (directory / "go.json").exists():
                save(directory / "go.json", {"abort": True})
            if not (directory / "release.json").exists():
                save(directory / "release.json", {"release": True})
            if process.poll() is None:
                process.wait(timeout=TIMEOUT + 90)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=Path, required=True)
    worker(parser.parse_args().worker.resolve())
