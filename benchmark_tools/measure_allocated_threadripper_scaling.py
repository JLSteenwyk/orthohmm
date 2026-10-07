"""Versioned allocation-aware collector; no scientific or timing admission.

The old collector is hash-bound historical evidence and remains unchanged.
Keep its interval, completion, host, CPU and memory semantics, but obtain CPU
IDs from this running step instead of assuming that OS CPUs 0-31 are allocated.
"""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.measure_threadripper_scaling import TIMEOUT, validate, read_job_memory
from benchmark_tools.measure_native_interval_step import run_command, step_memory
from benchmark_tools.measure_native_hierarchy_step import interval_point
from benchmark_tools.measure_native_lineage_step import evaluate as evaluate_lineage
from benchmark_tools.measure_native_root_context import read_point as read_root_point, evaluate, lineage_identity
from benchmark_tools.probe_host_counters import snapshot
from benchmark_tools.probe_dgx_step_separation import save, wait_file
from benchmark_tools.command_host_monitor import HostMonitor
from benchmark_tools.periodic_host_observer import PeriodicHostObserver
from benchmark_tools.observe_threadripper_process_identity import enriched_snapshot
from benchmark_tools.slurm_resource_snapshot import scoped_path
from benchmark_tools.observe_thread_affinity import observe
from benchmark_tools.disk_observation_sequence import DiskObservations
from benchmark_tools.report_finalization import observe as observe_reporting
from benchmark_tools.native_completion import completion_evidence
from benchmark_tools.native_factorial_allocated_placement import bind, validate as validate_placement, sources as placement_sources
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.verify_lineage_native_provenance import same


SCHEMA = "allocated_threadripper_scaling_v1"
COMMAND_SCHEMA = "allocated_threadripper_command_v1"


def sources():
    return [record(__file__), *placement_sources(), record(Path(__file__).with_name("measure_threadripper_scaling.py"))]


def command_record(command):
    validate(command, 32, TIMEOUT, 1.)
    return dict(schema=COMMAND_SCHEMA, command=command, cpus=32, timeout_s=TIMEOUT,
                interval_s=1., sources=sources(), cpu_policy="actual_allocated_physical_cores_v1")


def validate_command(value):
    expected = command_record(value["command"])
    if not same(value, expected):
        raise ValueError("Versioned command, CPU policy, resources or sources differ")
    for item in value["sources"]:
        check(item)


def selected_ids(values):
    if (type(values) is not list or len(values) != 32
            or any(type(cpu) is not int or not 0 <= cpu < 192 for cpu in values)
            or values != sorted(set(values))):
        raise ValueError("Require 32 distinct selected OS CPUs")
    return values


def read_point(pid, membership, job_id, failure_path, allowed_cpus):
    allowed_cpus = selected_ids(allowed_cpus)
    point = read_root_point(pid, membership, job_id, failure_path)
    scope = Path("/sys/fs/cgroup") / str(scoped_path(membership, job_id)).lstrip("/")
    while scope.name != "user":
        if scope.name.startswith("step_") or scope == scope.parent:
            raise ValueError("Require a Slurm user subtree")
        scope = scope.parent
    point["thread_affinity"] = observe(scope, allowed_cpus)
    return point


def worker(directory):
    plan = json.loads((directory / "command.json").read_text())
    validate_command(plan)
    job = int(os.environ["SLURM_JOB_ID"])
    placement = bind(job)
    save(directory / "ready.json", dict(pid=os.getpid(),
        cgroup=placement["bound"]["cgroup"], placement=placement["bound"],
        allocated_placement=placement))
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
            monitor_host=True, host_interval_s=30., *, release_guard=None):
    validate(command, cpus, timeout_s, interval_s)
    if release_guard is not None and not callable(release_guard):
        raise ValueError("Release guard must be callable")
    if (type(memory_bytes) is not int or memory_bytes != 128*1024**3
            or monitor_host is not True or type(host_interval_s) not in (int, float) or host_interval_s != 30.):
        raise ValueError("Require fixed memory and observer settings")
    if (type(job_id) is not int or job_id <= 0 or os.environ.get("SLURM_JOB_ID") != str(job_id)
            or os.environ.get("SLURM_CPUS_PER_TASK") != "64"
            or os.environ.get("SLURM_MEM_PER_NODE") != "131072" or os.uname().nodename != "bizon"):
        raise ValueError("Require matching Threadripper64-slot/128GiB allocation")
    directory = Path(directory).absolute()
    directory.mkdir(exist_ok=False)
    save(directory / "command.json", command_record(command))
    launched = ["srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1", "--cpus-per-task=64",
                "--cpu-bind=cores", sys.executable, "-B", str(Path(__file__).resolve()), "--worker", str(directory)]
    with (directory / "step.log").open("x") as log, (directory / "host_processes.jsonl").open("x") as host_log:
        process = subprocess.Popen(launched, stdout=log, stderr=subprocess.STDOUT)
        periodic_host = None
        try:
            ready = wait_file(directory / "ready.json")
            allowed = validate_placement(ready["allocated_placement"], job_id)
            if (not same(ready["placement"], ready["allocated_placement"]["bound"])
                    or ready["pid"] != ready["placement"]["pid"]
                    or ready["cgroup"] != ready["placement"]["cgroup"]):
                raise ValueError("Ready worker differs from allocated placement")
            scope = scoped_path(ready["cgroup"], job_id)
            job_scope = next(parent for parent in scope.parents if parent.name == f"job_{job_id}")
            host = HostMonitor(host_log, str(job_scope), sample_fn=enriched_snapshot)
            host.observe()
            periodic_host = PeriodicHostObserver(host, host_interval_s)
            periodic_host.start(anchor=host.last_started)
            job_memory_before = read_job_memory(job_scope)
            save(directory / "job_memory_before.json", job_memory_before)
            points = DiskObservations(directory)
            if release_guard is not None:
                release_guard(directory)
                checked_at = time.monotonic()
            points.append(read_point(ready["pid"], ready["cgroup"], job_id, directory / "failed_point.json", allowed))
            if release_guard is not None and time.monotonic() - checked_at > 1:
                save(directory / "release_freshness_failed.json", dict(status="release_guard_stale"))
                raise ValueError("Release guard stale after initial observation")
            save(directory / "go.json", {"go": True})
            start = time.monotonic()
            index = 0
            while True:
                index += 1
                time.sleep(max(0, start + index * interval_s - time.monotonic()))
                completed = (directory / "done.json").exists()
                if process.poll() is not None:
                    raise RuntimeError("Worker exited before final observation")
                points.append(read_point(ready["pid"], ready["cgroup"], job_id, directory / "failed_point.json", allowed))
                if completed:
                    break
                if time.monotonic() - start > TIMEOUT + 30:
                    raise TimeoutError("Native command exceeded timeout and cleanup allowance")
            done = json.loads((directory / "done.json").read_text())
            completion = completion_evidence(points[0]["thread_affinity"], points[-1]["thread_affinity"],
                                             ready["pid"], done["finished_ns"])
            save(directory / "native_completion.json", completion)
            if completion["errors"]:
                raise ValueError("Native subtree completion unverified; retain failed attempt")
            host_summary = periodic_host.finish(done["started_ns"] / 1e9, done["finished_ns"] / 1e9)
            save(directory / "host_process_summary.json", host_summary)
            memory = step_memory(interval_point(points[-1], job_id))
            save(directory / "step_memory.json", memory)
            job_memory_after = read_job_memory(job_scope)
            save(directory / "job_memory_after.json", job_memory_after)
            save(directory / "release.json", {"release": True})
            if process.wait(timeout=45) != 0:
                raise RuntimeError("Native step wrapper failed")
            with observe_reporting(directory, job_id, job_scope, read_job_memory):
                measured = dict(status="command_exited_zero" if done["exit_code"] == 0 else "command_failed",
                    schema=SCHEMA, sources=sources(), native=done,
                    native_wall_s=(done["finished_ns"]-done["started_ns"])/1e9,
                    native_completion=completion, job_id=job_id, launched=launched,
                    placement=ready["placement"], allocated_placement=ready["allocated_placement"],
                    point_records=points.records(), step_memory=memory, host_process_observation=host_summary,
                    job_memory=dict(before=job_memory_before, after=job_memory_after),
                    screening=evaluate_lineage(points, done, job_id), scientific_timings_admitted=False,
                    controlled_workload_verified=False, publication_ready=False,
                    limitations=["Allocation-aware engineering collector, not scientific authorization or timing admission.",
                        "CPU IDs and NUMA placement depend on the actual scheduler allocation.",
                        "Thread affinity is periodic evidence, not a hard CPU quota or continuous guarantee.",
                        "Timings on a shared host may be affected by unknown, tool-dependent contention.",
                        "Native-step and job peaks overlap; no addition, baseline or overhead subtraction.",
                        "Report finalization separately extends job memory evidence through reporting.",
                        "Host observations miss short-lived work and do not establish host isolation.",
                        "Disk-backed raw observations retain all flags; long-run observer overhead remains unvalidated."])
                save(directory / "lineage_report.json", measured)
                save(directory / "root_context_report.json", dict(status="native_root_context_measured", job_id=job_id,
                    native_wall_s=measured["native_wall_s"], context=evaluate(points, job_id),
                    lineage_report=lineage_identity(directory), scientific_timings_admitted=False,
                    environmental_validity_established=False))
            return measured
        finally:
            try:
                if periodic_host is not None:
                    periodic_host.close()
            finally:
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
