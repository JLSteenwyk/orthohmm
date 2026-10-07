"""Bind only this worker to 32 physical cores in its actual Slurm step.

New placement evidence is intentionally separate from the frozen 0-31 route.
It neither changes scheduler policy nor authorizes scientific execution.
"""

import os
from pathlib import Path
import time

from benchmark_tools.native_factorial_cpu_selection import choose, host_topology
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.probe_threadripper_allocation import inspect, validate as validate_allocation
from benchmark_tools.slurm_resource_snapshot import cpus, scoped_path
from benchmark_tools.verify_lineage_native_provenance import same


SCHEMA = "native_factorial_allocated_placement_v1"


def sources():
    root = Path(__file__).resolve().parent
    return [record(root / name) for name in (
        "native_factorial_allocated_placement.py", "native_factorial_cpu_selection.py",
        "probe_threadripper_allocation.py", "slurm_resource_snapshot.py")]


def validate(value, job_id):
    if type(job_id) is not int or job_id <= 0 or value.get("schema") != SCHEMA:
        raise ValueError("Require versioned allocated placement and positive job")
    before, after, topology = value["allocated"], value["bound"], value["allocated_topology"]
    selection = choose(before["affinity"], topology)
    if not same(selection, value["selection"]):
        raise ValueError("Allocated-core selection does not reproduce")
    if (any(type(value[k]) is not int or value[k] <= 0 for k in ("started_ns", "finished_ns"))
            or value["started_ns"] > value["finished_ns"]
            or type(before["pid"]) is not int or before["pid"] <= 0
            or before["pid"] != after["pid"]
            or value["affinity_changed"] is not True
            or any(value[k] is not False for k in (
                "scientific_execution_authorized", "scientific_timings_admitted", "publication_ready"))):
        raise ValueError("Invalid placement window, worker identity or admission")
    for key in ("host", "pid", "cgroup", "ancestors", "slurm"):
        if not same(before[key], after[key]):
            raise ValueError("Allocation context changed during own-worker binding")
    if (after["affinity"] != selection["native_cpu_ids"]
            or before["slurm"]["SLURM_JOB_ID"] != str(job_id)):
        raise ValueError("Bound CPUs or job differ from actual allocation")
    expected_before = sorted(
        [{k: row[k] for k in ("cpu", "package", "core")} for row in topology],
        key=lambda row: row["cpu"])
    expected_after = [row for row in expected_before if row["cpu"] in after["affinity"]]
    if not same(before["topology"], expected_before) or not same(after["topology"], expected_after):
        raise ValueError("Probe and selection topology disagree")
    projected = dict(before, affinity=after["affinity"], topology=expected_after)
    validate_allocation(projected, after, 64)
    scope = scoped_path(before["cgroup"], job_id)
    step = scope.parts[scope.parts.index(f"job_{job_id}") + 1]
    if step in ("step_batch", "step_extern"):
        raise ValueError("Binding must occur in a separate native step")
    root = Path("/sys/fs/cgroup")
    expected_paths = []
    current = root / str(scope).lstrip("/")
    while current != root:
        expected_paths.append(str(current))
        current = current.parent
        if not current.is_relative_to(root):
            raise ValueError("Invalid cgroup ancestry")
    if [row["path"] for row in before["ancestors"]] != expected_paths:
        raise ValueError("Missing, duplicated or reordered cgroup ancestors")
    for row in before["ancestors"]:
        if not set(before["affinity"]) <= set(cpus(row["cpuset.cpus.effective"])):
            raise ValueError("Allocated CPUs outside effective ancestor cpuset")
    expected_sources = sources()
    if not same(value["sources"], expected_sources):
        raise ValueError("Placement source bindings differ")
    for item in value["sources"]:
        check(item)
    return selection["native_cpu_ids"]


def bind(job_id):
    started = time.monotonic_ns()
    before = inspect()
    topology = host_topology(before["affinity"])
    selection = choose(before["affinity"], topology)
    # PID zero means this process only; no other job/process is re-affinitized.
    os.sched_setaffinity(0, selection["native_cpu_ids"])
    after = inspect()
    value = dict(schema=SCHEMA, allocated=before, bound=after,
        allocated_topology=topology, selection=selection, affinity_changed=True,
        started_ns=started, finished_ns=time.monotonic_ns(), sources=sources(),
        scientific_execution_authorized=False, scientific_timings_admitted=False,
        publication_ready=False)
    validate(value, job_id)
    return value
