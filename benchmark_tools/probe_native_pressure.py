"""Read pressure from a separate native Slurm step, with outer host brackets."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_dgx_pressure import parse_pressure
from benchmark_tools.probe_dgx_step_separation import save, wait_file, scope_parts, validate_scopes
from benchmark_tools.probe_host_counters import snapshot, summarize

RESOURCES = ("cpu", "memory", "io")


def read_point(pid, membership, job, *, proc_root=Path("/proc"), group_root=Path("/sys/fs/cgroup"),
               host_reader=snapshot):
    member_path = proc_root / str(pid) / "cgroup"
    if member_path.read_text() != membership:
        raise ValueError("Native process membership changed")
    before = host_reader()
    native = dict(raw=dict(cgroup_membership=membership))
    validate_scopes(before, native, job)
    scope = str(Path(*scope_parts(native, job)))
    directory = group_root / scope.lstrip("/")
    for path in (directory, *directory.parents):
        if path == group_root.parent:
            break
        if path.is_symlink():
            raise ValueError("Symlink in native cgroup path")
    stat = directory.stat()
    identity = [stat.st_dev, stat.st_ino]
    fields = {}
    for resource in RESOURCES:
        start = time.monotonic_ns()
        raw = (directory / f"{resource}.pressure").read_text()
        end = time.monotonic_ns()
        fields[resource] = dict(raw=raw, started_ns=start, finished_ns=end,
                                totals=parse_pressure(raw, resource, system_scope=False))
    after_stat = directory.stat()
    if identity != [after_stat.st_dev, after_stat.st_ino] or member_path.read_text() != membership:
        raise ValueError("Native scope disappeared or identity changed")
    after = host_reader()
    result = dict(host=[before, after], native_membership=membership, scope=scope,
                  scope_identity=identity, pressure=fields, ticks=os.sysconf("SC_CLK_TCK"))
    validate(result, job)
    return result


def validate(point, job):
    before, after = point["host"]
    summarize(before, after, point["ticks"])
    native = dict(raw=dict(cgroup_membership=point["native_membership"]))
    for host in (before, after):
        validate_scopes(host, native, job)
    if point["scope"] != str(Path(*scope_parts(native, job))):
        raise ValueError("Pressure scope is not the native step")
    identity = point["scope_identity"]
    if len(identity) != 2 or any(type(v) is not int or v < 0 for v in identity):
        raise ValueError("Invalid cgroup identity")
    if set(point["pressure"]) != set(RESOURCES):
        raise ValueError("Incomplete native pressure fields")
    previous = before["finished_monotonic_ns"]
    for resource in RESOURCES:
        field = point["pressure"][resource]
        if (type(field["started_ns"]) is not int or type(field["finished_ns"]) is not int
                or not previous <= field["started_ns"] <= field["finished_ns"]):
            raise ValueError("Native pressure reads not enclosed or ordered")
        if field["totals"] != parse_pressure(field["raw"], resource, system_scope=False):
            raise ValueError("Native PSI totals disagree with raw evidence")
        previous = field["finished_ns"]
    if previous > after["started_monotonic_ns"]:
        raise ValueError("Native pressure exceeds host bracket")


def compare(left, right, job):
    for point in (left, right):
        validate(point, job)
    for key in ("native_membership", "scope", "scope_identity", "ticks"):
        if left[key] != right[key]:
            raise ValueError("Native pressure identity changed between observations")
    summarize(left["host"][0], right["host"][1], left["ticks"])
    if left["host"][1]["finished_monotonic_ns"] >= right["host"][0]["started_monotonic_ns"]:
        raise ValueError("Overlapping pressure observations")
    delta = {}
    for resource in RESOURCES:
        a, b = left["pressure"][resource], right["pressure"][resource]
        values = {k: b["totals"][k] - a["totals"][k] for k in ("some", "full")}
        if any(v < 0 for v in values.values()):
            raise ValueError("Native pressure counter decreased")
        delta[resource] = values
    return dict(status="native_step_pressure_diagnostic", native_stall_usec=delta,
                scientific_timings_admitted=False, controlled_workload_verified=False,
                limitations=["Native-step PSI includes all descendants and wrappers, not only one executable.",
                    "Host and native pressure cannot be subtracted to identify foreign interference.",
                    "Reads are non-atomic; no causal attribution, exclusion threshold or timing correction."])


def run(output):
    if os.environ.get("SLURM_CPUS_PER_TASK") != "2":
        raise ValueError("Require a scheduled two-CPU diagnostic allocation")
    job = int(os.environ["SLURM_JOB_ID"])
    output.mkdir(exist_ok=False)
    helper = Path(__file__).with_name("probe_dgx_step_separation.py")
    sources = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (
        Path(__file__), helper, helper.with_name("probe_host_counters.py"),
        helper.with_name("audit_dgx_pressure.py"), helper.with_name("summarize_frontier_overhead.py"))}
    worker_dir = output / "worker"
    worker_dir.mkdir()
    argv = ["srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1", "--cpus-per-task=1",
            sys.executable, "-B", str(helper), "--worker", str(worker_dir)]
    with (output / "step.log").open("x") as log:
        worker = subprocess.Popen(argv, stdout=log, stderr=subprocess.STDOUT)
        try:
            ready = wait_file(worker_dir / "ready.json")
            membership = wait_file(worker_dir / "native_before.json")["raw"]["cgroup_membership"]
            before = read_point(ready["pid"], membership, job)
            save(output / "before.json", before)
            time.sleep(1)
            after = read_point(ready["pid"], membership, job)
            save(output / "after.json", after)
            result = compare(before, after, job)
        finally:
            save(worker_dir / "release.json", {"release": True})
            code = worker.wait(timeout=60)
        if code != 0:
            raise ValueError("Native diagnostic worker failed")
    for name, digest in sources.items():
        if hashlib.sha256(helper.with_name(name).read_bytes()).hexdigest() != digest:
            raise ValueError("Diagnostic source changed")
    save(output / "result.json", dict(job_id=job, sources=sources, argv=argv, worker_exit_code=code,
         points=[before, after], result=result, publication_ready=False,
         limitations=["Sleeping-step interface check, not contention calibration or inference timing."]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args().output.resolve())
