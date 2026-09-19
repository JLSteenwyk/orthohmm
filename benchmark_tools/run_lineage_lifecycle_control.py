"""Finite owned-service lifecycle controls, not comparative tool timings."""

import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys
import time

from benchmark_tools.probe_cgroup_lineage import snapshot, compare, LineageSnapshotError
from benchmark_tools.probe_dgx_step_separation import burn, save


def membership():
    raw = Path("/proc/self/cgroup").read_text().strip()
    if not raw.startswith("0::/") or "\n" in raw:
        raise ValueError("Require unified cgroup membership")
    return raw[3:]


def assess(manager, target, load, unit, removed, ancestor_result, outside_result, cpu):
    path = PurePosixPath(load["cgroup"].strip().removeprefix("0::"))
    base = PurePosixPath(manager)
    if (not load["cgroup"].startswith("0::/") or base not in path.parents
            or path.name != unit or PurePosixPath(target) in path.parents
            or not .75 <= load["cpu_seconds"] < .85
            or load["affinity"] != [cpu] or not removed):
        raise ValueError("Unexpected service membership, affinity, CPU duration or lifecycle")
    retained = ancestor_result["aggregate_deltas"][-1]["cpu_usec"] / 1e6
    outside = outside_result["root_minus_target_cpu_usec"] / 1e6
    return dict(descendant_cpu_retention_response=retained >= .5,
                outside_target_response=outside >= .5,
                manager_cpu_s=retained, signed_outside_target_cpu_s=outside,
                scientific_timings_admitted=False)


def trial(directory, index, job, manager, target, cpu):
    directory.mkdir()
    root = Path("/sys/fs/cgroup")
    unit = f"orthohmm-lineage-{job}-{index}.service"
    # Runtime and affinity limits apply even though this owned service is outside Slurm.
    command = ["systemd-run", "--user", "--quiet", "--wait", "--pipe", "--collect",
        "--unit=" + unit, "--property=RuntimeMaxSec=10s", "--property=MemoryMax=256M",
        "--property=TasksMax=8", "taskset", "-c", str(cpu), sys.executable, "-B",
        "-m", "benchmark_tools.run_lineage_lifecycle_control", "--burn"]
    # systemd services do not inherit the submitter's working directory.
    command.insert(7, "--working-directory=" + str(Path.cwd()))
    save(directory / "command.json", dict(command=command, target=target, manager=manager, cpu=cpu))
    before_manager = snapshot(root, manager)
    before_target = snapshot(root, target)
    save(directory / "before.json", dict(manager=before_manager, target=before_target))
    completed = subprocess.run(command, text=True, capture_output=True, timeout=25)
    save(directory / "service.json", dict(returncode=completed.returncode,
         stdout=completed.stdout, stderr=completed.stderr))
    completed.check_returncode()
    load = json.loads(completed.stdout)
    service = root / load["cgroup"].strip().removeprefix("0::/")
    # Validate the returned path before using it for lifecycle observation.
    path = PurePosixPath(load["cgroup"].strip().removeprefix("0::"))
    if PurePosixPath(manager) not in path.parents or path.name != unit or ".." in path.parts:
        raise ValueError("Unexpected owned service scope")
    observations = []
    deadline = time.monotonic() + 15
    while True:
        exists = service.exists()
        observations.append(dict(monotonic_ns=time.monotonic_ns(), exists=exists))
        if not exists or time.monotonic() >= deadline:
            break
        time.sleep(.1)
    save(directory / "removal.json", dict(scope=str(path), observations=observations))
    after_manager = snapshot(root, manager)
    after_target = snapshot(root, target)
    save(directory / "after.json", dict(manager=after_manager, target=after_target))
    if membership() != target:
        raise ValueError("Observer migrated between snapshots")
    manager_result = compare(before_manager, after_manager)
    target_result = compare(before_target, after_target)
    result = assess(manager, target, load, unit, not exists, manager_result, target_result, cpu)
    return dict(status="lifecycle_control_evaluated", index=index, load=load,
                manager=manager_result, outside=target_result, result=result)


def run(output):
    if (os.uname().nodename != "spark-7ff0" or os.environ.get("SLURM_CPUS_PER_TASK") != "2"
            or len(os.sched_getaffinity(0)) != 2):
        raise ValueError("Require DGX two-CPU engineering allocation")
    job = int(os.environ["SLURM_JOB_ID"])
    target = membership()
    if f"job_{job}" not in PurePosixPath(target).parts:
        raise ValueError("Observer is not in this Slurm job")
    manager = subprocess.check_output(["systemctl", "--user", "show", "--property=ControlGroup",
                                       "--value"], text=True, timeout=10).strip()
    output.mkdir(exist_ok=False)
    sources = list(Path(__file__).parent.glob("*.py"))
    sources.append(Path(__file__).parent / "results/LINEAGE_LIFECYCLE_CONTROL_PROTOCOL_20260919.md")
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    save(output / "identity.json", dict(job=job, target=target, manager=manager,
        affinity=sorted(os.sched_getaffinity(0)), sources=hashes,
        uname=list(os.uname()), python=sys.version,
        executables={name: dict(path=str(Path(shutil.which(name)).resolve()),
            sha256=hashlib.sha256(Path(shutil.which(name)).read_bytes()).hexdigest())
            for name in (sys.executable, "systemd-run", "systemctl", "taskset")}))
    rows = []
    for index in range(3):
        directory = output / f"trial_{index}"
        try:
            row = trial(directory, index, job, manager, target, min(os.sched_getaffinity(0)))
        except Exception as error:
            row = dict(status="failed", index=index, error=str(error), error_type=type(error).__name__)
            if isinstance(error, LineageSnapshotError):
                row["partial_snapshot"] = error.evidence
        save(directory / "result.json", row)
        rows.append(row)
    unchanged = hashes == {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    report = dict(status="lifecycle_controls_finished", trials=rows, sources_unchanged=unchanged,
        all_controls_met=unchanged and all(r["status"] == "lifecycle_control_evaluated"
            and r["result"]["descendant_cpu_retention_response"]
            and r["result"]["outside_target_response"] for r in rows),
        scientific_timings_admitted=False, publication_ready=False)
    save(output / "report.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--burn", action="store_true")
    group.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.burn:
        print(json.dumps(dict(burn(), affinity=sorted(os.sched_getaffinity(0)))))
    else:
        result = run(args.output.absolute())
        raise SystemExit(0 if result["all_controls_met"] else 1)
