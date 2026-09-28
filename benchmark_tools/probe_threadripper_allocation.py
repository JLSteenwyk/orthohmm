"""Read-only placement/limit probe; never authorizes benchmark timing."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def inspect():
    affinity = sorted(os.sched_getaffinity(0))
    topology = []
    for cpu in affinity:
        root = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        topology.append(dict(cpu=cpu, package=int((root / "physical_package_id").read_text()),
                             core=int((root / "core_id").read_text())))
    raw = Path("/proc/self/cgroup").read_text()
    rows = [line[3:] for line in raw.splitlines() if line.startswith("0::")]
    if len(rows) != 1:
        raise ValueError("Require unified cgroup identity")
    relative = Path(rows[0])
    if not relative.is_absolute() or ".." in relative.parts:
        raise ValueError("Invalid cgroup path")
    root = Path("/sys/fs/cgroup")
    current = root / str(relative).lstrip("/")
    limits = []
    while current != root:
        limits.append(dict(path=str(current), **{
            name: (current / name).read_text().strip()
            for name in ("cpuset.cpus.effective", "memory.max", "memory.swap.max", "cpu.max")
        }))
        current = current.parent
    return dict(host=os.uname().nodename, pid=os.getpid(), affinity=affinity,
                topology=topology, cgroup=raw, ancestors=limits,
                slurm={k: os.environ.get(k) for k in
                       ("SLURM_JOB_ID", "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE")})


def validate(parent, child):
    for sample in (parent, child):
        if (sample["host"] != "bizon" or len(sample["affinity"]) != 32
                or len(set(sample["affinity"])) != 32
                or sorted(t["cpu"] for t in sample["topology"]) != sample["affinity"]
                or len({(t["package"], t["core"]) for t in sample["topology"]}) != 32):
            raise ValueError("Require 32 distinct physical-core placements on bizon")
        if (sample["slurm"]["SLURM_CPUS_PER_TASK"] != "32"
                or sample["slurm"]["SLURM_MEM_PER_NODE"] != "131072"):
            raise ValueError("Requested allocation differs")
        caps = [int(row["memory.max"]) for row in sample["ancestors"]
                if row["memory.max"] != "max"]
        if not caps or min(caps) != 128 * 1024**3:
            raise ValueError("Require effective 128 GiB RAM cap")
    if (parent["affinity"] != child["affinity"] or parent["cgroup"] != child["cgroup"]
            or parent["slurm"] != child["slurm"]):
        raise ValueError("Child did not inherit allocation and placement")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.snapshot:
        print(json.dumps(inspect()))
        return
    if args.output is None or args.output.exists():
        raise ValueError("Require a fresh output path")
    parent = inspect()
    time.sleep(2)
    child = json.loads(subprocess.check_output(
        [sys.executable, str(Path(__file__).resolve()), "--snapshot"], text=True, timeout=30))
    report = dict(parent=parent, child=child, scientific_timings_admitted=False,
                  limitations=["No native inference or collector tested.",
                               "CPU quota, swap and all raw ancestor limits are recorded separately.",
                               "No quiet-host or full-run placement guarantee."])
    try:
        validate(parent, child)
        report["status"] = "placement_and_ram_cap_verified"
    except ValueError as error:
        report.update(status="allocation_probe_failed", error=str(error))
    with args.output.open("x") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    if report["status"] != "placement_and_ram_cap_verified":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
