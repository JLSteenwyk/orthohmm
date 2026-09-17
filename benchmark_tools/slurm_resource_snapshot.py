"""Read scoped cgroup-v2 and process RSS evidence without changing job limits."""

import argparse
from datetime import datetime, timezone
import json
import hashlib
from pathlib import Path, PurePosixPath
import re
import time

import psutil


def scoped_path(text, job_id):
    rows = [line.split(":", 2) for line in text.splitlines() if line.startswith("0::")]
    if len(rows) != 1 or job_id <= 0:
        raise ValueError("Require one unified cgroup and a positive Slurm job ID")
    path = PurePosixPath(rows[0][2])
    if not path.is_absolute() or ".." in path.parts or path.parts.count(f"job_{job_id}") != 1:
        raise ValueError("Process is not within the requested Slurm job")
    index = path.parts.index(f"job_{job_id}")
    if len(path.parts) <= index + 1 or not path.parts[index + 1].startswith("step_"):
        raise ValueError("Require an identifiable Slurm step subtree")
    return path


def counters(text):
    values = {}
    for line in text.splitlines():
        fields = line.split()
        if len(fields) != 2 or fields[0] in values or not fields[1].isdigit():
            raise ValueError("Malformed or duplicate cgroup counter")
        values[fields[0]] = int(fields[1])
    return values


def cpus(text):
    result = set()
    for field in text.strip().split(","):
        if not re.fullmatch(r"\d+(?:-\d+)?", field):
            raise ValueError("Malformed cpuset")
        bounds = [int(value) for value in field.split("-")]
        start, end = bounds[0], bounds[-1]
        if end < start or end - start > 100000:
            raise ValueError("Invalid cpuset range")
        result.update(range(start, end + 1))
    return sorted(result)


def metrics(directory):
    names = ("cpu.stat", "memory.current", "memory.peak", "memory.stat", "memory.events", "cpuset.cpus.effective", "memory.max")
    raw = {name: (directory / name).read_text() for name in names}
    cpu = counters(raw["cpu.stat"])
    if not {"usage_usec", "user_usec", "system_usec"} <= cpu.keys():
        raise ValueError("Missing cumulative CPU accounting")
    for name in ("memory.current", "memory.peak"):
        if not raw[name].strip().isdigit():
            raise ValueError("Missing numeric cgroup memory accounting")
    limit = raw["memory.max"].strip()
    if limit != "max" and not limit.isdigit():
        raise ValueError("Invalid memory limit")
    return {"raw": raw, "cpu": cpu, "memory_current_bytes": int(raw["memory.current"]),
            "memory_peak_since_creation_or_reset_bytes": int(raw["memory.peak"]),
            "memory_stat": counters(raw["memory.stat"]), "memory_events": counters(raw["memory.events"]),
            "effective_cpus": cpus(raw["cpuset.cpus.effective"]), "memory_max_bytes": None if limit == "max" else int(limit)}


def snapshot(pid, job_id):
    started = time.monotonic()
    cgroup_text = Path(f"/proc/{pid}/cgroup").read_text()
    scope = scoped_path(cgroup_text, job_id)
    directory = Path("/sys/fs/cgroup") / str(scope).lstrip("/")
    result = {"status": "resource_snapshot", "timestamp_utc": datetime.now(timezone.utc).isoformat(),
              "source": {"path": str(Path(__file__).resolve()), "sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
              "psutil_version": psutil.__version__,
              "pid": pid, "job_id": job_id, "proc_cgroup": cgroup_text, "scope": str(scope),
              "metrics": metrics(directory), "processes": [], "sampling_errors": [], "ancestor_limits": []}
    current = directory
    while True:
        value = (current / "memory.max").read_text().strip()
        if value != "max" and not value.isdigit():
            raise ValueError("Invalid ancestor memory limit")
        result["ancestor_limits"].append({"path": str(current), "memory_max": value})
        if current.name == f"job_{job_id}":
            break
        current = current.parent
    identifiers = set()
    for path in [directory / "cgroup.procs", *directory.glob("**/cgroup.procs")]:
        identifiers.update(int(value) for value in path.read_text().split())
    for identifier in sorted(identifiers):
        try:
            process = psutil.Process(identifier)
            created = process.create_time()
            rss = process.memory_info().rss
            membership = scoped_path(Path(f"/proc/{identifier}/cgroup").read_text(), job_id)
            if not membership.is_relative_to(scope) or not process.is_running():
                raise ValueError("Process left the sampled subtree or exited")
            result["processes"].append({"pid": identifier, "created": created, "rss_bytes": rss,
                                        "name": process.name(), "cgroup": str(membership)})
        except (psutil.Error, OSError, ValueError) as error:
            result["sampling_errors"].append({"pid": identifier, "error_type": type(error).__name__, "error": str(error)})
    if scoped_path(Path(f"/proc/{pid}/cgroup").read_text(), job_id) != scope:
        raise ValueError("Anchor process changed cgroup during collection")
    result.update(sampled_sum_process_rss_bytes=sum(p["rss_bytes"] for p in result["processes"]),
                  collection_seconds=time.monotonic() - started,
                  limitations=["Read-only snapshot, not a completed timing run or historical peak-RSS reconstruction.",
                    "Cgroup memory includes charged file cache and kernel memory; it is not process RSS.",
                    "Cgroup peak is since creation or its last reset; this collector never resets it. CPU values are cumulative microseconds.",
                    "Process RSS is a non-atomic sample and shared pages may be counted more than once; short-lived processes may be missed.",
                    "Sampling errors make the RSS sum incomplete. This does not verify host quietness or exclusivity.",
                    "Scope is the anchor process subtree, not necessarily the entire job. Ancestor limits shown only through the job root."])
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pid", type=int, required=True)
    parser.add_argument("--job-id", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = snapshot(args.pid, args.job_id)
    with args.output.open("x") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
