"""Bounded cgroup accounting control for preparation-owned tmpfs input pages."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.slurm_resource_snapshot import scoped_path, counters
from benchmark_tools.snapshot_orthohmm_input_order import record

SIZE = 64 * 1024**2


def scope(job, kind):
    membership = Path("/proc/self/cgroup").read_text()
    path = scoped_path(membership, job)
    name = f"job_{job}" if kind == "job" else next(p.name for p in path.parents if p.name.startswith("step_"))
    path = next(p for p in (path, *path.parents) if p.name == name)
    return Path("/sys/fs/cgroup") / str(path).lstrip("/")


def snapshot(path):
    start = time.monotonic_ns()
    raw = {name: (path / name).read_text() for name in (
        "memory.current", "memory.peak", "memory.max", "memory.swap.max", "memory.stat", "memory.events")}
    return dict(scope=str(path), started_ns=start, finished_ns=time.monotonic_ns(), raw=raw)


def read_payload(path, output):
    job = int(os.environ["SLURM_JOB_ID"])
    target = scope(job, "step")
    before = snapshot(target)
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024**2), b""):
            digest.update(block)
            size += len(block)
    after = snapshot(target)
    save(output, dict(before=before, after=after, bytes_read=size, sha256=digest.hexdigest()))


def evaluate(before, prepared, native, expected_sha):
    name = Path(before["scope"]).name
    if before["scope"] != prepared["scope"] or not name.startswith("job_") or not name[4:].isdigit():
        raise ValueError("Changed job scope")
    parent = Path(before["scope"])
    child = Path(native["before"]["scope"])
    if (child.parent != parent or not child.name.startswith("step_")
            or child.name in ("step_batch", "step_extern") or native["after"]["scope"] != str(child)):
        raise ValueError("Require distinct native step within the same job")
    if native["bytes_read"] != SIZE or native["sha256"] != expected_sha:
        raise ValueError("Native read did not reproduce prepared bytes")
    observations = (before, prepared, native["before"], native["after"])
    for row in observations:
        raw = row["raw"]
        if int(raw["memory.current"]) > int(raw["memory.peak"]):
            raise ValueError("Invalid memory gauges")
    for a, b in zip(observations, observations[1:]):
        if a["finished_ns"] > b["started_ns"]:
            raise ValueError("Memory observation windows overlap")
    job_delta = counters(prepared["raw"]["memory.stat"])["shmem"] - counters(before["raw"]["memory.stat"])["shmem"]
    native_shmem = counters(native["after"]["raw"]["memory.stat"])["shmem"]
    if job_delta < SIZE or native_shmem >= SIZE // 2:
        raise ValueError("Expected preparation-owned charge was not observed")
    return dict(prepared_job_shmem_delta_bytes=job_delta, native_step_shmem_after_bytes=native_shmem,
                status="preparation_owned_tmpfs_charge_observed", scientific_timings_admitted=False)


def run(directory):
    directory = directory.absolute()
    directory.mkdir(exist_ok=False)
    job = int(os.environ["SLURM_JOB_ID"])
    payload = Path(f"/dev/shm/orthohmm_memory_charge_{job}.bin")
    if payload.exists() or payload.is_symlink():
        raise FileExistsError(payload)
    before = snapshot(scope(job, "job"))
    save(directory / "before.json", before)
    with payload.open("xb") as handle:
        digest = hashlib.sha256()
        block = b"X" * 1024**2
        for _ in range(64):
            handle.write(block)
            digest.update(block)
    try:
        time.sleep(1)
        prepared = snapshot(scope(job, "job"))
        save(directory / "prepared.json", prepared)
        command = ["srun", "--exclusive", "--exact", "-N1", "-n1", "-c1",
                   sys.executable, "-B", "-m", "benchmark_tools.probe_tmpfs_memory_charge",
                   "--read", str(payload), "--output", str(directory / "native.json")]
        save(directory / "command.json", dict(argv=command, payload_bytes=SIZE, sha256=digest.hexdigest()))
        with (directory / "native.log").open("x") as log:
            completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=60)
        save(directory / "exit.json", dict(exit_code=completed.returncode))
        if completed.returncode:
            raise ValueError("Native memory control failed")
        native = json.loads((directory / "native.json").read_text())
        result = evaluate(before, prepared, native, digest.hexdigest())
        save(directory / "result.json", dict(result, job_id=job, source=record(__file__),
            limitations=["Single 64-MiB read control, not full pipeline memory use.",
                         "Job memory includes preparation and observer; native-step peak alone excludes these input pages.",
                         "Counters are non-atomic observations; do not sum scope peaks or subtract this control from timings."]))
    finally:
        payload.unlink()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path)
    parser.add_argument("--read", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.read:
        read_payload(args.read, args.output)
    else:
        run(args.directory)
