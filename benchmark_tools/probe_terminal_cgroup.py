"""One owned Slurm diagnostic of read-only cgroup descriptors across job exit."""

import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess
import time


def worker(directory):
    marker = directory / "worker.json"
    marker.write_text(json.dumps(dict(pid=os.getpid(), job=os.environ["SLURM_JOB_ID"],
                                    membership=Path("/proc/self/cgroup").read_text())) + "\n")
    time.sleep(10)
    data = bytearray(32 * 1024 * 1024)
    for i in range(0, len(data), 4096):
        data[i] = 1
    time.sleep(3)
    del data
    time.sleep(2)


def job_scope(membership, job):
    paths = [line[3:] for line in membership.splitlines() if line.startswith("0::")]
    if len(paths) != 1:
        raise ValueError("Require one unified cgroup")
    path = Path(paths[0])
    if not path.is_absolute() or ".." in path.parts or path.parts.count(f"job_{job}") != 1:
        raise ValueError("Wrong job scope")
    scope = next(p for p in path.parents if p.name == f"job_{job}")
    return Path("/sys/fs/cgroup") / str(scope).lstrip("/")


def probe(directory):
    directory.mkdir(parents=True, exist_ok=False)
    command = ["sbatch", "--parsable", "--partition=all", "--nodelist=bizon", "--nodes=1",
               "--ntasks=1", "--cpus-per-task=1", "--mem=256M", "--time=00:02:00",
               "--job-name=orthohmm-final-counter-probe", f"--output={directory}/slurm.log",
               "--wrap=" + shlex.join(["/usr/bin/python3", "-I", "-B", str(Path(__file__).resolve()),
                                        "--worker", "--output", str(directory)])]
    submission = subprocess.run(command, capture_output=True, text=True, check=True, timeout=30)
    job = int(submission.stdout.strip().split(";")[0])
    (directory / "submission.json").write_text(json.dumps(dict(command=command, job_id=job,
        stdout=submission.stdout, stderr=submission.stderr), indent=2) + "\n")
    print(json.dumps(dict(job_id=job, status="submitted_diagnostic_only")), flush=True)
    deadline = time.monotonic() + 180
    marker = directory / "worker.json"
    while not marker.exists():
        if time.monotonic() > deadline:
            raise TimeoutError("Marker unavailable; inspect existing job, do not restart")
        time.sleep(.1)
    info = json.loads(marker.read_text())
    if info["job"] != str(job):
        raise ValueError("Worker job mismatch")
    scope = job_scope(info["membership"], job)
    names = ("cgroup.events", "memory.peak", "memory.current", "cpu.stat")
    descriptors = {}
    rows = []
    try:
        for name in names:
            descriptors[name] = os.open(scope / name, os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW)
        identity = scope.stat()
        while time.monotonic() < deadline:
            row = dict(started_ns=time.monotonic_ns(), raw={}, errors={})
            for name, fd in descriptors.items():
                try:
                    row["raw"][name] = os.pread(fd, 16384, 0).decode()
                except OSError as error:
                    row["errors"][name] = dict(errno=error.errno, message=str(error))
            try:
                now = scope.stat()
                row["same_scope_identity"] = (now.st_dev, now.st_ino) == (identity.st_dev, identity.st_ino)
            except OSError as error:
                row["scope_stat_error"] = dict(errno=error.errno, message=str(error))
            row["finished_ns"] = time.monotonic_ns()
            rows.append(row)
            if row["errors"] or "scope_stat_error" in row:
                break
            time.sleep(.01)
    finally:
        for fd in descriptors.values():
            os.close(fd)
    (directory / "observations.json").write_text(json.dumps(dict(job_id=job, scope=str(scope), rows=rows), indent=2) + "\n")
    # Missing queue membership is checked against accounting, not inferred terminal.
    for _ in range(60):
        terminal = subprocess.run(["sacct", "-j", str(job), "--parsable2", "--noheader",
                                   "--format=JobIDRaw,State,ExitCode"], capture_output=True, text=True, check=True)
        allocation = [line.split("|") for line in terminal.stdout.splitlines() if line.split("|")[0] == str(job)]
        if allocation and allocation[0][1] in {"COMPLETED", "FAILED", "TIMEOUT", "CANCELLED", "OUT_OF_MEMORY", "NODE_FAIL"}:
            break
        time.sleep(1)
    result = dict(job_id=job, scheduler_stdout=terminal.stdout, scheduler_stderr=terminal.stderr,
        observed_empty_rows=sum("populated 0" in r["raw"].get("cgroup.events", "") and not r["errors"] for r in rows),
        observations=len(rows), last_observation=rows[-1] if rows else None,
        scientific_timings_admitted=False, complete_job_accounting_verified=False,
        limitations=["One small diagnostic; reads are non-atomic and a removed cgroup may make descriptors unreadable.",
                     "An empty observation alone does not establish final counter stability or production readiness."])
    (directory / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    (worker if args.worker else probe)(args.output.resolve())
