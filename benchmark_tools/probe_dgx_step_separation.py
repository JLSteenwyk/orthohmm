"""Engineering test of completed sibling CPU work outside a sleeping Slurm step."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent))
from probe_host_counters import snapshot, summarize, parse_group


def save(path, value):
    pending = path.with_name(path.name + ".pending")
    with pending.open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
    os.link(pending, path)
    pending.unlink()


def wait_file(path, seconds=45):
    deadline = time.monotonic() + seconds
    while not path.exists():
        if time.monotonic() >= deadline:
            raise TimeoutError(str(path))
        time.sleep(.02)
    return json.loads(path.read_text())


def worker(directory):
    first = snapshot()
    # Publish ready only after closing the complete snapshot file.
    save(directory / "native_before.json", first)
    save(directory / "ready.json", {"pid": os.getpid()})
    wait_file(directory / "release.json")
    save(directory / "native_after.json", snapshot())


def burn():
    start = time.process_time()
    value = 1
    while time.process_time() - start < .75:
        for _ in range(1000):
            value = (value * 1664525 + 1013904223) % 4294967296
    return dict(pid=os.getpid(), cgroup=Path("/proc/self/cgroup").read_text(),
                cpu_seconds=time.process_time() - start, checksum=value)


def scope_parts(sample, job):
    path = Path(parse_group(sample["raw"]["cgroup_membership"]))
    if path.parts.count(f"job_{job}") != 1:
        raise ValueError("Wrong Slurm job scope")
    index = path.parts.index(f"job_{job}")
    if len(path.parts) <= index + 1 or not path.parts[index + 1].startswith("step_"):
        raise ValueError("Missing step scope")
    return path.parts[:index + 2]


def validate_scopes(parent, native, job):
    batch, child = scope_parts(parent, job), scope_parts(native, job)
    if batch[-1] != "step_batch" or child[-1] in ("step_batch", "step_extern") or batch[:-1] != child[:-1]:
        raise ValueError("Native and observer are not separate steps of the same job")
    return {"observer_step": "/".join(batch), "native_step": "/".join(child)}


def run(output):
    job = int(os.environ["SLURM_JOB_ID"])
    if os.uname().nodename != "spark-7ff0" or os.environ.get("SLURM_CPUS_PER_TASK") != "2":
        raise ValueError("Require scheduled two-CPU engineering task on spark-7ff0")
    output.mkdir(exist_ok=False)
    source = Path(__file__).resolve()
    sources = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
               for p in (source, source.with_name("probe_host_counters.py"))}
    trials = []
    for mode in ("quiet", "burst"):
        directory = output / mode
        directory.mkdir()
        command = ["srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1", "--cpus-per-task=1",
                   sys.executable, "-B", str(source), "--worker", str(directory)]
        with (directory / "native.log").open("x") as log:
            process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
            try:
                ready = wait_file(directory / "ready.json")
                native_before = json.loads((directory / "native_before.json").read_text())
                before = snapshot()
                scopes = validate_scopes(before, native_before, job)
                child = None
                if mode == "burst":
                    finished = subprocess.run([sys.executable, "-B", str(source), "--burn"],
                        check=True, capture_output=True, text=True, timeout=15)
                    child = json.loads(finished.stdout)
                    if parse_group(child["cgroup"]) != parse_group(before["raw"]["cgroup_membership"]):
                        raise ValueError("Burst escaped batch scope")
                    if not .75 <= child["cpu_seconds"] < 5:
                        raise ValueError("Unexpected burst CPU consumption")
                else:
                    time.sleep(1)
                after = snapshot()
                save(directory / "release.json", {"release": True})
                if process.wait(timeout=45) != 0:
                    raise RuntimeError("Native step failed")
                native_after = json.loads((directory / "native_after.json").read_text())
                validate_scopes(after, native_after, job)
                summary = summarize(before, after, os.sysconf("SC_CLK_TCK"))
                if mode == "burst" and summary["accounted_host_busy_cpu_s"] < .5:
                    raise ValueError("Host counters did not retain the known completed CPU burst")
                trials.append(dict(mode=mode, command=command, ready=ready, scopes=scopes, burst=child,
                    snapshots=[before, after], native_snapshots=[native_before, native_after], summary=summary,
                    native_exit_code=process.returncode))
            finally:
                if process.poll() is None:
                    # Release our own waiting worker even when observation fails.
                    if not (directory / "release.json").exists():
                        save(directory / "release.json", {"release": True, "failed_observation": True})
                    process.wait(timeout=60)
    if sources != {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                   for p in (source, source.with_name("probe_host_counters.py"))}:
        raise ValueError("Probe sources changed")
    report = dict(status="step_separation_engineering_complete", job_id=job, sources=sources, trials=trials,
        publication_ready=False, controlled_workload_verified=False,
        limitations=["One quiet and one known CPU-burst trial, not general calibration or a benchmark.",
            "Burst is outside native step but inside the same allocation, not an unrelated user's job.",
            "Host totals include observer and other work; no foreign-load estimate or bound is inferred.",
            "Native worker sleeps; compute/process-heavy overhead and memory-peak behavior remain untested."])
    save(output / "report.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--worker", type=Path)
    mode.add_argument("--burn", action="store_true")
    mode.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.worker:
        worker(args.worker)
    elif args.burn:
        print(json.dumps(burn()))
    else:
        run(args.output.absolute())
