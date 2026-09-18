"""Fixed-work process-churn engineering probe, not scientific timing admission."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent))
from probe_host_counters import snapshot, summarize, parse_group
from probe_dgx_step_separation import save, wait_file, validate_scopes

BUFFER_BYTES = 16 * 1024 * 1024
ROUNDS = 256
CHILDREN = 4
MODES = (False, True, True, False, False, True)
INTERVAL = .2


def checksum(rounds=ROUNDS, size=BUFFER_BYTES):
    payload = b"x" * size
    result = hashlib.sha256()
    for _ in range(rounds):
        result.update(hashlib.sha256(payload).digest())
    return result.hexdigest()


def expected_checksum():
    return hashlib.sha256(hashlib.sha256(b"x" * BUFFER_BYTES).digest() * ROUNDS).hexdigest()


def load():
    start = time.process_time()
    value = checksum()
    return dict(checksum=value, cpu_s=time.process_time() - start,
                cgroup=Path("/proc/self/cgroup").read_text(), pid=os.getpid())


def cpu_usage(sample):
    rows = dict(row.split() for row in sample["optional"]["cgroup_cpu.stat"].splitlines())
    return int(rows["usage_usec"]) / 1e6


def validate_trial(trial, job):
    before, after = trial["native_snapshots"]
    if before["errors"] or after["errors"] or any(s["errors"] for s in trial["snapshots"]):
        raise ValueError("Counter read failure")
    summary = summarize(before, after, trial["clock_ticks_per_second"])
    for sample in trial["snapshots"]:
        validate_scopes(sample, before, job)
    start, end = trial["work_started_ns"], trial["work_finished_ns"]
    if not before["finished_monotonic_ns"] < start < end < after["started_monotonic_ns"]:
        raise ValueError("Native snapshots do not bracket work")
    if len(trial["children"]) != CHILDREN or trial["native_exit_code"] != 0:
        raise ValueError("Incomplete fixed work")
    if trial["monitor"] and len(trial["snapshots"]) < 2:
        raise ValueError("Insufficient monitored snapshots")
    if not trial["monitor"] and len(trial["snapshots"]) != 1:
        raise ValueError("Unexpected unmonitored snapshots")
    expected = expected_checksum()
    for child in trial["children"]:
        if child["checksum"] != expected or child["cpu_s"] <= 0:
            raise ValueError("Invalid fixed-work result")
        if parse_group(child["cgroup"]) != parse_group(before["raw"]["cgroup_membership"]):
            raise ValueError("Child escaped native scope")
    delta = cpu_usage(after) - cpu_usage(before)
    child_cpu = sum(c["cpu_s"] for c in trial["children"])
    if delta < child_cpu * .95:
        raise ValueError("Native cgroup did not retain child CPU")
    peak = int(after["optional"]["cgroup_memory.peak"])
    if peak < BUFFER_BYTES:
        raise ValueError("Native peak below known buffer allocation")
    return dict(work_wall_s=(end - start) / 1e9, child_cpu_s=child_cpu,
                native_cgroup_cpu_s=delta, native_memory_peak_bytes=peak,
                host_summary=summary)


def worker(directory):
    save(directory / "native_before.json", snapshot())
    save(directory / "ready.json", {"pid": os.getpid()})
    wait_file(directory / "go.json")
    started = time.monotonic_ns()
    children = []
    for _ in range(CHILDREN):
        completed = subprocess.run([sys.executable, "-I", "-B", str(Path(__file__).resolve()), "--load"],
                                   check=True, capture_output=True, text=True, timeout=90)
        children.append(json.loads(completed.stdout))
    finished = time.monotonic_ns()
    save(directory / "native_after.json", snapshot())
    save(directory / "done.json", dict(children=children, work_started_ns=started, work_finished_ns=finished))


def run(output):
    job = int(os.environ["SLURM_JOB_ID"])
    if os.uname().nodename != "spark-7ff0" or os.environ.get("SLURM_CPUS_PER_TASK") != "2":
        raise ValueError("Require scheduled two-CPU task on spark-7ff0")
    output.mkdir(exist_ok=False)
    source = Path(__file__).resolve()
    sources = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in
               (source, source.with_name("probe_host_counters.py"), source.with_name("probe_dgx_step_separation.py"))}
    trials = []
    for index, monitor in enumerate(MODES):
        directory = output / str(index)
        directory.mkdir()
        command = ["srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1", "--cpus-per-task=1",
                   sys.executable, "-I", "-B", str(source), "--worker", str(directory)]
        with (directory / "native.log").open("x") as log:
            process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
            try:
                wait_file(directory / "ready.json")
                before = json.loads((directory / "native_before.json").read_text())
                samples = [snapshot()]
                validate_scopes(samples[0], before, job)
                observer_start = time.process_time()
                save(directory / "go.json", {"go": True})
                deadline = time.monotonic() + 400
                while process.poll() is None:
                    if time.monotonic() > deadline:
                        raise TimeoutError("Fixed-work step exceeded bound")
                    time.sleep(INTERVAL)
                    if monitor:
                        samples.append(snapshot())
                observer_cpu = time.process_time() - observer_start
                if process.returncode != 0:
                    raise RuntimeError("Native fixed-work step failed")
                trial = json.loads((directory / "done.json").read_text())
                trial.update(monitor=monitor, command=command, native_exit_code=process.returncode,
                    observer_cpu_s=observer_cpu, snapshots=samples, clock_ticks_per_second=os.sysconf("SC_CLK_TCK"),
                    native_snapshots=[before, json.loads((directory / "native_after.json").read_text())])
                trial["summary"] = validate_trial(trial, job)
                save(directory / "trial.json", trial)
                trials.append(trial)
            finally:
                if process.poll() is None:
                    if not (directory / "go.json").exists():
                        save(directory / "go.json", {"go": True, "failed_observation": True})
                    process.wait(timeout=410)
    for name, digest in sources.items():
        if hashlib.sha256(source.with_name(name).read_bytes()).hexdigest() != digest:
            raise ValueError("Source changed")
    wall = {str(mode): statistics.median(t["summary"]["work_wall_s"] for t in trials if t["monitor"] == mode)
            for mode in (False, True)}
    report = dict(status="fixed_work_counter_engineering_complete", job_id=job, sources=sources, trials=trials,
        median_work_wall_s=wall, descriptive_monitored_wall_change_percent=100 * (wall["True"] / wall["False"] - 1),
        publication_ready=False, controlled_workload_verified=False,
        limitations=["Six fixed-order trials on one-CPU hashing workloads, not native orthology inference.",
            "Wall change is descriptive, without confidence intervals or an overhead acceptance threshold.",
            "Native memory.peak is cgroup memory, not RSS or a calibrated peak accuracy measurement.",
            "Counters do not prove exclusivity; no host-minus-native foreign-load estimate is admitted.",
            "Does not upgrade the historical 27 timing runs or authorize a new scaling panel."])
    save(output / "report.json", report)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--load", action="store_true")
    mode.add_argument("--worker", type=Path)
    mode.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.load:
        print(json.dumps(load()))
    elif args.worker:
        worker(args.worker)
    else:
        run(args.output.absolute())
