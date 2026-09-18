"""Prospective outer-host/inner-native CPU controls on the dedicated DGX."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.probe_host_counters import snapshot, parse_group
from benchmark_tools.probe_dgx_step_separation import save, wait_file
from benchmark_tools.screen_bracketed_cpu import screen

MODES = (("quiet", 0), ("completed_burst", 1), ("sustained", 3))


def burn(seconds):
    if seconds not in (1, 2, 3):
        raise ValueError("Only frozen engineering CPU durations")
    start = time.process_time()
    value = 1
    while time.process_time() - start < seconds:
        for _ in range(1000):
            value = (value * 1664525 + 1013904223) % 4294967296
    return dict(cpu_s=time.process_time() - start, checksum=value,
                cgroup=Path("/proc/self/cgroup").read_text(), pid=os.getpid())


def worker(directory):
    save(directory / "ready.json", {"pid": os.getpid(), "cgroup": Path("/proc/self/cgroup").read_text()})
    wait_file(directory / "go.json")
    before = snapshot()
    started = time.monotonic_ns()
    load = burn(2)
    finished = time.monotonic_ns()
    after = snapshot()
    save(directory / "done.json", dict(snapshots=[before, after], load=load,
        work_started_ns=started, work_finished_ns=finished))


def evaluate(host_before, host_after, native, sibling, job, ticks):
    if not 2 <= native["load"]["cpu_s"] < 2.1:
        raise ValueError("Unexpected native positive-control CPU duration")
    before, after = native["snapshots"]
    if parse_group(native["load"]["cgroup"]) != parse_group(before["raw"]["cgroup_membership"]):
        raise ValueError("Native load changed scope")
    if sibling is not None:
        if not sibling["requested_cpu_s"] <= sibling["cpu_s"] < sibling["requested_cpu_s"] + .1:
            raise ValueError("Unexpected sibling control CPU duration")
        if parse_group(sibling["cgroup"]) != parse_group(host_before["raw"]["cgroup_membership"]):
            raise ValueError("Sibling load escaped observer batch scope")
    return screen(host_before, before, after, host_after, job, ticks,
                  native["work_started_ns"], native["work_finished_ns"])


def run(output):
    job = int(os.environ["SLURM_JOB_ID"])
    if os.uname().nodename != "spark-7ff0" or os.environ.get("SLURM_CPUS_PER_TASK") != "2":
        raise ValueError("Require dedicated DGX two-CPU engineering task")
    output.mkdir(exist_ok=False)
    source = Path(__file__).resolve()
    names = (source.name, "probe_host_counters.py", "probe_dgx_step_separation.py", "screen_bracketed_cpu.py")
    hashes = {name: hashlib.sha256(source.with_name(name).read_bytes()).hexdigest() for name in names}
    trials = []
    for mode, seconds in MODES:
        directory = output / mode
        directory.mkdir()
        command = ["srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1", "--cpus-per-task=1",
                   sys.executable, "-I", "-B", str(source), "--worker", str(directory)]
        with (directory / "step.log").open("x") as log:
            process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
            try:
                ready = wait_file(directory / "ready.json")
                host_before = snapshot()
                save(directory / "go.json", {"go": True})
                sibling = None
                if seconds:
                    child = subprocess.run([sys.executable, "-I", "-B", str(source), "--burn", str(seconds)],
                                           check=True, capture_output=True, text=True, timeout=30)
                    sibling = dict(json.loads(child.stdout), requested_cpu_s=seconds)
                if process.wait(timeout=45) != 0:
                    raise RuntimeError("Native control failed")
                host_after = snapshot()
                native = json.loads((directory / "done.json").read_text())
                ticks = os.sysconf("SC_CLK_TCK")
                result = evaluate(host_before, host_after, native, sibling, job, ticks)
                trial = dict(mode=mode, ready=ready, command=command, native=native, sibling=sibling,
                    host_snapshots=[host_before, host_after], clock_ticks_per_second=ticks, screen=result,
                    expected_screen_pass=not bool(seconds), expectation_met=result["screen_passed"] == (not bool(seconds)))
                save(directory / "trial.json", trial)
                trials.append(trial)
            finally:
                if process.poll() is None:
                    if not (directory / "go.json").exists():
                        save(directory / "go.json", {"go": True, "failed_observation": True})
                    process.wait(timeout=60)
    if hashes != {name: hashlib.sha256(source.with_name(name).read_bytes()).hexdigest() for name in names}:
        raise ValueError("Probe source changed")
    report = dict(status="bracketed_cpu_controls_evaluated", job_id=job, sources=hashes, trials=trials,
        all_control_expectations_met=all(t["expectation_met"] for t in trials),
        controlled_workload_verified=False, scientific_timings_admitted=False, publication_ready=False,
        limitations=["Three fixed-order controls, not general calibration or repeated native pipeline timing.",
            "Positive loads are sibling children in our allocation, not unrelated user jobs.",
            "The sustained load may extend past native work; outer-window residual retains that excess deliberately.",
            "CPU-only controls do not test I/O, memory, thermal interference or long-run burst averaging."])
    save(output / "report.json", report)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--worker", type=Path)
    mode.add_argument("--burn", type=int)
    mode.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.worker:
        worker(args.worker)
    elif args.burn is not None:
        print(json.dumps(burn(args.burn)))
    else:
        run(args.output.absolute())
