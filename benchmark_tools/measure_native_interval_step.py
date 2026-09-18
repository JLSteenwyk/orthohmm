"""Complete-command interval collector for native engineering smokes."""

import argparse
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.measure_native_counter_step import validate
from benchmark_tools.probe_host_counters import snapshot
from benchmark_tools.probe_dgx_step_separation import save, wait_file
from benchmark_tools.probe_interval_cpu import read_point, interval
from benchmark_tools.screen_bracketed_cpu import screen


def run_command(command, log, timeout_s):
    with subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, start_new_session=True) as process:
        try:
            return process.wait(timeout=timeout_s), False
        except subprocess.TimeoutExpired:
            # This group belongs to our command, not the observer or other jobs.
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait()
            return 124, True


def worker(directory):
    plan = json.loads((directory / "command.json").read_text())
    validate(plan["command"], plan["cpus"], plan["timeout_s"], plan["interval_s"])
    save(directory / "ready.json", dict(pid=os.getpid(), cgroup=Path("/proc/self/cgroup").read_text()))
    wait_file(directory / "go.json")
    before = snapshot()
    started = time.monotonic_ns()
    with (directory / "native.log").open("x") as log:
        code, timed_out = run_command(plan["command"], log, plan["timeout_s"])
    finished = time.monotonic_ns()
    after = snapshot()
    save(directory / "done.json", dict(exit_code=code, timed_out=timed_out,
        started_ns=started, finished_ns=finished, snapshots=[before, after]))
    # Keep the native step aggregate available for the observer's final reads.
    wait_file(directory / "release.json")


def evaluate(points, done, job):
    if len(points) < 2:
        raise ValueError("Missing start or final observation")
    if not (points[0]["host"][1]["finished_monotonic_ns"] < done["started_ns"]
            < done["finished_ns"] < points[-1]["host"][0]["started_monotonic_ns"]):
        raise ValueError("Observation points do not enclose complete command")
    intervals = [interval(a, b, job) for a, b in zip(points, points[1:])]
    before, after = done["snapshots"]
    whole = screen(points[0]["host"][0], before, after, points[-1]["host"][1],
                   job, points[0]["ticks_per_second"], done["started_ns"], done["finished_ns"])
    return dict(intervals=intervals, whole_command_screen=whole,
        flagged_intervals=[i for i, item in enumerate(intervals) if not item["screen_passed"]],
        scientific_timings_admitted=False, controlled_workload_verified=False)


def step_memory(point):
    scope = point["native_cpu_scope"]
    directory = Path("/sys/fs/cgroup") / scope.lstrip("/")
    started = time.monotonic_ns()
    raw, errors = {}, []
    for name in ("memory.current", "memory.peak", "memory.events"):
        try:
            raw[name] = (directory / name).read_text()
        except OSError as error:
            errors.append(dict(field=name, type=type(error).__name__, errno=error.errno))
    return dict(scope=scope, started_ns=started, finished_ns=time.monotonic_ns(), raw=raw, errors=errors)


def measure(command, directory, job_id, cpus, memory_bytes, timeout_s, interval_s,
            monitor_host=True, host_interval_s=30.):
    validate(command, cpus, timeout_s, interval_s)
    if int(os.environ["SLURM_JOB_ID"]) != job_id or os.environ.get("SLURM_CPUS_PER_TASK") != "20":
        raise ValueError("Wrong scheduler allocation")
    if os.uname().nodename != "spark-7ff0" or memory_bytes != 96 * 1024 ** 3:
        raise ValueError("Wrong host or memory request")
    directory = Path(directory).absolute()
    directory.mkdir(exist_ok=False)
    save(directory / "command.json", dict(command=command, cpus=cpus, timeout_s=timeout_s, interval_s=interval_s))
    launched = ["srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1", "--cpus-per-task=20",
                sys.executable, "-B", str(Path(__file__).resolve()), "--worker", str(directory)]
    with (directory / "step.log").open("x") as log:
        process = subprocess.Popen(launched, stdout=log, stderr=subprocess.STDOUT)
        try:
            ready = wait_file(directory / "ready.json")
            points = [read_point(ready, job_id)]
            save(directory / "point_0000.json", points[0])
            save(directory / "go.json", {"go": True})
            start = time.monotonic()
            index = 0
            while True:
                index += 1
                time.sleep(max(0, start + index * interval_s - time.monotonic()))
                # Observe completion BEFORE the point to guarantee a post-command boundary.
                completed = (directory / "done.json").exists()
                if process.poll() is not None:
                    raise RuntimeError("Worker exited before final observation")
                points.append(read_point(ready, job_id))
                save(directory / f"point_{index:04d}.json", points[-1])
                if completed:
                    break
                if time.monotonic() - start > timeout_s + 30:
                    raise TimeoutError("Native command exceeded timeout and cleanup allowance")
            done = json.loads((directory / "done.json").read_text())
            memory = step_memory(points[-1])
            save(directory / "step_memory.json", memory)
            save(directory / "release.json", {"release": True})
            if process.wait(timeout=45) != 0:
                raise RuntimeError("Native step wrapper failed")
            screening = evaluate(points, done, job_id)
            result = dict(status="command_exited_zero" if done["exit_code"] == 0 else "command_failed",
                native=done, native_wall_s=(done["finished_ns"] - done["started_ns"]) / 1e9,
                job_id=job_id, launched=launched, points=points, step_memory=memory,
                screening=screening, scientific_timings_admitted=False,
                controlled_workload_verified=False, publication_ready=False,
                limitations=["Complete-command engineering smoke, not controlled comparative timing.",
                    "Observer shares allocation CPUs with the native step; scopes are separate.",
                    "Intervals use step aggregate CPU; whole-command native reads use worker cgroup including descendants.",
                    "Outer host intervals overlap: do not sum residuals or correct native runtime.",
                    "Final step memory peak includes wrapper/startup and cache, not process RSS.",
                    "Accounting/overhead, non-CPU interference and scientific inclusion policy remain unresolved."])
            save(directory / "interval_report.json", result)
            return result
        finally:
            for name in ("go.json", "release.json"):
                if not (directory / name).exists():
                    save(directory / name, {"cleanup": True})
            if process.poll() is None:
                process.wait(timeout=timeout_s + 90)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=Path, required=True)
    worker(parser.parse_args().worker)
