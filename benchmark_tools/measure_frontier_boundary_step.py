"""Boundary-only control for incremental periodic frontier observation cost."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.measure_native_counter_step import validate
from benchmark_tools.measure_native_interval_step import worker, step_memory
from benchmark_tools.measure_native_frontier_step import read_frontier_point, validate_point, interval_point
from benchmark_tools.probe_dgx_cpu_hierarchy import compare
from benchmark_tools.probe_cgroup_frontier import compare as compare_frontier
from benchmark_tools.probe_dgx_step_separation import save, wait_file
from benchmark_tools.screen_bracketed_cpu import screen


def evaluate(points, done, job):
    if len(points) != 2:
        raise ValueError("Boundary control requires exactly two observations")
    for point in points:
        validate_point(point, job)
    if not (points[0]["host"][1]["finished_monotonic_ns"] < done["started_ns"]
            < done["finished_ns"] < points[-1]["host"][0]["started_monotonic_ns"]):
        raise ValueError("Boundary observations do not enclose complete command")
    before, after = done["snapshots"]
    whole = screen(points[0]["host"][0], before, after, points[-1]["host"][1],
                   job, points[0]["ticks"], done["started_ns"], done["finished_ns"])
    return dict(whole_command_screen=whole,
                hierarchy_boundary=compare(points[0], points[1], job),
                frontier_boundary=compare_frontier(points[0]["frontier"], points[1]["frontier"]),
                interval_screening_available=False, flagged_intervals=None,
                scientific_timings_admitted=False, controlled_workload_verified=False,
                limitations=["Boundary-only observations cannot detect localized competing activity.",
                             "Same worker and completion polling; estimates incremental periodic observation cost only.",
                             "No quiet-host assertion, overhead subtraction, or scientific timing admission."])


def measure(command, directory, job_id, cpus, memory_bytes, timeout_s, interval_s,
            monitor_host=True, host_interval_s=30.):
    validate(command, cpus, timeout_s, interval_s)
    if (int(os.environ["SLURM_JOB_ID"]) != job_id or os.environ.get("SLURM_CPUS_PER_TASK") != "20"
            or os.environ.get("SLURM_MEM_PER_NODE") != "98304"
            or os.uname().nodename != "spark-7ff0" or memory_bytes != 96 * 1024 ** 3):
        raise ValueError("Require matching DGX20CPU/96GiB allocation")
    directory = Path(directory).absolute()
    directory.mkdir(exist_ok=False)
    save(directory / "command.json", dict(command=command, cpus=cpus, timeout_s=timeout_s, interval_s=interval_s))
    launched = ["srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1", "--cpus-per-task=20",
                sys.executable, "-B", str(Path(__file__).resolve()), "--worker", str(directory)]
    with (directory / "step.log").open("x") as log:
        process = subprocess.Popen(launched, stdout=log, stderr=subprocess.STDOUT)
        try:
            ready = wait_file(directory / "ready.json")
            points = [read_frontier_point(ready["pid"], ready["cgroup"], job_id)]
            save(directory / "point_0000.json", points[0])
            save(directory / "go.json", {"go": True})
            start = time.monotonic()
            index = 0
            while True:
                index += 1
                time.sleep(max(0, start + index * interval_s - time.monotonic()))
                completed = (directory / "done.json").exists()
                if process.poll() is not None:
                    raise RuntimeError("Worker exited before final frontier observation")
                if completed:
                    points.append(read_frontier_point(ready["pid"], ready["cgroup"], job_id))
                    save(directory / f"point_{index:04d}.json", points[-1])
                    break
                if time.monotonic() - start > timeout_s + 30:
                    raise TimeoutError("Native command exceeded timeout and cleanup allowance")
            done = json.loads((directory / "done.json").read_text())
            memory = step_memory(interval_point(points[-1], job_id))
            save(directory / "step_memory.json", memory)
            save(directory / "release.json", {"release": True})
            if process.wait(timeout=45) != 0:
                raise RuntimeError("Native step wrapper failed")
            result = dict(status="command_exited_zero" if done["exit_code"] == 0 else "command_failed",
                native=done, native_wall_s=(done["finished_ns"] - done["started_ns"]) / 1e9,
                job_id=job_id, launched=launched, points=points, step_memory=memory,
                screening=evaluate(points, done, job_id), scientific_timings_admitted=False,
                controlled_workload_verified=False, publication_ready=False,
                limitations=["Boundary-only incremental-overhead control, not controlled comparative timing.",
                    "No periodic interval screening is available; absence of flags is not a quiet-host result.",
                    "Host/parent/step reads are non-atomic; counters may have different accounting delays.",
                    "Native-step memory peak includes wrappers and cache, not maximum process RSS.",
                    "Additional counter reads can affect observation overhead; no overhead correction is assumed.",
                    "Non-CPU interference and scientific inclusion policy remain unresolved."])
            save(directory / "boundary_report.json", result)
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
