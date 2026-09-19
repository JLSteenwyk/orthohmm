"""Complete native measurement with a separate outside-job cgroup frontier."""

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
from benchmark_tools.measure_native_hierarchy_step import interval_point, evaluate as hierarchy_evaluate
from benchmark_tools.probe_dgx_cpu_hierarchy import read_point as read_hierarchy, usage
from benchmark_tools.probe_dgx_step_separation import save, wait_file
from benchmark_tools.probe_host_counters import snapshot as host_snapshot, summarize
from benchmark_tools.probe_cgroup_frontier import snapshot as frontier_snapshot, validate as validate_frontier, compare as compare_frontier
from benchmark_tools.probe_cgroup_frontier import FrontierSnapshotError


def validate_point(point, job):
    interval_point(point, job)
    frontier = point["frontier"]
    validate_frontier(frontier)
    old_after = point["hierarchy_host_after"]
    new_after = point["host"][1]
    summarize(old_after, new_after, point["ticks"])
    if (frontier["target"] != point["parent"][0]["scope"]
            or frontier["boot_id"] != new_after["raw"]["boot_id"].strip()):
        raise ValueError("Frontier target or boot differs from native observation")
    target = next(row for row in frontier["rows"] if row["scope"] == frontier["target"])
    if usage(target["raw"]) < usage(point["parent"][1]["raw"]):
        raise ValueError("Job CPU decreased before frontier read")
    if not (point["parent"][1]["finished_ns"] <= old_after["started_monotonic_ns"]
            <= old_after["finished_monotonic_ns"] <= frontier["root"][0]["started_ns"]
            <= frontier["root"][1]["finished_ns"] <= new_after["started_monotonic_ns"]):
        raise ValueError("Frontier counters are not enclosed by host brackets")


def read_frontier_point(pid, membership, job, failure_path=None):
    point = read_hierarchy(pid, membership, job)
    point["hierarchy_host_after"] = point["host"][1]
    try:
        point["frontier"] = frontier_snapshot(Path("/sys/fs/cgroup"), point["parent"][0]["scope"])
    except FrontierSnapshotError as error:
        if failure_path is not None:
            save(failure_path, dict(status="invalid_frontier_observation", hierarchy=point,
                                   frontier=error.evidence, scientific_timings_admitted=False))
        raise
    point["host"][1] = host_snapshot()
    if Path(f"/proc/{pid}/cgroup").read_text() != membership:
        raise ValueError("Native process disappeared or changed scope")
    validate_point(point, job)
    return point


def evaluate(points, done, job):
    for point in points:
        validate_point(point, job)
    result = hierarchy_evaluate(points, done, job)
    result["frontier_intervals"] = [compare_frontier(a["frontier"], b["frontier"])
                                    for a, b in zip(points, points[1:])]
    result["frontier_limitations"] = [
        "Frontier and native counters have different read windows; do not subtract them as synchronized measurements.",
        "Outside-scope CPU is descriptive accounting, not causal interference or a native wall-time correction.",
        "Expanded host brackets include the additional collector reads; original thresholds remain unchanged.",
        "Collector overhead and non-CPU isolation remain unvalidated."
    ]
    return result


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
            points = [read_frontier_point(ready["pid"], ready["cgroup"], job_id, directory / "failed_point.json")]
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
                points.append(read_frontier_point(ready["pid"], ready["cgroup"], job_id, directory / "failed_point.json"))
                save(directory / f"point_{index:04d}.json", points[-1])
                if completed:
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
                limitations=["Complete-command frontier engineering test, not controlled comparative timing.",
                    "Original CPU thresholds retained; hierarchy residuals do not alter flags or correct wall time.",
                    "Host/parent/step reads are non-atomic; counters may have different accounting delays.",
                    "Native-step memory peak includes wrappers and cache, not maximum process RSS.",
                    "Additional counter reads can affect observation overhead; no overhead correction is assumed.",
                    "Non-CPU interference and scientific inclusion policy remain unresolved."])
            save(directory / "frontier_report.json", result)
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
