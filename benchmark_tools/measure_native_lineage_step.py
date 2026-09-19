"""Prospective native collector with aggregate lineage CPU and unchanged screens."""

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
from benchmark_tools.probe_host_counters import snapshot as host_snapshot, summarize
from benchmark_tools.probe_cgroup_lineage import snapshot as lineage_snapshot, validate as validate_lineage
from benchmark_tools.probe_cgroup_lineage import compare as compare_lineage, LineageSnapshotError
from benchmark_tools.probe_native_pressure import read_point as read_pressure, validate as validate_pressure
from benchmark_tools.probe_native_pressure import compare as compare_pressure
from benchmark_tools.probe_interval_cpu import interval
from benchmark_tools.probe_dgx_step_separation import save, wait_file


def validate_point(point, job):
    if point.get("schema") != "native_lineage_v1" or "frontier" in point:
        raise ValueError("Require distinct native lineage schema, not frontier evidence")
    if any(key not in point for key in ("lineage", "native_pressure", "hierarchy_host_after")):
        raise ValueError("Incomplete lineage observation")
    interval_point(point, job)
    lineage = point["lineage"]
    validate_lineage(lineage)
    old, new = point["hierarchy_host_after"], point["host"][1]
    summarize(old, new, point["ticks"])
    if (lineage["target"] != point["parent"][0]["scope"]
            or lineage["boot_before"] != new["raw"]["boot_id"].strip()
            or usage(lineage["rows"][-1]["raw"]) < usage(point["parent"][1]["raw"])):
        raise ValueError("Lineage target, boot or parent counter differs")
    pressure = point["native_pressure"]
    validate_pressure(pressure, job)
    summarize(old, pressure["host"][1], point["ticks"])
    if pressure["native_membership"] != point["native_membership"] or pressure["ticks"] != point["ticks"]:
        raise ValueError("Pressure scope differs from native observation")
    if not (point["parent"][1]["finished_ns"] <= old["started_monotonic_ns"]
            <= old["finished_monotonic_ns"] <= lineage["rows"][0]["started_ns"]
            <= lineage["rows"][-1]["finished_ns"] <= pressure["host"][0]["started_monotonic_ns"]
            <= pressure["host"][1]["finished_monotonic_ns"] <= new["started_monotonic_ns"]):
        raise ValueError("Lineage and pressure reads are not enclosed and ordered")
    return interval_point({**point, "host": [point["host"][0], old]}, job)


def read_point(pid, membership, job, failure_path=None):
    point = dict(schema="native_lineage_v1")
    try:
        point.update(read_hierarchy(pid, membership, job))
        point["hierarchy_host_after"] = point["host"][1]
        point["lineage"] = lineage_snapshot(Path("/sys/fs/cgroup"), point["parent"][0]["scope"])
        point["native_pressure"] = read_pressure(pid, membership, job)
        point["host"][1] = host_snapshot()
        if Path(f"/proc/{pid}/cgroup").read_text() != membership:
            raise ValueError("Native process disappeared or changed scope")
        validate_point(point, job)
    except (OSError, ValueError, KeyError) as error:
        if failure_path is not None:
            failure = dict(status="invalid_native_lineage_observation", preceding_point=point,
                error=str(error), error_type=type(error).__name__, scientific_timings_admitted=False)
            if isinstance(error, LineageSnapshotError):
                failure["partial_lineage"] = error.evidence
            save(failure_path, failure)
        raise
    return point


def compare(left, right, job, *, enforce_gap=True):
    narrow_left, narrow_right = (validate_point(p, job) for p in (left, right))
    outer = interval(interval_point(left, job), interval_point(right, job), job, enforce_gap=enforce_gap)
    narrow = interval(narrow_left, narrow_right, job, enforce_gap=enforce_gap)
    if outer["native_cpu_s"] != narrow["native_cpu_s"] or outer["wall_s"] != narrow["wall_s"]:
        raise ValueError("Paired native counter measurements differ")
    return dict(status="native_lineage_interval_v1", outer=outer, narrow=narrow,
        lineage=compare_lineage(left["lineage"], right["lineage"]),
        native_pressure=compare_pressure(left["native_pressure"], right["native_pressure"], job),
        scientific_timings_admitted=False, controlled_workload_verified=False)


def evaluate(points, done, job):
    for point in points:
        validate_point(point, job)
    original = hierarchy_evaluate(points, done, job)
    pairs = [compare(a, b, job) for a, b in zip(points, points[1:])]
    if [p["outer"] for p in pairs] != original["original_threshold_screen"]["intervals"]:
        raise ValueError("Lineage outer intervals differ from original evaluator")
    return dict(schema="native_lineage_screening_v1", original_screening=original,
        intervals=pairs, narrow_flagged_intervals=[i for i, p in enumerate(pairs) if not p["narrow"]["screen_passed"]],
        observation_window=compare(points[0], points[-1], job, enforce_gap=False),
        scientific_timings_admitted=False, controlled_workload_verified=False,
        limitations=["Both original CPU thresholds retained; neither screen admits scientific timing.",
            "Lineage signed differences are descriptive, non-atomic and not foreign-CPU bounds.",
            "Observation endpoints include wrapper activity, not exact native command boundaries.",
            "No frontier substitution in historical reports, timing correction or overhead admission."])


def measure(command, directory, job_id, cpus, memory_bytes, timeout_s, interval_s,
            monitor_host=True, host_interval_s=30., *, point_reader=None):
    validate(command, cpus, timeout_s, interval_s)
    reader = read_point if point_reader is None else point_reader
    if not callable(reader):
        raise ValueError("Point reader must be callable")
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
            points = [reader(ready["pid"], ready["cgroup"], job_id, directory / "failed_point.json")]
            save(directory / "point_0000.json", points[0])
            save(directory / "go.json", {"go": True})
            start = time.monotonic()
            index = 0
            while True:
                index += 1
                time.sleep(max(0, start + index * interval_s - time.monotonic()))
                completed = (directory / "done.json").exists()
                if process.poll() is not None:
                    raise RuntimeError("Worker exited before final lineage observation")
                points.append(reader(ready["pid"], ready["cgroup"], job_id, directory / "failed_point.json"))
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
                limitations=["Complete-command lineage engineering measurement, not controlled comparative timing.",
                             "Original screens and failures retained; no threshold changes or wall-time correction.",
                             "Native-step memory peak includes wrappers and cache, not maximum process RSS.",
                             "Full-node behavior, overhead and scientific inclusion remain unvalidated."])
            save(directory / "lineage_report.json", result)
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
