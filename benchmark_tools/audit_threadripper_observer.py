"""Audit a bounded descendant-accounting diagnostic, never admit tool timings."""

import argparse
import json
import math
from pathlib import Path

from benchmark_tools.derive_threadripper_resources import derive
from benchmark_tools.run_integrated_full_job import record, save
from benchmark_tools.run_integrated_publication_workflow import check

DIMENSIONS = dict(workers=32, threads_per_worker=4, arena_mib_per_worker=64,
                  duration_seconds=30)
LIMITS = dict(min_complete_samples=5, max_point_seconds=1., max_sample_gap_seconds=1.5,
              cpu_rounding_tolerance_seconds=.01, cpu_excess_absolute_seconds=10.,
              cpu_excess_fraction=.05)


def verify_protocol(path, digest):
    if record(path)["sha256"] != digest:
        raise ValueError("Calibration protocol digest differs")
    value = json.loads(path.read_text())
    if (value["schema"] != "threadripper_observer_calibration_v1"
            or value["dimensions"] != DIMENSIONS or value["limits"] != LIMITS
            or value["production_timing"] is not False):
        raise ValueError("Calibration protocol differs")
    for pin in [value["resource_protocol"], value["interpreter"], *value["sources"]]:
        check(pin)
    required = (record(Path(__file__)), record(Path(__file__).with_name("calibrate_threadripper_observer.py")))
    # Compare full pins, not merely path labels supplied by the protocol.
    if not all(pin in value["sources"] for pin in required):
        raise ValueError("Calibration sources omitted")
    return value


def positive_integer(value):
    return type(value) is int and value > 0


def evaluate(witness, scoped, points, done, *, dimensions=DIMENSIONS):
    if any(witness[k] != v for k, v in dimensions.items()):
        raise ValueError("Witness workload dimensions differ")
    rows = witness["witnesses"]
    if len(rows) != dimensions["workers"] or witness["scientific_timings_admitted"] is not False:
        raise ValueError("Incomplete worker witnesses")
    expected, cpu_seconds, allocation_bytes = {}, 0., 0
    scope = scoped["cpu"]["scope"]
    for index, row in enumerate(rows):
        identity = row["identity"]
        membership = f"0::{scope}\n"
        if (identity["index"] != index or identity["membership"] != membership
                or row["final_membership"] != membership
                or identity["affinity"] != list(range(32))
                or row["final_affinity"] != identity["affinity"]
                or len(identity["threads"]) != dimensions["threads_per_worker"]
                or identity["arena_bytes"] != dimensions["arena_mib_per_worker"] * 1024**2):
            raise ValueError("Worker identity, scope, affinity or allocation differs")
        if (not positive_integer(row["started_ns"]) or not positive_integer(row["finished_ns"])
                or not done["started_ns"] <= row["started_ns"] < row["finished_ns"] <= done["finished_ns"]
                or (row["finished_ns"] - row["started_ns"]) / 1e9 < dimensions["duration_seconds"]):
            raise ValueError("Invalid worker interval")
        allocation_bytes += identity["arena_bytes"]
        for key in ("user_seconds", "system_seconds"):
            value = row[key]
            if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
                raise ValueError("Invalid worker CPU witness")
            cpu_seconds += value
        for thread in [dict(tid=identity["pid"], start_ticks=identity["start_ticks"],
                            membership=membership, affinity=identity["affinity"]), *identity["threads"]]:
            if (not positive_integer(thread["tid"]) or not positive_integer(thread["start_ticks"])
                    or thread["tid"] in expected or thread["membership"] != membership
                    or thread["affinity"] != identity["affinity"]):
                raise ValueError("Invalid or duplicated thread witness")
            expected[thread["tid"]] = thread
    common_start = max(row["started_ns"] for row in rows)
    common_finish = min(row["finished_ns"] for row in rows)
    costs, interior_times, complete = [], [], 0
    thread_errors = []
    for point in points:
        observed = point["thread_affinity"]
        start = point["host"][0]["started_monotonic_ns"]
        finish = observed["finished_ns"]
        if not positive_integer(start) or not positive_integer(finish) or finish < start:
            raise ValueError("Invalid point read interval")
        costs.append((finish - start) / 1e9)
        if common_start <= start <= finish <= common_finish:
            interior_times.append(start / 1e9)
            tids = {row["tid"]: row for row in observed["threads"]}
            valid = (observed["status"] == "observed_within_affinity" and not observed["errors"]
                     and not observed["violating_tids"]
                     and len(tids) == len(observed["threads"])
                     and all(tid in tids and tids[tid]["start_ticks"] == thread["start_ticks"]
                             and tids[tid]["cgroup"] == thread["membership"]
                             and tids[tid]["affinity"] == thread["affinity"]
                             and not tids[tid]["outside_cpus"] for tid, thread in expected.items()))
            if valid:
                complete += 1
            else:
                thread_errors.append(start)
    excess = scoped["primary"]["cpu_seconds"] - cpu_seconds
    gap = max((b - a for a, b in zip(interior_times, interior_times[1:])), default=math.inf)
    checks = dict(native_exited_zero=scoped["native_outcome"] == "exited_zero"
                  and scoped["native_exit_code"] == 0,
                  all_worker_threads_sampled=complete >= LIMITS["min_complete_samples"] and not thread_errors,
                  point_cost_within_cadence=bool(costs) and max(costs) <= LIMITS["max_point_seconds"],
                  sample_gap_within_limit=gap <= LIMITS["max_sample_gap_seconds"],
                  cpu_contains_worker_witness=excess >= -LIMITS["cpu_rounding_tolerance_seconds"],
                  cpu_excess_bounded=excess <= max(LIMITS["cpu_excess_absolute_seconds"],
                                                 LIMITS["cpu_excess_fraction"] * cpu_seconds),
                  peak_contains_allocations=scoped["primary"]["peak_memory_bytes"] >= allocation_bytes)
    return dict(status="calibration_checks_passed" if all(checks.values()) else "calibration_checks_failed",
                checks=checks, expected_threads=len(expected), complete_interior_samples=complete,
                incomplete_interior_samples=thread_errors, point_cost_seconds=costs,
                maximum_interior_sample_gap_seconds=gap if math.isfinite(gap) else None,
                witnessed_worker_cpu_seconds=cpu_seconds, native_cpu_excess_seconds=excess,
                simultaneously_retained_allocation_bytes=allocation_bytes,
                common_worker_interval_ns=[common_start, common_finish],
                primary=scoped["primary"], controlled_timing_admitted=False,
                full_run_containment_verified=False, slowdown_overhead_validated=False,
                publication_ready=False,
                limitations=["Python threads share a GIL; this tests 32 process workers and thread enumeration, not 128-core scaling.",
                    "CPU witnesses cover worker intervals only; the uncorrected native bracket also includes startup and wrappers.",
                    "Read durations measure collection cost, not causal workload slowdown or observer CPU.",
                    "Known allocations provide a peak-memory lower bound, not exact process RSS or an overhead-subtracted peak.",
                    "Periodic evidence is not continuous containment or quiet-host certification."])


def audit(directory, protocol_path, digest, output):
    if output.exists():
        raise FileExistsError(output)
    protocol = verify_protocol(protocol_path, digest)
    started_path = directory / "started.json"
    started = json.loads(started_path.read_text())
    launcher = record(Path(__file__).with_name("calibrate_threadripper_observer.py"))
    expected_command = [protocol["interpreter"]["path"], "-B", launcher["path"],
                        "--workload", "--output", str(directory / "witnesses")]
    if (started["protocol"] != record(protocol_path) or started["source"] != launcher
            or started["command"] != expected_command or started["production_timing"] is not False):
        raise ValueError("Launch and prospective protocol disagree")
    witness_path = directory / "witnesses/complete.json"
    witness = json.loads(witness_path.read_text())
    if witness["source"] != launcher:
        raise ValueError("Witness source differs")
    measurement = directory / "measurement"
    if json.loads((measurement / "command.json").read_text())["command"] != expected_command:
        raise ValueError("Measured command differs")
    resource_output = output.with_name(output.name + ".scoped.json")
    scoped = derive(measurement, started["job_id"], Path(protocol["resource_protocol"]["path"]),
                    protocol["resource_protocol"]["sha256"], resource_output)
    points = [json.loads(Path(pin["path"]).read_text())
              for pin in scoped["evidence"] if Path(pin["path"]).name.startswith("point_")]
    points.sort(key=lambda point: point["thread_affinity"]["started_ns"])
    result = evaluate(witness, scoped, points, json.loads((measurement / "done.json").read_text()))
    result.update(job_id=started["job_id"], protocol=record(protocol_path), source=record(__file__),
                  evidence=[record(started_path), record(witness_path), record(resource_output), *scoped["evidence"]])
    for pin in result["evidence"]:
        check(pin)
    verify_protocol(protocol_path, digest)
    save(output, result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--protocol-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.directory.resolve(), args.protocol.resolve(), args.protocol_sha256, args.output.resolve())
    print(result["status"])
    raise SystemExit(0 if result["status"] == "calibration_checks_passed" else 1)
