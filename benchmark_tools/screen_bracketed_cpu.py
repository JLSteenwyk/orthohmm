"""Operational CPU screen requiring outer host and inner native read windows.

This is not a statistical interference bound or scientific timing admission.
"""

import math
import hashlib
import json
from pathlib import Path

from benchmark_tools.probe_host_counters import summarize
from benchmark_tools.probe_dgx_step_separation import validate_scopes

MAX_RESIDUAL_CORES = .25
NEGATIVE_TOLERANCE_CPU_S = .5


def usage(sample):
    rows = [row.split() for row in sample["optional"]["cgroup_cpu.stat"].splitlines()]
    if any(len(row) != 2 for row in rows):
        raise ValueError("Malformed cgroup CPU counter")
    keys = [row[0] for row in rows]
    if len(set(keys)) != len(keys) or keys.count("usage_usec") != 1:
        raise ValueError("Missing or duplicate CPU usage field")
    values = {key: int(value) for key, value in rows}
    if any(value < 0 for value in values.values()):
        raise ValueError("Negative cgroup CPU counter")
    return values["usage_usec"]


def screen(host_before, native_before, native_after, host_after, job_id, ticks_per_second,
           work_started_ns, work_finished_ns):
    samples = (host_before, native_before, native_after, host_after)
    if any(sample["errors"] for sample in samples):
        raise ValueError("Counter read failures cannot support screening")
    for key in ("boot_id", "online_cpus"):
        if len({sample["raw"][key] for sample in samples}) != 1:
            raise ValueError("Host identity/topology changed")
    host = summarize(host_before, host_after, ticks_per_second)
    summarize(native_before, native_after, ticks_per_second)
    validate_scopes(host_before, native_before, job_id)
    validate_scopes(host_after, native_after, job_id)
    if not (host_before["finished_monotonic_ns"] <= native_before["started_monotonic_ns"]
            and native_after["finished_monotonic_ns"] <= host_after["started_monotonic_ns"]):
        raise ValueError("Host read window must fully bracket native read window")
    if (type(work_started_ns) is not int or type(work_finished_ns) is not int
            or not native_before["finished_monotonic_ns"] < work_started_ns < work_finished_ns < native_after["started_monotonic_ns"]):
        raise ValueError("Native reads must bracket command work")
    native_usec = usage(native_after) - usage(native_before)
    if native_usec < 0:
        raise ValueError("Native CPU counter decreased")
    native_cpu = native_usec / 1e6
    work_s = (work_finished_ns - work_started_ns) / 1e9
    residual = host["accounted_host_busy_cpu_s"] - native_cpu
    cores = residual / work_s
    if not all(math.isfinite(v) for v in (native_cpu, work_s, residual, cores)):
        raise ValueError("Nonfinite counter arithmetic")
    reasons = []
    if residual < -NEGATIVE_TOLERANCE_CPU_S:
        reasons.append("negative_accounting_discrepancy")
    if cores > MAX_RESIDUAL_CORES:
        reasons.append("excess_unassigned_cpu")
    if host["cpu_delta_ticks"]["steal"]:
        reasons.append("host_steal_time")
    return dict(status="bracketed_cpu_operational_screen", screen_passed=not reasons, reasons=reasons,
        host_busy_cpu_s=host["accounted_host_busy_cpu_s"], native_cpu_s=native_cpu,
        signed_unassigned_cpu_s=residual, signed_unassigned_average_cores=cores, work_wall_s=work_s,
        maximum_residual_cores=MAX_RESIDUAL_CORES, negative_tolerance_cpu_s=NEGATIVE_TOLERANCE_CPU_S,
        controlled_workload_verified=False, scientific_timings_admitted=False,
        limitations=["Operational thresholds, not calibrated probabilistic or deterministic interference bounds.",
            "Unassigned CPU includes observer, kernel and launch/exit work; it is not identified foreign CPU.",
            "Counter accounting granularity/delay and non-CPU interference remain unresolved.",
            "A run-wide average can hide concentrated bursts; interval-level monitoring remains necessary.",
            "No subtraction from native wall time or correction of performance results is permitted."])


def audit_existing(archive):
    trials = []
    for index in range(3):
        directory = archive / "counter_native_smoke_v1" / f"run_{index:02d}" / "measurement"
        report_path, stream_path = directory / "counter_report.json", directory / "observer.jsonl"
        sources = []
        for path in (report_path, stream_path):
            data = path.read_bytes()
            sources.append(dict(path=str(path), bytes=len(data), sha256=hashlib.sha256(data).hexdigest()))
        report = json.loads(report_path.read_text())
        observations = [json.loads(line) for line in stream_path.read_text().splitlines()]
        if len(observations) != report["observer_snapshots"]:
            raise ValueError("Observer stream length differs")
        native_before, native_after = report["native_snapshots"]
        done = report["native"]
        snapshots = [observations[0], native_before, native_after, observations[-1]]
        try:
            screen(*snapshots, report["job_id"], 100, done["started_ns"], done["finished_ns"])
        except ValueError as error:
            reason = str(error)
            if reason != "Host read window must fully bracket native read window":
                raise
        else:
            raise ValueError("Historical smoke unexpectedly passed bracketing")
        for item in sources:
            if hashlib.sha256(Path(item["path"]).read_bytes()).hexdigest() != item["sha256"]:
                raise ValueError("Evidence changed during audit")
        trials.append(dict(index=index, sources=sources, snapshots=snapshots, job_id=report["job_id"],
            work_started_ns=done["started_ns"], work_finished_ns=done["finished_ns"],
            left_overlap_ns=observations[0]["finished_monotonic_ns"] - native_before["started_monotonic_ns"],
            rejection=reason))
    return dict(status="historical_native_windows_rejected", trials=trials,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        clock_ticks_per_second=100, controlled_workload_verified=False, scientific_timings_admitted=False,
        limitations=["Actual recorded boundaries, not a contention estimate.",
                     "No historical timing is changed, corrected or newly admitted."])


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-existing", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit_existing(args.audit_existing.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
