"""Replay quiet-control counters and describe CPU fields without attribution."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.diagnose_dgx_interval_residuals import decompose, delta_stat
from benchmark_tools.measure_native_hierarchy_step import evaluate, interval_point
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check

INPUT_SHA = "29d0268a44d77c23e4822a2a062cc9d1aa4fffe88279885e82aab0a32835441d"


def diagnose(measurement):
    points, job = measurement["points"], measurement["job_id"]
    replay = evaluate(points, measurement["native"], job)
    if replay != measurement["screening"]:
        raise ValueError("Stored screening differs from counter replay")
    rows = []
    for index, (left, right) in enumerate(zip(points, points[1:])):
        row = decompose(interval_point(left, job), interval_point(right, job), job)
        parent = delta_stat(left["parent"][0]["raw"], right["parent"][1]["raw"])
        fields = row["host_outer_field_cpu_s"]
        # These accounting systems differ; subtraction is diagnostic, not attribution.
        residual = dict(
            user_and_nice=fields["user"] + fields["nice"] - parent["user_usec"],
            system=fields["system"] - parent["system_usec"],
            irq_and_softirq=fields["irq"] + fields["softirq"],
        )
        residual_sum = sum(residual.values())
        direct = row["host_outer_busy_cpu_s"] - parent["usage_usec"]
        accounting_gap = parent["usage_usec"] - parent["user_usec"] - parent["system_usec"]
        if abs(residual_sum - direct - accounting_gap) > 1e-8:
            raise ValueError("CPU field decomposition identity failed")
        rows.append(dict(index=index, **row, job_outer_delta_cpu_s=parent,
                         diagnostic_field_residual_cpu_s=residual,
                         job_usage_minus_user_system_cpu_s=accounting_gap,
                         host_outer_minus_job_outer_cpu_s=direct))
    return dict(job_id=job, screening=replay, intervals=rows)


def run(source, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    identity = record(source)
    if identity["sha256"] != INPUT_SHA:
        raise ValueError("Changed retained quiet-control evidence")
    report = json.loads(source.read_text())
    if len(report["runs"]) != 3:
        raise ValueError("Require all three retained methods")
    rows = [dict(run_index=i, **diagnose(row["verification"]["measurement"]))
            for i, row in enumerate(report["runs"])]
    check(identity)
    result = dict(status="retrospective_quiet_cpu_field_decomposition", input=identity,
        source=record(__file__), runs=rows,
        helpers=[record(Path(__file__).with_name(name)) for name in (
            "diagnose_dgx_interval_residuals.py", "measure_native_hierarchy_step.py",
            "measure_native_interval_step.py", "probe_dgx_cpu_hierarchy.py",
            "probe_host_counters.py", "probe_interval_cpu.py", "screen_bracketed_cpu.py",
            "probe_dgx_step_separation.py", "prepare_ob_candidate_neighborhood.py")],
        scientific_timings_admitted=False, controlled_workload_verified=False,
        limitations=[
            "Retrospective accounting diagnostic across all three methods, not a timing correction.",
            "Host and cgroup counters differ in scope, quantization and update latency.",
            "Interrupt time cannot automatically be assigned to unrelated work.",
            "User/system residuals do not identify a process or prove external interference.",
            "Signed residuals retained; overlapping outer windows must not be summed.",
            "No threshold changes, new benchmark executions or scientific timing admission."])
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    run(args.source.resolve(), args.output.absolute())
