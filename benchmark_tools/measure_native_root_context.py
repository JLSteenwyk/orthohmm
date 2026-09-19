"""Add separately timed root context to the existing lineage engineering collector."""

from pathlib import Path

from benchmark_tools import measure_native_lineage_step as lineage
from benchmark_tools import probe_root_cpu_context as context
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.probe_dgx_step_separation import save


def validate_point(point, job):
    lineage.validate_point(point, job)
    value = point["root_context"]
    context.validate(value)
    if (value["boot_before"] != point["lineage"]["boot_before"]
            or value["ticks"] != point["ticks"]
            or value["host_before"]["started_ns"] < point["host"][1]["finished_monotonic_ns"]):
        raise ValueError("Supplementary context boot, tick rate or read order differs")


def read_point(pid, membership, job, failure_path=None):
    point = lineage.read_point(pid, membership, job, failure_path)
    try:
        point["root_context"] = context.snapshot()
        validate_point(point, job)
    except (OSError, ValueError, KeyError) as error:
        if failure_path is not None:
            failure = dict(status="invalid_supplementary_root_context", preceding_point=point,
                           error=str(error), error_type=type(error).__name__, scientific_timings_admitted=False)
            if isinstance(error, context.RootContextError):
                failure["partial_context"] = error.evidence
            save(Path(failure_path).with_name("failed_root_context.json"), failure)
        raise
    return point


def evaluate(points, job):
    if len(points) < 2:
        raise ValueError("Require at least two supplementary observations")
    for point in points:
        validate_point(point, job)
    for a, b in zip(points, points[1:]):
        if a["root_context"]["host_after"]["finished_ns"] >= b["host"][0]["started_monotonic_ns"]:
            raise ValueError("Supplementary context overlaps next native observation")
    return dict(schema="supplementary_root_context_v1", observations=len(points),
        intervals=[context.compare(a["root_context"], b["root_context"]) for a, b in zip(points, points[1:])],
        observation_window=context.compare(points[0]["root_context"], points[-1]["root_context"]),
        scientific_timings_admitted=False, environmental_validity_established=False,
        limitations=["Supplementary reads follow each completed native observation; windows are distinct.",
            "Original lineage screens and native timer are unchanged; added observer work can still affect workload runtime.",
            "No CPU attribution, overhead correction, scientific timing admission or exhaustive root-child accounting."])


def lineage_identity(directory):
    value = record(Path(directory).absolute() / "lineage_report.json")
    return dict(value, path="lineage_report.json")


def measure(command, directory, job_id, cpus=20, memory_bytes=96 * 1024**3, timeout_s=60, interval_s=1.):
    measured = lineage.measure(command, directory, job_id, cpus, memory_bytes, timeout_s, interval_s,
                               point_reader=read_point)
    result = dict(status="native_root_context_measured", job_id=job_id,
                  native_wall_s=measured["native_wall_s"], context=evaluate(measured["points"], job_id),
                  lineage_report=lineage_identity(directory),
                  scientific_timings_admitted=False, environmental_validity_established=False)
    save(Path(directory).absolute() / "root_context_report.json", result)
    return measured, result


def measure_native_run(command, directory, job_id, cpus, memory_bytes, timeout_s, interval_s,
                       monitor_host=True, host_interval_s=30.):
    """Keep the native workflow's dictionary contract and retain both raw reports."""
    if (type(cpus) is not int or cpus != 20 or type(memory_bytes) is not int
            or memory_bytes != 96 * 1024**3 or timeout_s != 900 or interval_s != 1.
            or type(interval_s) is bool or monitor_host is not True or host_interval_s != 30.):
        raise ValueError("Require frozen native root-context collection settings")
    measured, _ = measure(command, directory, job_id, cpus, memory_bytes, timeout_s, interval_s)
    return measured
