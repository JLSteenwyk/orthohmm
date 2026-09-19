"""Prospective paired CPU brackets; never changes original measurement gates."""

from benchmark_tools.measure_native_frontier_step import read_frontier_point, validate_point
from benchmark_tools.measure_native_hierarchy_step import interval_point
from benchmark_tools.probe_host_counters import summarize
from benchmark_tools.probe_interval_cpu import interval
from benchmark_tools.probe_native_pressure import compare as compare_pressure
from benchmark_tools.probe_cgroup_frontier import compare as compare_frontier


def validate(point, job):
    validate_point(point, job)
    if "native_pressure" not in point:
        raise ValueError("Require native pressure alongside CPU brackets")
    earlier, later = point["hierarchy_host_after"], point["host"][1]
    summarize(earlier, later, point["ticks"])
    return interval_point({**point, "host": [point["host"][0], earlier]}, job)


def read_point(pid, membership, job, failure_path=None):
    point = read_frontier_point(pid, membership, job, failure_path, native_pressure=True)
    validate(point, job)
    return point


def compare(left, right, job, *, enforce_gap=True):
    narrow_left, narrow_right = (validate(p, job) for p in (left, right))
    # These two residuals share native counters, timestamps, and the left host read.
    # Frontier/pressure reads extend only the outer right host endpoint.
    outer = interval(interval_point(left, job), interval_point(right, job), job, enforce_gap=enforce_gap)
    narrow = interval(narrow_left, narrow_right, job, enforce_gap=enforce_gap)
    if outer["native_cpu_s"] != narrow["native_cpu_s"] or outer["wall_s"] != narrow["wall_s"]:
        raise ValueError("Paired native measurements differ")
    return dict(status="dual_bracket_cpu_diagnostic", outer=outer, narrow=narrow,
                native_pressure=compare_pressure(left["native_pressure"], right["native_pressure"], job),
                frontier=compare_frontier(left["frontier"], right["frontier"]),
                scientific_timings_admitted=False, controlled_workload_verified=False,
                limitations=["Both original thresholds retained; narrow is prospective diagnostic only.",
                             "Non-atomic readings and accounting delay remain; neither residual identifies foreign CPU.",
                             "No wall-time subtraction, overhead admission or scientific panel authorization."])
