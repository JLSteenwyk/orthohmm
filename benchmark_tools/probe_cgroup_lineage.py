"""Prospective aggregate CPU observations independent of sibling cgroup names."""

from pathlib import Path

from benchmark_tools.probe_cgroup_frontier import read_counter, validate_scope
from benchmark_tools.probe_dgx_cpu_hierarchy import usage


def scopes(target):
    path = validate_scope(target)
    return [str(p) for p in reversed(path.parents)] + [str(path)]


def identities(root, names):
    result = {}
    for name in names:
        directory = root / name.lstrip("/")
        if directory.is_symlink() or not directory.is_dir():
            raise ValueError("Missing or symlinked lineage scope")
        stat = directory.stat()
        result[name] = [stat.st_dev, stat.st_ino]
    return result


class LineageSnapshotError(ValueError):
    def __init__(self, error, evidence):
        super().__init__(str(error))
        self.evidence = dict(evidence, status="invalid_lineage_snapshot",
                             error_type=type(error).__name__, error=str(error),
                             scientific_timings_admitted=False)


def snapshot(root, target, boot_path=Path("/proc/sys/kernel/random/boot_id"), *, after_read=None):
    point = dict(status="aggregate_lineage_snapshot", target=target, rows=[],
                 identities_before=None, identities_after=None)
    try:
        names = scopes(target)
        point["boot_before"] = boot_path.read_text().strip()
        point["identities_before"] = identities(root, names)
        for name in names:
            point["rows"].append(read_counter(root, name))
            if after_read is not None:
                after_read(name)
        point["identities_after"] = identities(root, names)
        point["boot_after"] = boot_path.read_text().strip()
        validate(point)
    except (OSError, ValueError, KeyError) as error:
        raise LineageSnapshotError(error, point) from error
    return point


def validate(point):
    names = scopes(point["target"])
    if point["status"] != "aggregate_lineage_snapshot":
        raise ValueError("Invalid lineage snapshot status")
    if not point["boot_before"] or point["boot_before"] != point["boot_after"]:
        raise ValueError("Missing or changed boot identity")
    before, after = point["identities_before"], point["identities_after"]
    if set(before) != set(names) or before != after:
        raise ValueError("Lineage identity changed or coverage is incomplete")
    for identity in before.values():
        if (not isinstance(identity, list) or len(identity) != 2
                or any(type(v) is not int or v < 0 for v in identity)):
            raise ValueError("Invalid filesystem identity")
    if len({tuple(v) for v in before.values()}) != len(names):
        raise ValueError("Lineage scopes alias the same identity")
    if [row["scope"] for row in point["rows"]] != names:
        raise ValueError("Missing, duplicate or unordered lineage counters")
    previous = -1
    for row in point["rows"]:
        start, end = row["started_ns"], row["finished_ns"]
        if (type(start) is not int or type(end) is not int
                or not previous <= start <= end):
            raise ValueError("Invalid lineage read order")
        usage(row["raw"])
        previous = end


def compare(left, right):
    for point in (left, right):
        validate(point)
    if (left["target"] != right["target"] or left["boot_before"] != right["boot_before"]
            or left["identities_before"] != right["identities_before"]):
        raise ValueError("Boot, target or lineage identity changed between observations")
    if left["rows"][-1]["finished_ns"] >= right["rows"][0]["started_ns"]:
        raise ValueError("Lineage observation windows overlap")
    deltas = []
    for a, b in zip(left["rows"], right["rows"]):
        delta = usage(b["raw"]) - usage(a["raw"])
        if delta < 0:
            raise ValueError("Lineage CPU counter decreased")
        deltas.append(dict(scope=a["scope"], cpu_usec=delta,
                           read_midpoint_elapsed_s=((b["started_ns"] + b["finished_ns"])
                             - (a["started_ns"] + a["finished_ns"])) / 2e9))
    # Adjacent differences telescope; overlapping ancestor totals must not be summed.
    complements = [dict(ancestor=a["scope"], excluded_child=b["scope"],
                        signed_cpu_usec=a["cpu_usec"] - b["cpu_usec"])
                   for a, b in zip(deltas, deltas[1:])]
    return dict(status="aggregate_lineage_cpu_diagnostic", aggregate_deltas=deltas,
                signed_complements=complements,
                root_minus_target_cpu_usec=deltas[0]["cpu_usec"] - deltas[-1]["cpu_usec"],
                scientific_timings_admitted=False, controlled_workload_verified=False,
                limitations=[
                    "Aggregate ancestor counters are inclusive and overlapping; never sum them.",
                    "Signed adjacent differences have non-atomic read windows and accounting delay; they are not bounds.",
                    "No sibling inventory is needed, but no individual transient service is identified.",
                    "Counter stability cannot prove process membership stability or absence of non-CPU interference.",
                    "Kernel transient-descendant accounting, collector overhead and scientific inclusion remain unvalidated.",
                    "No original flags are removed, thresholds changed, or native wall times corrected."])
