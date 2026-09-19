"""Read-only root/slice CPU context; signed differences are not causal bounds."""

import os
from pathlib import Path
import time

from benchmark_tools.probe_cgroup_lineage import identities
from benchmark_tools.probe_cgroup_frontier import read_counter
from benchmark_tools.probe_dgx_cpu_hierarchy import usage
from benchmark_tools.probe_host_counters import parse_cpu

SCOPES = ("/", "/system.slice", "/user.slice", "/init.scope")


def read_text(path):
    started = time.monotonic_ns()
    raw = path.read_text()
    return dict(started_ns=started, finished_ns=time.monotonic_ns(), raw=raw)


def pids(raw):
    tokens = raw.split()
    if any(not t.isascii() or not t.isdigit() or int(t) <= 0 for t in tokens):
        raise ValueError("Invalid root membership PID")
    return sorted(set(map(int, tokens)))


class RootContextError(ValueError):
    def __init__(self, error, evidence):
        super().__init__(str(error))
        self.evidence = dict(evidence, status="invalid_root_cpu_context", error=str(error),
                             error_type=type(error).__name__, scientific_timings_admitted=False)


def snapshot(root=Path("/sys/fs/cgroup"), proc=Path("/proc")):
    point = dict(status="root_cpu_context_v1", rows=[])
    try:
        point["ticks"] = os.sysconf("SC_CLK_TCK")
        point["boot_before"] = (proc / "sys/kernel/random/boot_id").read_text().strip()
        point["identities_before"] = identities(root, SCOPES)
        point["host_before"] = read_text(proc / "stat")
        point["members_before"] = read_text(root / "cgroup.procs")
        for scope in SCOPES:
            point["rows"].append(read_counter(root, scope))
        point["members_after"] = read_text(root / "cgroup.procs")
        point["host_after"] = read_text(proc / "stat")
        point["identities_after"] = identities(root, SCOPES)
        point["boot_after"] = (proc / "sys/kernel/random/boot_id").read_text().strip()
        validate(point)
    except (OSError, ValueError, KeyError) as error:
        raise RootContextError(error, point) from error
    return point


def validate(point):
    if point["status"] != "root_cpu_context_v1" or not point["boot_before"]:
        raise ValueError("Invalid root context identity")
    if point["boot_before"] != point["boot_after"]:
        raise ValueError("Boot changed during context read")
    if type(point["ticks"]) is not int or point["ticks"] <= 0:
        raise ValueError("Invalid host clock tick frequency")
    before, after = point["identities_before"], point["identities_after"]
    if set(before) != set(SCOPES) or before != after:
        raise ValueError("Missing or changed scope identity")
    if (any(not isinstance(v, list) or len(v) != 2 or any(type(n) is not int or n < 0 for n in v)
            for v in before.values()) or len({tuple(v) for v in before.values()}) != len(SCOPES)):
        raise ValueError("Invalid or aliased scope identity")
    if [r["scope"] for r in point["rows"]] != list(SCOPES):
        raise ValueError("Missing, duplicate or unordered scope counters")
    sequence = [point["host_before"], point["members_before"], *point["rows"],
                point["members_after"], point["host_after"]]
    previous = -1
    for item in sequence:
        start, end = item["started_ns"], item["finished_ns"]
        if type(start) is not int or type(end) is not int or not previous <= start <= end:
            raise ValueError("Invalid context read order")
        previous = end
    for row in point["rows"]:
        usage(row["raw"])
    for key in ("members_before", "members_after"):
        pids(point[key]["raw"])
    a, b = (parse_cpu(point[key]["raw"]) for key in ("host_before", "host_after"))
    if any(b[key] < a[key] for key in a):
        raise ValueError("Host CPU counter decreased during context read")


def compare(left, right):
    for point in (left, right):
        validate(point)
    if (left["boot_before"] != right["boot_before"] or left["ticks"] != right["ticks"]
            or left["identities_before"] != right["identities_before"]
            or left["host_after"]["finished_ns"] >= right["host_before"]["started_ns"]):
        raise ValueError("Changed context identity or overlapping observations")
    deltas = {a["scope"]: usage(b["raw"]) - usage(a["raw"])
              for a, b in zip(left["rows"], right["rows"])}
    if min(deltas.values()) < 0:
        raise ValueError("Scope CPU counter decreased")
    a, b = (parse_cpu(p[key]["raw"]) for p, key in ((left, "host_before"), (right, "host_after")))
    host = {key: b[key] - a[key] for key in a}
    if min(host.values()) < 0:
        raise ValueError("Host CPU counter decreased")
    memberships = [pids(p[key]["raw"]) for p in (left, right) for key in ("members_before", "members_after")]
    return dict(status="root_cpu_context_comparison", scope_cpu_usec=deltas,
        root_minus_system_cpu_usec=deltas["/"] - deltas["/system.slice"],
        root_minus_three_named_children_cpu_usec=deltas["/"] - sum(deltas[s] for s in SCOPES[1:]),
        enclosing_host_category_ticks=host, ticks=left["ticks"],
        root_membership_snapshots=memberships,
        observed_root_membership_changed=any(m != memberships[0] for m in memberships[1:]),
        scientific_timings_admitted=False, environmental_validity_established=False,
        limitations=["Root and non-root counters have distinct accounting semantics and non-atomic read windows.",
            "Three named children are not an exhaustive root-child inventory; residual is not root-task CPU.",
            "PID membership does not measure task CPU, identify PID reuse, or exclude unobserved transient tasks.",
            "Host categories have an enclosing window, not exact scope-counter boundaries; guest fields overlap user/nice.",
            "No flags removed, thresholds changed, CPU attributed or scientific timings admitted."])
