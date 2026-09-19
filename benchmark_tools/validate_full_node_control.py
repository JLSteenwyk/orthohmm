"""Independent workload-witness checks for the prospective full-node controls."""

import math


def require(condition, message):
    if not condition:
        raise ValueError(message)


def integer(value):
    return type(value) is int and value > 0


def nonnegative(value):
    return type(value) in (int, float) and math.isfinite(value) and value >= 0


def validate_work(row):
    require(integer(row["pid"]) and integer(row["parent_pid"]), "Invalid process identity")
    require(integer(row["started_ns"]) and integer(row["finished_ns"])
            and row["finished_ns"] - row["started_ns"] >= 20_000_000_000,
            "Work shorter than frozen duration")
    for name in ("self_cpu_s", "waited_child_user_s", "waited_child_system_s"):
        require(nonnegative(row[name]), "Invalid CPU witness")
    require(row["self_cpu_s"] > 0, "No measured worker CPU")
    require(row["affinity"] == row["final_affinity"]
            and row["membership"] == row["final_membership"], "Worker scope changed")
    require(type(row["creations"]) is int and row["creations"] >= 0, "Invalid creation count")
    require(row["creation_cap_reached"] is False, "Creation cap reached")


def validate_witnesses(mode, ready, done, native, native_membership,
                       observer_membership, competitor_ready=None, competitor=None):
    require(mode in {"steady", "churn", "contended"}, "Unknown control condition")
    work_mode = "churn" if mode == "churn" else "steady"
    cpus = done["cpus"]
    require(len(cpus) == 20 and all(type(cpu) is int and cpu >= 0 for cpu in cpus)
            and cpus == sorted(set(cpus)), "Require 20 distinct ordered CPUs")
    require(ready["cpus"] == cpus and ready["mode"] == done["mode"] == work_mode,
            "Workload design differs")
    require(done["duration_s"] == 20 and done["creation_cap"] == 200000,
            "Workload bounds differ")
    require(integer(done["pid"]) and ready["pid"] == done["pid"], "Parent identity differs")
    require(ready["membership"] == done["membership"] == native_membership
            and observer_membership != native_membership, "Native/observer scope differs")
    require(type(native["exit_code"]) is int and native["exit_code"] == 0
            and native["timed_out"] is False, "Native command failed")
    require(integer(native["started_ns"]) and integer(native["finished_ns"])
            and native["started_ns"] <= done["started_ns"] < done["finished_ns"] <= native["finished_ns"],
            "Parent work outside command boundaries")
    require(len(ready["workers"]) == len(done["workers"]) == 20, "Missing workers")
    pids = []
    for cpu, before, row in zip(cpus, ready["workers"], done["workers"]):
        validate_work(row)
        require(set(before) == {"pid", "parent_pid", "cpu", "allowed", "affinity", "membership"},
                "Incomplete worker readiness")
        require(all(row[key] == value for key, value in before.items()), "Readiness identity differs")
        require(row["parent_pid"] == done["pid"] and row["pid"] != done["pid"], "Wrong parent")
        require(row["cpu"] == cpu and row["affinity"] == [cpu] and row["allowed"] == cpus,
                "Worker affinity differs")
        require(row["membership"] == native_membership, "Wrong native membership")
        require(row["mode"] == work_mode and row["duration_s"] == 20
                and row["creation_cap"] == 200000, "Worker design differs")
        require(done["started_ns"] <= row["started_ns"] < row["finished_ns"] <= done["finished_ns"],
                "Worker outside parent boundaries")
        require((0 < row["creations"] < 200000) if work_mode == "churn" else row["creations"] == 0,
                "Creation count disagrees with condition")
        pids.append(row["pid"])
    require(len(set(pids)) == 20, "Duplicate worker PID")
    require(set(done["statuses"]) == {str(pid) for pid in pids}
            and all(type(status) is int and status == 0 for status in done["statuses"].values()),
            "Missing or failed worker exit status")
    start = max(row["started_ns"] for row in done["workers"])
    finish = min(row["finished_ns"] for row in done["workers"])
    require(finish - start >= 19_500_000_000, "Insufficient common worker overlap")
    overlap = None
    if mode == "contended":
        require(competitor_ready is not None and competitor is not None, "Missing competitor")
        validate_work(competitor)
        require(set(competitor_ready) == {"pid", "parent_pid", "affinity", "membership"},
                "Incomplete competitor readiness")
        require(all(competitor[key] == value for key, value in competitor_ready.items()),
                "Competitor readiness differs")
        require(competitor["membership"] == observer_membership
                and competitor["affinity"] == [min(cpus)], "Wrong competitor scope/affinity")
        require(competitor["pid"] not in pids + [done["pid"]], "Competitor PID reused")
        require(competitor["creations"] == 0 and 5 <= competitor["self_cpu_s"] <= 21,
                "Invalid competitor CPU dose")
        start, finish = max(start, competitor["started_ns"]), min(finish, competitor["finished_ns"])
        overlap = (finish - start) / 1e9
        require(overlap >= 19.5, "Insufficient competitor overlap")
    else:
        require(competitor_ready is None and competitor is None, "Unexpected competitor")
    return dict(status="workload_witnesses_validated", common_started_ns=start,
                common_finished_ns=finish, competitor_overlap_s=overlap,
                scientific_timings_admitted=False)


def common_intervals(points, validation):
    """Use the entire enclosing original host window, not a midpoint shortcut."""
    start, finish = validation["common_started_ns"], validation["common_finished_ns"]
    return [i for i, (a, b) in enumerate(zip(points, points[1:]))
            if start <= a["host"][0]["started_monotonic_ns"]
            and b["host"][1]["finished_monotonic_ns"] <= finish]
