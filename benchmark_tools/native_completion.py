"""Check retained cgroup thread inventories at native command boundaries."""


def completion_evidence(before, after, anchor_pid, command_finished_ns):
    errors = []
    if type(anchor_pid) is not int or anchor_pid <= 0:
        raise ValueError("Require a positive anchor PID")
    for label, observation in (("before", before), ("after", after)):
        if (observation.get("status") != "observed_within_affinity"
                or observation.get("errors") != [] or observation.get("violating_tids") != []):
            errors.append(f"{label}: incomplete or invalid thread observation")
        if (observation.get("initial_tids") != [anchor_pid]
                or observation.get("final_tids") != [anchor_pid]):
            errors.append(f"{label}: native user subtree is not anchor-only")
        threads = observation.get("threads", [])
        if (len(threads) != 1 or threads[0].get("tid") != anchor_pid
                or type(threads[0].get("start_ticks")) is not int):
            errors.append(f"{label}: missing unique anchor identity")
    if not errors and before["threads"][0]["start_ticks"] != after["threads"][0]["start_ticks"]:
        errors.append("Anchor identity changed")
    if not before.get("scope") or before.get("scope") != after.get("scope"):
        errors.append("Native user subtree changed or missing")
    times = (before.get("started_ns"), before.get("finished_ns"), command_finished_ns,
             after.get("started_ns"), after.get("finished_ns"))
    if (any(type(value) is not int or value < 0 for value in times)
            or not times[0] <= times[1] < times[2] < times[3] <= times[4]):
        errors.append("Thread observations do not bracket command completion")
    return dict(schema="native_completion_v1", anchor_pid=anchor_pid,
        command_finished_ns=command_finished_ns, errors=errors,
        status="anchor_only_at_boundaries" if not errors else "native_completion_unverified",
        before=before, after=after, scientific_timings_admitted=False,
        limitations=["Non-atomic boundary observations, not continuous containment or whole-job teardown accounting.",
                     "Does not exclude work migrated outside the observed subtree or prove descendants were reaped."])
