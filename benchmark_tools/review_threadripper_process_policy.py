"""Compare observed process identities with an explicitly reviewed background inventory.

This is one prerequisite for environmental review, never a quiet-host certificate.
It does not create policy entries from the observations it is asked to assess.
"""

import math
from pathlib import PurePosixPath

from benchmark_tools.observe_host_competition import analyze


def number(value, positive=False):
    if (type(value) not in (int, float) or not math.isfinite(value)
            or (value <= 0 if positive else value < 0)):
        raise ValueError("Invalid process counter or timestamp")
    return value


def group(value):
    if not isinstance(value, str):
        raise ValueError("Invalid cgroup path")
    path = PurePosixPath(value)
    if not path.is_absolute() or value.startswith("//") or ".." in path.parts or str(path) != value:
        raise ValueError("Require canonical absolute cgroup path")
    return path


def identity(row):
    if type(row["pid"]) is not int or row["pid"] <= 0:
        raise ValueError("Invalid process PID")
    number(row["created"], positive=True)
    group(row["cgroup"])
    if not isinstance(row["name"], str) or not row["name"].strip():
        raise ValueError("Missing process name")
    return row["pid"], row["created"], row["cgroup"], row["name"]


def inventory(sample):
    start = number(sample["started_monotonic_s"])
    finish = number(sample["finished_monotonic_s"])
    if finish < start or not isinstance(sample["errors"], list):
        raise ValueError("Invalid snapshot bounds or errors")
    rows = sample["processes"]
    if not isinstance(rows, list) or not rows:
        raise ValueError("Require nonempty observed process inventory")
    indexed = {}
    for row in rows:
        key = identity(row)
        if key[0] in indexed:
            raise ValueError("Duplicate process PID")
        for field in ("user_s", "system_s"):
            number(row[field])
        when = number(row["observed_monotonic_s"])
        if not start <= when <= finish:
            raise ValueError("Process observation outside snapshot bounds")
        indexed[key[0]] = row
    return indexed


def review(policy, before, after, *, boot_id, job_scope, observer_pid):
    """Classify two snapshots against same-boot, exact-identity policy entries.

    The caller must separately verify policy/evidence hashes and the review's
    factual basis. The observer must be inside the job, not exempted by PID.
    """
    if (policy.get("schema") != "threadripper_process_policy_v1"
            or not isinstance(boot_id, str) or not boot_id.strip()
            or policy.get("boot_id") != boot_id
            or not isinstance(policy.get("review_reference"), str)
            or not policy["review_reference"].strip()):
        raise ValueError("Require explicit same-boot process policy review")
    scope = group(job_scope)
    if scope == PurePosixPath("/"):
        raise ValueError("Cannot exempt the entire host")
    if type(observer_pid) is not int or observer_pid <= 0:
        raise ValueError("Invalid observer PID")
    if not isinstance(policy.get("ordinary_processes"), list):
        raise ValueError("Require an explicit ordinary-process list")
    allowed = {}
    for row in policy["ordinary_processes"]:
        key = identity(row)
        if key[0] in allowed or PurePosixPath(key[2]).is_relative_to(scope):
            raise ValueError("Duplicate or in-job background policy entry")
        if (row.get("classification") != "ordinary_background"
                or not isinstance(row.get("reason"), str) or not row["reason"].strip()):
            raise ValueError("Require explicit ordinary-background classification and reason")
        allowed[key[0]] = key
    a, b = inventory(before), inventory(after)
    if before["finished_monotonic_s"] >= after["started_monotonic_s"]:
        raise ValueError("Require distinct ordered snapshots")
    for rows in (a, b):
        if (observer_pid not in rows
                or not group(rows[observer_pid]["cgroup"]).is_relative_to(scope)):
            raise ValueError("Observer must be observed inside this job")
    if identity(a[observer_pid]) != identity(b[observer_pid]):
        raise ValueError("Observer identity changed")

    unreviewed, changed, missing = [], [], []
    for label, rows in (("before", a), ("after", b)):
        for pid, row in sorted(rows.items()):
            key = identity(row)
            if not group(row["cgroup"]).is_relative_to(scope) and allowed.get(pid) != key:
                unreviewed.append(dict(snapshot=label, pid=pid, created=key[1],
                                       cgroup=key[2], name=key[3]))
        missing.extend(dict(snapshot=label, pid=pid) for pid in sorted(allowed)
                       if pid not in rows or identity(rows[pid]) != allowed[pid])
    for pid in sorted(set(a) & set(b)):
        # Do not ignore a process migrating into the job between observations.
        if identity(a[pid]) != identity(b[pid]):
            changed.append(dict(pid=pid, reason="identity_or_membership_changed"))
        elif any(b[pid][k] < a[pid][k] for k in ("user_s", "system_s")):
            changed.append(dict(pid=pid, reason="cpu_counter_decreased"))
    diagnostic = analyze(before, after, str(scope), observer_pid)
    unresolved = bool(unreviewed or changed or missing or diagnostic["sampling_error_count"]
                      or diagnostic["unmatched_foreign_processes"] or diagnostic["uncertain_processes"])
    return dict(schema="threadripper_process_policy_review_v1",
        status="unresolved_process_inventory" if unresolved else "reviewed_process_inventory_observed",
        process_policy_matched=not unresolved, controlled_workload_verified=False,
        scientific_timings_admitted=False, boot_id=boot_id, job_scope=str(scope),
        review_reference=policy["review_reference"], unreviewed_processes=unreviewed,
        changed_processes=changed, missing_policy_processes=missing, cpu_diagnostic=diagnostic,
        limitations=["Exact process identities do not attest executables, service configuration or intent.",
                     "The caller must verify frozen policy provenance and the underlying review evidence.",
                     "Two snapshots can miss short-lived work and do not establish quiet I/O or memory bandwidth.",
                     "Policy matching does not authorize native release or replace whole-run observation.",
                     "The CPU diagnostic threshold is not a timing admission rule."])
