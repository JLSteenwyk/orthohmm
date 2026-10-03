"""Compare observed process identities with an explicitly reviewed background inventory.

This is one prerequisite for environmental review, never a quiet-host certificate.
It does not create policy entries from the observations it is asked to assess.
"""

import math
from pathlib import PurePosixPath

from benchmark_tools.observe_host_competition import analyze

SHARED_SCOPE = "shared_host_matched_resources"


def shared_environment(policy):
    if policy.get("schema") != "threadripper_environment_policy_v3":
        return False
    if (policy.get("execution_scope") != SHARED_SCOPE
            or policy.get("foreign_cpu_role") != "diagnostic_only"
            or policy.get("native_pressure_role") != "diagnostic_only"
            or policy.get("preflight_pressure_role") != "diagnostic_only"):
        raise ValueError("Require explicit shared-host diagnostic contention roles")
    return True


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


def kernel_type(row, sample, boot):
    if (sample.get("schema") != "threadripper_typed_process_snapshot_v1"
            or sample.get("boot_id") != boot or "kernel_identity_error" in row):
        raise ValueError("Require same-boot typed process observations")
    proof = row.get("kernel_identity", {})
    if (type(proof.get("pid")) is not int or proof["pid"] != row["pid"]
            or type(proof.get("tgid")) is not int or proof["tgid"] != row["pid"]
            or type(proof.get("kthread")) is not int or proof["kthread"] not in (0, 1)):
        raise ValueError("Missing or inconsistent kernel identity evidence")
    start = number(proof["started_monotonic_s"])
    finish = number(proof["finished_monotonic_s"])
    if not row["observed_monotonic_s"] <= start <= finish <= sample["finished_monotonic_s"]:
        raise ValueError("Kernel identity observation outside snapshot bounds")
    return proof["kthread"]


def review(policy, before, after, *, boot_id, job_scope, observer_pid):
    """Classify two snapshots against same-boot, exact-identity policy entries.

    The caller must separately verify policy/evidence hashes and the review's
    factual basis. The observer must be inside the job, not exempted by PID.
    """
    if policy.get("schema") == "threadripper_process_policy_v3":
        return review_shared(policy, before, after, boot_id=boot_id,
                             job_scope=job_scope, observer_pid=observer_pid)
    typed = policy.get("schema") == "threadripper_process_policy_v2"
    if (policy.get("schema") not in {"threadripper_process_policy_v1", "threadripper_process_policy_v2"}
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
    allowed, kernel_pids = {}, set()
    for row in policy["ordinary_processes"]:
        key = identity(row)
        if key[0] in allowed or PurePosixPath(key[2]).is_relative_to(scope):
            raise ValueError("Duplicate or in-job background policy entry")
        classification = row.get("classification")
        classes = {"ordinary_background", "reviewed_kernel_thread"} if typed else {"ordinary_background"}
        if (classification not in classes
                or not isinstance(row.get("reason"), str) or not row["reason"].strip()):
            raise ValueError("Require explicit ordinary-background classification and reason")
        allowed[key[0]] = key
        if classification == "reviewed_kernel_thread":
            kernel_pids.add(key[0])
    a, b = inventory(before), inventory(after)
    types = {}
    if typed:
        for label, sample, rows in (("before", before, a), ("after", after, b)):
            for pid, row in rows.items():
                types[label, pid] = kernel_type(row, sample, boot_id)
    def matches(label, pid, row):
        key = identity(row)
        expected = allowed.get(pid)
        if typed and pid in allowed and types[label, pid] != int(pid in kernel_pids):
            return False
        return key[:3] == expected[:3] if pid in kernel_pids else key == expected
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
            if not group(row["cgroup"]).is_relative_to(scope) and not matches(label, pid, row):
                unreviewed.append(dict(snapshot=label, pid=pid, created=key[1],
                                       cgroup=key[2], name=key[3]))
        missing.extend(dict(snapshot=label, pid=pid) for pid in sorted(allowed)
                       if pid not in rows or not matches(label, pid, rows[pid]))
    kernel_name_changes = []
    in_job_changes = []
    for pid in sorted(set(a) & set(b)):
        # Do not ignore a process migrating into the job between observations.
        name_only_kernel = (pid in kernel_pids and matches("before", pid, a[pid])
                            and matches("after", pid, b[pid]))
        within_job = (typed and pid != observer_pid
                      and group(a[pid]["cgroup"]).is_relative_to(scope)
                      and group(b[pid]["cgroup"]).is_relative_to(scope))
        changed_identity = identity(a[pid]) != identity(b[pid])
        if within_job and changed_identity:
            in_job_changes.append(dict(pid=pid, before=list(identity(a[pid])),
                                       after=list(identity(b[pid]))))
        if changed_identity and not name_only_kernel and not within_job:
            changed.append(dict(pid=pid, reason="identity_or_membership_changed"))
        if within_job and b[pid]["created"] < a[pid]["created"]:
            changed.append(dict(pid=pid, reason="process_creation_time_decreased"))
        if name_only_kernel and a[pid]["name"] != b[pid]["name"]:
            kernel_name_changes.append(dict(pid=pid, before=a[pid]["name"], after=b[pid]["name"]))
        replaced_within_job = within_job and b[pid]["created"] > a[pid]["created"]
        if not replaced_within_job and any(b[pid][k] < a[pid][k] for k in ("user_s", "system_s")):
            changed.append(dict(pid=pid, reason="cpu_counter_decreased"))
    diagnostic = analyze(before, after, str(scope), observer_pid)
    unresolved = bool(unreviewed or changed or missing or diagnostic["sampling_error_count"]
                      or diagnostic["unmatched_foreign_processes"] or diagnostic["uncertain_processes"])
    result = dict(schema="threadripper_process_policy_review_v2" if typed else "threadripper_process_policy_review_v1",
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
    if typed:
        result["observed_in_job_identity_changes"] = in_job_changes
        result["reviewed_kernel_name_changes"] = kernel_name_changes
        result["limitations"].append("Verified kernel type permits reviewed name changes only; CPU and other checks remain.")
        result["limitations"].append("Non-observer transitions observed wholly inside the job are retained, not foreign competition; movement between samples remains unobserved.")
    return result


def review_shared(policy, before, after, *, boot_id, job_scope, observer_pid):
    """Validate observation/attribution while retaining ordinary outside churn."""
    if (policy.get("schema") != "threadripper_process_policy_v3"
            or policy.get("execution_scope") != SHARED_SCOPE
            or not isinstance(boot_id, str) or not boot_id.strip()
            or policy.get("boot_id") != boot_id
            or not isinstance(policy.get("review_reference"), str)
            or not policy["review_reference"].strip()):
        raise ValueError("Require explicit same-boot shared-host process policy")
    scope = group(job_scope)
    if scope == PurePosixPath("/") or type(observer_pid) is not int or observer_pid <= 0:
        raise ValueError("Require scoped job and observer identity")
    a, b = inventory(before), inventory(after)
    if before["finished_monotonic_s"] >= after["started_monotonic_s"]:
        raise ValueError("Require distinct ordered snapshots")
    for sample, rows in ((before, a), (after, b)):
        if (sample.get("boot_id") != boot_id or observer_pid not in rows
                or not group(rows[observer_pid]["cgroup"]).is_relative_to(scope)
                or kernel_type(rows[observer_pid], sample, boot_id) != 0):
            raise ValueError("Observer boot, job membership or kernel identity differs")
    if identity(a[observer_pid]) != identity(b[observer_pid]):
        raise ValueError("Observer identity changed")
    failures = []
    for pid in sorted(set(a) & set(b)):
        first, last = a[pid], b[pid]
        inside = [group(row["cgroup"]).is_relative_to(scope) for row in (first, last)]
        if any(inside) and inside[0] != inside[1]:
            failures.append(dict(pid=pid, reason="job_membership_changed"))
        if all(inside):
            if (last["created"] < first["created"] or
                    last["created"] == first["created"] and
                    any(last[k] < first[k] for k in ("user_s", "system_s"))):
                failures.append(dict(pid=pid, reason="in_job_counter_or_identity_decreased"))
    diagnostic = analyze(before, after, str(scope), observer_pid)
    # Outside collection gaps/churn limit contention estimates; native accounting
    # and the observer's identity are validated separately, never inferred quiet.
    return dict(schema="threadripper_shared_process_review_v1",
        status="shared_host_observation_valid" if not failures else "invalid_shared_host_observation",
        process_policy_matched=not failures, controlled_workload_verified=False,
        execution_scope=SHARED_SCOPE, foreign_cpu_used_for_eligibility=False,
        scientific_timings_admitted=False, boot_id=boot_id, job_scope=str(scope),
        review_reference=policy["review_reference"], changed_processes=failures,
        cpu_diagnostic=diagnostic,
        limitations=["Background activity and churn are annotated, not an isolation certificate.",
            "Outside sampling errors and unmatched processes make CPU estimates incomplete.",
            "Native cgroup accounting, resources and runtime integrity require separate checks.",
            "Contention distortion is unknown and can be method-dependent."])
