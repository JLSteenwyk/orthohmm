"""Determine panel position from recorded attempts; never submit or authorize work."""

from benchmark_tools.prepare_scaling_inputs import planned_runs

LIVE = {"PENDING", "RUNNING", "CONFIGURING", "COMPLETING", "SUSPENDED", "RESIZING"}
TERMINAL = {"COMPLETED", "FAILED", "TIMEOUT", "CANCELLED", "OUT_OF_MEMORY",
            "NODE_FAIL", "PREEMPTED", "BOOT_FAIL", "DEADLINE", "REVOKED"}
NATIVE_OUTCOMES = {"exited_zero", "exited_nonzero", "timed_out"}


def position(runs, attempts):
    """Require a contiguous one-attempt prefix and explicit post-run review.

    The caller must independently validate scheduler observations and review
    evidence. This function checks their state-machine consistency only.
    """
    expected = planned_runs()
    return _position(runs, attempts, expected)


def overhead_position(tasks, attempts):
    """Keep the 54 engineering identities separate and stop on native failure."""
    from benchmark_tools.prepare_threadripper_overhead import arm_order

    expected = []
    for run in planned_runs():
        for arm in arm_order(run["method"], run["proteomes"], run["repeat"]):
            expected.append(dict(index=len(expected), pair=run["index"], arm=arm,
                method=run["method"], proteomes=run["proteomes"], repeat=run["repeat"]))
    return _position(tasks, attempts, expected, stop_on_native_failure=True)


def _position(runs, attempts, expected, *, stop_on_native_failure=False):
    identities = [{k: r[k] for k in expected[0]} for r in runs]
    if identities != expected or any(
            type(r[k]) is not type(e[k]) for r, e in zip(identities, expected) for k in e):
        raise ValueError("Require frozen panel identities and order")
    if not isinstance(attempts, list) or len(attempts) > len(runs):
        raise ValueError("Require one-attempt panel prefix")
    jobs = set()
    reviewed = []
    for index, attempt in enumerate(attempts):
        if type(attempt.get("index")) is not int or attempt["index"] != index:
            raise ValueError("Attempts skip, reorder or repeat a panel identity")
        job = attempt.get("job_id")
        if type(job) is not int or job <= 0 or job in jobs:
            raise ValueError("Require unique positive scheduler job IDs")
        jobs.add(job)
        state = attempt.get("scheduler_state")
        if state not in LIVE | TERMINAL:
            raise ValueError("Unknown scheduler state; obtain authoritative status")
        review = attempt.get("review")
        outcome = attempt.get("native_outcome")
        if outcome is not None and outcome not in NATIVE_OUTCOMES:
            raise ValueError("Unknown native outcome")
        status = None
        if state in LIVE:
            if review is not None:
                raise ValueError("Cannot admit review of a live scheduler job")
            status = "wait_for_existing_job"
        elif review is None:
            status = "terminal_attempt_requires_review"
        else:
            if not isinstance(review, dict) or set(review) != {
                    "runtime", "environment", "resources", "outputs_or_failure"}:
                raise ValueError("Require complete post-run review fields")
            if any(v not in {"passed", "failed", "unresolved"} for v in review.values()):
                raise ValueError("Unknown post-run review decision")
            if state != "COMPLETED":
                status = "infrastructure_failure_requires_resolution"
            elif not all(v == "passed" for v in review.values()):
                status = "post_run_review_prevents_continuation"
            elif outcome is None:
                raise ValueError("Reviewed attempt lacks native outcome")
            elif stop_on_native_failure and outcome != "exited_zero":
                status = "native_failure_prevents_continuation"
            else:
                reviewed.append(dict(index=index, job_id=job, native_outcome=outcome))
        if status:
            if index != len(attempts) - 1:
                raise ValueError("Later attempt follows an unreviewed or ineligible predecessor")
            return result(status, index, reviewed, job)
    if len(attempts) == len(runs):
        return result("all_attempts_reviewed", None, reviewed)
    return result("next_identity_requires_preflight", len(attempts), reviewed)


def result(status, index, reviewed, job=None):
    return dict(status=status, index=index, existing_job_id=job,
                reviewed_attempts=reviewed, automatic_retry=False,
                scientific_execution_authorized=False, scientific_timings_admitted=False,
                limitations=["Recorded state consistency only, not live scheduler or evidence verification.",
                             "Next identity still requires pinned recipe, allocation and fresh quiet-host preflight.",
                             "All attempts reviewed does not imply all native runs succeeded or publication readiness."])
