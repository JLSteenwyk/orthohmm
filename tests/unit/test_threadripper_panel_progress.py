from copy import deepcopy

import pytest

from benchmark_tools.prepare_scaling_inputs import planned_runs
from benchmark_tools.threadripper_panel_progress import overhead_position, position
from benchmark_tools.prepare_threadripper_overhead import arm_order


def attempt(index=0, outcome="exited_zero"):
    return dict(index=index, job_id=100+index, scheduler_state="COMPLETED",
                native_outcome=outcome, review={k: "passed" for k in
                ("runtime", "environment", "resources", "outputs_or_failure")})


def test_empty_panel_requires_preflight():
    r = position(planned_runs(), [])
    assert r["index"] == 0
    assert r["status"] == "next_identity_requires_preflight"
    assert r["scientific_execution_authorized"] is False


@pytest.mark.parametrize("state", ["PENDING", "RUNNING", "COMPLETING", "SUSPENDED"])
def test_live_job_waits_even_with_native_exit(state):
    a = attempt()
    a.update(scheduler_state=state, review=None)
    assert position(planned_runs(), [a])["status"] == "wait_for_existing_job"


def test_terminal_without_review_stops():
    a = attempt()
    a["review"] = None
    assert position(planned_runs(), [a])["status"] == "terminal_attempt_requires_review"


@pytest.mark.parametrize("outcome", ["exited_zero", "exited_nonzero", "timed_out"])
def test_reviewed_native_outcome_advances_without_retry(outcome):
    r = position(planned_runs(), [attempt(outcome=outcome)])
    assert r["index"] == 1
    assert r["reviewed_attempts"][0]["native_outcome"] == outcome
    assert r["automatic_retry"] is False


@pytest.mark.parametrize("state", ["FAILED", "TIMEOUT", "CANCELLED", "NODE_FAIL", "OUT_OF_MEMORY"])
def test_scheduler_failure_never_advances(state):
    a = attempt()
    a["scheduler_state"] = state
    assert position(planned_runs(), [a])["status"] == "infrastructure_failure_requires_resolution"


@pytest.mark.parametrize("field", ["runtime", "environment", "resources", "outputs_or_failure"])
@pytest.mark.parametrize("decision", ["failed", "unresolved"])
def test_adverse_review_stops(field, decision):
    a = attempt()
    a["review"][field] = decision
    assert position(planned_runs(), [a])["status"] == "post_run_review_prevents_continuation"


def test_cannot_skip_past_ineligible_attempt():
    a = attempt()
    a["review"] = None
    with pytest.raises(ValueError, match="predecessor"):
        position(planned_runs(), [a, attempt(1)])


@pytest.mark.parametrize("change", ["repeat", "skip", "duplicate_job", "bool_index", "unknown_state"])
def test_malformed_attempts_rejected(change):
    rows = [attempt(), attempt(1)]
    if change == "repeat": rows[1]["index"] = 0
    if change == "skip": rows[1]["index"] = 2
    if change == "duplicate_job": rows[1]["job_id"] = rows[0]["job_id"]
    if change == "bool_index": rows[0]["index"] = False
    if change == "unknown_state": rows[1]["scheduler_state"] = "MISSING"
    with pytest.raises(ValueError):
        position(planned_runs(), rows)


def test_all_reviewed_retains_failed_native_outcomes():
    r = position(planned_runs(), [attempt(i, "exited_nonzero") for i in range(27)])
    assert r["status"] == "all_attempts_reviewed"
    assert r["index"] is None
    assert len(r["reviewed_attempts"]) == 27
    assert r["scientific_timings_admitted"] is False


def test_frozen_order_and_types_required():
    rows = deepcopy(planned_runs())
    rows[0]["index"] = False
    with pytest.raises(ValueError, match="frozen"):
        position(rows, [])
    rows = planned_runs()[::-1]
    with pytest.raises(ValueError, match="frozen"):
        position(rows, [])


def overhead_tasks():
    tasks = []
    for row in planned_runs():
        for arm in arm_order(row["method"], row["proteomes"], row["repeat"]):
            tasks.append(dict(row, index=len(tasks), pair=row["index"], arm=arm))
    return tasks


def test_overhead_has_54_separate_identities_and_no_admission():
    tasks = overhead_tasks()
    assert overhead_position(tasks, [])['index'] == 0
    result = overhead_position(tasks, [attempt(i) for i in range(54)])
    assert result['status'] == 'all_attempts_reviewed'
    assert not result['scientific_timings_admitted']
    with pytest.raises(ValueError): position(tasks, [])
    with pytest.raises((KeyError, ValueError)): overhead_position(planned_runs(), [])


@pytest.mark.parametrize('outcome', ['exited_nonzero', 'timed_out'])
def test_overhead_stops_even_after_passed_failure_reviews(outcome):
    failed = attempt(outcome=outcome)
    assert overhead_position(overhead_tasks(), [failed])['status'] == 'native_failure_prevents_continuation'
    with pytest.raises(ValueError, match='predecessor'):
        overhead_position(overhead_tasks(), [failed, attempt(1)])


@pytest.mark.parametrize('field,value', [('pair', True), ('arm', 'other'), ('index', False)])
def test_overhead_identity_drift_rejected(field, value):
    tasks = overhead_tasks()
    tasks[0][field] = value
    with pytest.raises(ValueError, match='frozen'): overhead_position(tasks, [])
