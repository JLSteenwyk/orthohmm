import copy

import pytest

from benchmark_tools.admit_component_recovery import check_phase_links, require_terminal_job


def phases():
    candidate = {"original_execution_plan": "plan", "original_producer": "producer",
        "status": "isolated_component_candidates_exactly_reconstructed", "inputs_reverified_unchanged": True,
        "native_inference_rerun": False, "defaults_promoted": False, "goal_completed": False}
    plan = {"original_execution_plan": "plan", "original_failed_producer": "producer", "native_probe": "probe"}
    scored = {"status": "scored_recovery_pending_independent_admission", "native_inference_rerun": False,
              "defaults_promoted": False, "goal_completed": False}
    return candidate, plan, scored


def test_three_phases_must_share_exact_original_artifacts_without_modifying_them():
    items = phases()
    before = copy.deepcopy(items)
    check_phase_links("plan", "producer", "probe", *items)
    assert before == items


@pytest.mark.parametrize("phase,field", [(0, "original_execution_plan"), (0, "original_producer"),
    (1, "original_execution_plan"), (1, "original_failed_producer"), (1, "native_probe"), (0, "status"), (2, "status")])
def test_changed_artifact_links_or_incomplete_phases_cannot_be_admitted(phase, field):
    items = phases()
    items[phase][field] = "wrong"
    with pytest.raises(ValueError):
        check_phase_links("plan", "producer", "probe", *items)


@pytest.mark.parametrize("phase", [0, 2])
@pytest.mark.parametrize("field", ["native_inference_rerun", "defaults_promoted", "goal_completed"])
def test_scope_changes_cannot_be_admitted(phase, field):
    items = phases()
    items[phase][field] = True
    with pytest.raises(ValueError, match="scientific scope"):
        check_phase_links("plan", "producer", "probe", *items)


@pytest.mark.parametrize("field,value", [("State", "RUNNING"), ("State", "FAILED"), ("ExitCode", "0:11"),
    ("AllocCPUS", "1"), ("JobID", "different")])
def test_each_phase_requires_its_own_successful_terminal_allocation(field, value):
    row = {"JobID": "24220", "State": "COMPLETED", "ExitCode": "0:0", "AllocCPUS": "16"}
    row[field] = value
    text = "|".join(row) + "\n" + "|".join(row.values()) + "\n"
    with pytest.raises(ValueError):
        require_terminal_job(text, "24220", 16)


def test_successful_terminal_allocation_is_accepted():
    text = "JobID|State|ExitCode|AllocCPUS\n24220|COMPLETED|0:0|16\n"
    assert require_terminal_job(text, "24220", 16)["State"] == "COMPLETED"
