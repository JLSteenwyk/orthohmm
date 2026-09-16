from copy import deepcopy

import pytest

from benchmark_tools.validate_unconstrained_control import check_status


@pytest.mark.parametrize("change", [None, "job", "failed", "scored", "parent", "constraints", "method"])
def test_diagnostic_status_admission(change):
    cell = {"label": "p1_c1_r1_unconstrained_v2", "omitted_membership_constraints": "/frozen/constraints"}
    status = {"status": "finished_pending_native_validation", "failed_methods": [],
              "accuracy_evaluated": False, "native_outputs_validated": False,
              "dataset": cell["label"], "methods": {cell["label"]: {}},
              "provenance": {"slurm_job_id": "21290", "slurm_array_task_id": None, "parent_cell": "p1_c1_r1",
                             "omitted_membership_constraints": cell["omitted_membership_constraints"]}}
    status = deepcopy(status)
    if change == "job":
        status["provenance"]["slurm_job_id"] = "21289"
    elif change == "failed":
        status["failed_methods"] = [cell["label"]]
    elif change == "scored":
        status["accuracy_evaluated"] = True
    elif change == "parent":
        status["provenance"]["parent_cell"] = "p0_c0_r1"
    elif change == "constraints":
        status["provenance"]["omitted_membership_constraints"] = "/other"
    elif change == "method":
        status["methods"]["other"] = {}
    if change is None:
        check_status(status, cell)
    else:
        with pytest.raises(ValueError):
            check_status(status, cell)
