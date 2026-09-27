import copy

import pytest

from benchmark_tools.admit_qfo_canonical_assessment import verify_completion


def fixture():
    plan, command, source = {"sha256": "plan"}, ["frozen", "command"], {"sha256": "source"}
    preflight = dict(plan=plan, status="running", job_id="22336", source=source,
                     accuracy_admitted=False, command=command)
    report = dict(preflight, status="process_succeeded_pending_independent_admission", exit_code=0)
    return report, copy.deepcopy(preflight), plan, command, source


def test_completed_binding():
    verify_completion(*fixture())


@pytest.mark.parametrize("index,key,value", [(0,"job_id","other"), (0,"exit_code",1),
    (0,"status","running"), (0,"accuracy_admitted",True), (0,"plan",{}),
    (0,"source",{}), (0,"command",[]), (1,"status","complete"), (1,"extra",True)])
def test_wrong_completion_rejected(index, key, value):
    args = fixture()
    args[index][key] = value
    with pytest.raises(ValueError):
        verify_completion(*args)
