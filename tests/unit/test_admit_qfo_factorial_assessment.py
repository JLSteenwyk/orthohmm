from copy import deepcopy

import pytest

from benchmark_tools.admit_qfo_factorial_assessment import check_completion


def fixture():
    preflight = {"status": "running", "cell": "p0_c1_r0", "index": 2,
                 "job_id": "21681", "accuracy_admitted": False, "stage": {"filtered_pairs": "frozen"}}
    report = {**deepcopy(preflight), "status": "process_succeeded_pending_independent_admission", "exit_code": 0}
    scheduler = {"JobIDRaw": "21681", "State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon", "AllocCPUS": "8"}
    return report, preflight, scheduler


def test_exact_fresh_completion():
    check_completion(*fixture(), 2)


@pytest.mark.parametrize("key,value", [("State", "RUNNING"), ("ExitCode", "1:0"), ("JobIDRaw", "other"),
                                      ("NodeList", "spark-7ff0"), ("AllocCPUS", "32")])
def test_wrong_scheduler_rejected(key, value):
    report, preflight, scheduler = fixture()
    scheduler[key] = value
    with pytest.raises(ValueError):
        check_completion(report, preflight, scheduler, 2)


@pytest.mark.parametrize("key,value", [("status", "running"), ("exit_code", 1), ("accuracy_admitted", True),
                                      ("cell", "p1_c1_r0"), ("index", 6), ("stage", {"filtered_pairs": "other"})])
def test_wrong_report_rejected(key, value):
    report, preflight, scheduler = fixture()
    report[key] = value
    with pytest.raises(ValueError):
        check_completion(report, preflight, scheduler, 2)


def test_changed_preflight_rejected():
    report, preflight, scheduler = fixture()
    preflight["status"] = "complete"
    with pytest.raises(ValueError):
        check_completion(report, preflight, scheduler, 2)


@pytest.mark.parametrize("index", [0, 4])
def test_reused_assessment_cannot_pass_as_fresh(index):
    with pytest.raises(ValueError, match="Reused"):
        check_completion(*fixture(), index)
