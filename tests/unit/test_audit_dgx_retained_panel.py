import pytest

from benchmark_tools.audit_dgx_retained_panel import scheduler_rows


def fixture():
    snapshot = {"runs": [{"JobID": f"21656_{i}", "JobIDRaw": str(21657 + i) if i < 26 else "21656",
        "State": "COMPLETED" if i < 26 else "RUNNING", "ExitCode": "0:0", "AllocCPUS": "20",
        "ReqMem": "96G", "NodeList": "spark-7ff0"} for i in range(27)]}
    final = {"preceding_successful_controller_observation": {"exit_code": 0, "ArrayJobId": "21656",
        "ArrayTaskId": "26", "JobState": "COMPLETED", "ExitCode": "0:0", "JobId": "21656",
        "NodeList": "spark-7ff0", "NumCPUs": 20, "AllocTRES": "cpu=20,mem=96G,node=1,billing=20"}}
    return snapshot, final


def test_retained_provenance_keeps_final_controller_distinct():
    snapshot, final = fixture()
    rows = scheduler_rows(snapshot, final)
    assert len(rows) == 27
    assert rows["21656_0"]["origin"] == "retained_sacct_snapshot"
    assert rows["21656_26"]["origin"] == "transcribed_terminal_controller_fields_not_fresh_sacct"
    assert rows["21656_26"]["contract_row"] == ["21656_26", "21656", "COMPLETED", "0:0", "20", "96G", "spark-7ff0"]
    assert snapshot["runs"][-1]["State"] == "RUNNING"


@pytest.mark.parametrize("problem", ["missing", "duplicate", "failed", "final_failed", "task", "rawid", "cpu", "memory"])
def test_wrong_scheduler_evidence_rejected(problem):
    snapshot, final = fixture()
    observed = final["preceding_successful_controller_observation"]
    if problem == "missing":
        snapshot["runs"].pop()
    elif problem == "duplicate":
        snapshot["runs"][-1] = snapshot["runs"][0]
    elif problem == "failed":
        snapshot["runs"][0]["State"] = "FAILED"
    elif problem == "final_failed":
        observed["ExitCode"] = "1:0"
    elif problem == "task":
        observed["ArrayTaskId"] = "25"
    elif problem == "rawid":
        observed["JobId"] = "999"
    elif problem == "cpu":
        observed["NumCPUs"] = 32
    elif problem == "memory":
        snapshot["runs"][-1]["ReqMem"] = "64G"
    with pytest.raises(ValueError):
        scheduler_rows(snapshot, final)
