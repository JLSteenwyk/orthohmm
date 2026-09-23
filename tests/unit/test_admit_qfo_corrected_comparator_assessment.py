import pytest

from benchmark_tools.admit_qfo_corrected_comparator_assessment import validate_completion


def fixture():
    preflight = {"status": "running", "job_id": "1", "accuracy_admitted": False, "command": ["frozen"]}
    report = {**preflight, "status": "process_succeeded_pending_independent_admission", "exit_code": 0}
    scheduler = {"State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon", "AllocCPUS": "8", "JobIDRaw": "1"}
    return report, preflight, scheduler


def test_valid_terminal_assessment():
    validate_completion(*fixture())


@pytest.mark.parametrize("key,value", [("status", "running"), ("status", "failed"),
    ("exit_code", 1), ("accuracy_admitted", True), ("job_id", "other"), ("command", ["changed"])])
def test_reject_report_drift(key, value):
    args = fixture()
    args[0][key] = value
    with pytest.raises(ValueError):
        validate_completion(*args)


@pytest.mark.parametrize("key,value", [("State", "RUNNING"), ("ExitCode", "1:0"),
    ("NodeList", "other"), ("AllocCPUS", "2"), ("JobIDRaw", "2")])
def test_reject_scheduler(key, value):
    args = fixture()
    args[2][key] = value
    with pytest.raises(ValueError, match="scheduler"):
        validate_completion(*args)


def test_reject_preflight_status():
    args = fixture()
    args[1]["status"] = "prepared_unrun"
    with pytest.raises(ValueError, match="preflight"):
        validate_completion(*args)
def test_fastoma_replacement_executor_is_explicitly_pinned(tmp_path):
    from benchmark_tools.admit_qfo_corrected_comparator_assessment import executor_identity
    path, commit = executor_identity(tmp_path, "fastoma")
    assert path == tmp_path / "benchmarks/work/publication_qfo_corrected_fastoma_assessment_v2"
    assert commit == "0cc0a96c44e377f4e87a1e579ce12012d3478e4e"
