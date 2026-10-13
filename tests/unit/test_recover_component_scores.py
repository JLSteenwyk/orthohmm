import pytest

from benchmark_tools.recover_component_scores import require_completed_probe


def accounting(job="24219", state="COMPLETED", exit_code="0:0", cpus="1"):
    return f"JobID|State|ExitCode|AllocCPUS\n{job}|{state}|{exit_code}|{cpus}\n"


def test_clean_probe_requires_successful_terminal_one_cpu_record():
    assert require_completed_probe(accounting(), "24219")["State"] == "COMPLETED"


@pytest.mark.parametrize("kwargs", [{"job": "wrong"}, {"state": "RUNNING"}, {"state": "FAILED"},
    {"exit_code": "0:11"}, {"cpus": "16"}])
def test_other_jobs_failed_or_live_jobs_and_wrong_allocations_cannot_authorize_recovery(kwargs):
    with pytest.raises(ValueError):
        require_completed_probe(accounting(**kwargs), "24219")


def test_missing_or_duplicate_scheduler_record_cannot_authorize_recovery():
    with pytest.raises(ValueError):
        require_completed_probe("JobID|State|ExitCode|AllocCPUS\n", "24219")
    with pytest.raises(ValueError):
        require_completed_probe(accounting() + "24219|COMPLETED|0:0|1\n", "24219")
