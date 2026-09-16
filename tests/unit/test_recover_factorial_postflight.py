import pytest

from benchmark_tools.recover_factorial_postflight import require_known_failure


def evidence():
    accounting = "JobID|JobIDRaw|State|ExitCode|Elapsed\n21248_0|21249|FAILED|1:0|00:26:37\n"
    status = {"provenance": {"slurm_job_id": "21249", "slurm_array_task_id": "0"},
              "status": "finished_pending_native_validation", "failed_methods": [],
              "accuracy_evaluated": False, "native_outputs_validated": False}
    log = ('Traceback (most recent call last):\n  File "runner", line 117, in main\n'
           '    verify_environment(environment)\nValueError: Package inventory changed: orthofinder\n')
    return accounting, status, log


def test_known_failure_retains_scheduler_failure():
    accounting, status, log = evidence()
    assert require_known_failure(accounting, 0, status, log)["State"] == "FAILED"


@pytest.mark.parametrize("change", ["running", "duplicate", "other_exit", "job", "index", "method",
                                    "scored", "native", "preflight", "other_error", "extra_traceback"])
def test_rejects_other_failure_or_unprovenance(change):
    accounting, status, log = evidence()
    if change == "running":
        accounting = accounting.replace("FAILED", "RUNNING")
    elif change == "duplicate":
        accounting += accounting.splitlines()[1] + "\n"
    elif change == "other_exit":
        accounting = accounting.replace("1:0", "2:0")
    elif change == "job":
        status["provenance"]["slurm_job_id"] = "other"
    elif change == "index":
        status["provenance"]["slurm_array_task_id"] = "1"
    elif change == "method":
        status["failed_methods"] = ["p0_c0_r1"]
    elif change == "scored":
        status["accuracy_evaluated"] = True
    elif change == "native":
        status["native_outputs_validated"] = True
    elif change == "preflight":
        log = log.replace("line 117", "line 92")
    elif change == "other_error":
        log = log.replace("orthofinder", "orthohmm")
    else:
        log = "Traceback (most recent call last):\n" + log
    with pytest.raises(ValueError):
        require_known_failure(accounting, 0, status, log)
