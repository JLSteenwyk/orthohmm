from copy import deepcopy
import pytest

from benchmark_tools.recover_native_factorial_receipt_22427 import receipt_failure


def evidence():
    run = dict(output_root="/retained/run_00", index=0, cell="p0_c0_r0", native_order=["s.fa"])
    plan = dict(path="/plan", sha256="digest")
    metrics = dict(status="complete", counts=dict(genes=2), stages=dict(search={}))
    initial = dict(status="native_factorial_running", plan=plan, index=0, cell=run["cell"],
                   native_order=run["native_order"], automatic_retry=False, factors={})
    pending = dict(initial, status="native_factorial_completed_pending_output_review", counts=metrics["counts"], stages=["search"])
    log = "save(report, receipt)\nFileExistsError: [Errno 17] File exists: '/retained/run_00/native_execution.json.pending' -> '/retained/run_00/native_execution.json'\n"
    return initial, pending, metrics, run, plan, log


def test_exact_post_completion_failure_identified_without_rewriting():
    data = evidence()
    prior = deepcopy(data)
    receipt_failure(*data)
    assert data == prior


@pytest.mark.parametrize("change", ["running", "metrics", "counts", "stage", "identity", "order", "retry", "other_failure", "other_path"])
def test_incomplete_or_unrelated_failure_cannot_be_recovered(change):
    initial, pending, metrics, run, plan, log = evidence()
    if change == "running": pending["status"] = "native_factorial_running"
    elif change == "metrics": metrics["status"] = "failed"
    elif change == "counts": pending["counts"] = dict(genes=999)
    elif change == "stage": pending["stages"] = ["search", "phylogeny"]
    elif change == "identity": pending["index"] = 1
    elif change == "order": pending["native_order"] = ["other.fa"]
    elif change == "retry": initial["automatic_retry"] = pending["automatic_retry"] = True
    elif change == "other_failure": log = log.replace("FileExistsError", "RuntimeError")
    else: log = log.replace("/retained/run_00", "/another/run")
    with pytest.raises(ValueError):
        receipt_failure(initial, pending, metrics, run, plan, log)
