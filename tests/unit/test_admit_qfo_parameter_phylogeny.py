from copy import deepcopy

import pytest

from benchmark_tools import admit_qfo_parameter_phylogeny as module


def accounting(state="COMPLETED", exit_code="0:0", node="bizon", cpus="32"):
    return ("JobID|JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n"
            f"22034_0|23000|{state}|{exit_code}|01:00:00|{node}|{cpus}\n")


def test_completed_array_task_uses_raw_job_identity():
    row = module.completed_task(accounting(), 0)
    assert row["JobIDRaw"] == "23000"


def test_original_or_unrelated_array_cannot_admit_replacement():
    assert module.JOB == "22034"
    for job in ("21932", "22035"):
        with pytest.raises(ValueError, match="completed"):
            module.completed_task(accounting().replace("22034_0", job + "_0"), 0)


def test_native_admission_queries_replacement_before_reading_outputs(tmp_path, monkeypatch):
    def query(command, **kwargs):
        assert command[:3] == ["sacct", "-j", "22034"]
        return accounting(state="RUNNING")

    monkeypatch.setattr(module.subprocess, "check_output", query)
    monkeypatch.setattr(module, "verify_sources", lambda *a: pytest.fail("Read unfinished output"))
    with pytest.raises(ValueError, match="completed"):
        module.admit(tmp_path, 0)


@pytest.mark.parametrize("change", [{"state": "RUNNING"}, {"state": "FAILED"},
    {"exit_code": "1:0"}, {"node": "spark-7ff0"}, {"cpus": "2"}])
def test_wrong_accounting_rejected(change):
    with pytest.raises(ValueError, match="completed"):
        module.completed_task(accounting(**change), 0)


@pytest.mark.parametrize("index", [-1, 4, True, "0", 0.0])
def test_invalid_index(index):
    with pytest.raises(ValueError):
        module.completed_task(accounting(), index)


def test_duplicate_or_missing_task_rejected():
    with pytest.raises(ValueError):
        module.completed_task(accounting() + accounting().splitlines()[1] + "\n", 0)
    with pytest.raises(ValueError):
        module.completed_task(accounting(), 1)


def test_live_job_rejected_before_artifact_access(tmp_path, monkeypatch):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: accounting(state="RUNNING"))
    monkeypatch.setattr(module, "verify_sources", lambda *a: pytest.fail("Read unfinished output"))
    with pytest.raises(ValueError, match="completed"):
        module.admit(tmp_path, 0)


def execution_fixture():
    cell = {"label": "candidate_norm_low", "argv": ["python", "frozen.py"]}
    manifest = {"input_fastas": [{"path": "input.fa", "sha256": "abc", "bytes": 3}]}
    expected = {"job_id": "23000", "array_task_id": "0", "cell": cell}
    status = {"provenance": deepcopy(expected), "verified_inputs": {"status": "ready", "inputs": [
        {**manifest["input_fastas"][0], "absolute_path": "input.fa"}]},
        "dataset": cell["label"], "methods": {cell["label"]: {}},
        "status": "finished_pending_native_validation", "failed_methods": [],
        "accuracy_evaluated": False, "native_outputs_validated": False}
    post = {"status": "complete_pending_native_validation", "cell": cell,
            "accuracy_evaluated": False, "native_outputs_validated": False}
    return status, deepcopy(expected), post, expected, manifest, cell


def test_exact_unscored_execution_accepted():
    module.check_execution(*execution_fixture())


@pytest.mark.parametrize("problem", ["preflight", "provenance", "inputs", "dataset", "methods",
    "status", "failed", "accuracy", "native", "postflight", "post_cell", "post_scored"])
def test_inconsistent_execution_rejected(problem):
    status, pre, post, expected, manifest, cell = execution_fixture()
    if problem == "preflight":
        pre["job_id"] = "other"
    elif problem == "provenance":
        status["provenance"]["array_task_id"] = "1"
    elif problem == "inputs":
        status["verified_inputs"]["inputs"] = []
    elif problem == "dataset":
        status["dataset"] = "other"
    elif problem == "methods":
        status["methods"]["extra"] = {}
    elif problem == "status":
        status["status"] = "running"
    elif problem == "failed":
        status["failed_methods"] = [cell["label"]]
    elif problem == "accuracy":
        status["accuracy_evaluated"] = True
    elif problem == "native":
        status["native_outputs_validated"] = True
    elif problem == "postflight":
        post["status"] = "failed"
    elif problem == "post_cell":
        post["cell"] = {"label": "other"}
    else:
        post["accuracy_evaluated"] = True
    with pytest.raises(ValueError):
        module.check_execution(status, pre, post, expected, manifest, cell)


def test_native_adaptation_preserves_seed_and_other_arms():
    manifest = {"input_fastas": [{"path": "fasta"}], "candidate_arms": {
        "p1_c1": {"seed_partition": {"path": "seed"}, "candidate_partition": {"path": "baseline"},
                  "membership_constraints": {"path": "old_constraints"}},
        "p0_c0": {"candidate_partition": {"path": "other"}}}}
    before = deepcopy(manifest)
    arm = {"partition": {"path": "new_candidate"}, "constraints": {"path": "new_constraints"}}
    adapted = module.adapted_manifest(manifest, arm, {"path": "launcher"})
    assert adapted["candidate_arms"]["p1_c1"] == {
        "seed_partition": {"path": "seed"}, "candidate_partition": arm["partition"],
        "membership_constraints": arm["constraints"]}
    assert adapted["candidate_arms"]["p0_c0"] == before["candidate_arms"]["p0_c0"]
    assert adapted["fasta_inputs"] == manifest["input_fastas"]
    assert manifest == before
