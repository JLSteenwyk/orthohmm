from copy import deepcopy
from pathlib import Path

import pytest

from benchmark_tools.validate_ob_candidate_phylogeny import adapted_prepared, check_execution, completed_task, check_baseline


@pytest.mark.parametrize("state,code", [("COMPLETED", "0:0"), ("RUNNING", "0:0"), ("FAILED", "1:0"), ("COMPLETED", "1:0")])
def test_scheduler_requires_terminal_success(state, code):
    text = f"JobID|JobIDRaw|State|ExitCode|Elapsed\n21316_0|21317|{state}|{code}|00:07:35\n"
    if state == "COMPLETED" and code == "0:0":
        assert completed_task(text, 0)["JobIDRaw"] == "21317"
    else:
        with pytest.raises(ValueError):
            completed_task(text, 0)


def test_missing_and_duplicate_tasks_rejected():
    header = "JobID|JobIDRaw|State|ExitCode|Elapsed\n"
    row = "21316_0|21317|COMPLETED|0:0|00:07:35\n"
    for text in (header, header + row + row):
        with pytest.raises(ValueError):
            completed_task(text, 0)


def test_cpm_scheduler_cannot_admit_threshold_task():
    header = "JobID|JobIDRaw|State|ExitCode|Elapsed\n"
    threshold = "21316_0|21317|COMPLETED|0:0|00:07:35\n"
    cpm = "21324_0|21325|COMPLETED|0:0|00:07:35\n"
    assert completed_task(header + threshold + cpm, 0, cpm=True)["JobIDRaw"] == "21325"
    with pytest.raises(ValueError):
        completed_task(header + threshold, 0, cpm=True)


@pytest.mark.parametrize("label", ["cpm_low", "cpm_high"])
def test_cpm_adapter_retains_own_seed(label):
    prepared = {"candidate_arms": {"p1_c1": {"seed_partition": "baseline"}}}
    arm = {"label": label, "seed_partition": {"path": "own_seed"},
           "partition": {"path": "own_candidates"}, "constraints": {"path": "own_trace"}}
    result = adapted_prepared(prepared, arm)
    assert result["candidate_arms"]["p1_c1"] == {
        "seed_partition": arm["seed_partition"], "candidate_partition": arm["partition"],
        "membership_constraints": arm["constraints"]}
    assert prepared["candidate_arms"]["p1_c1"]["seed_partition"] == "baseline"


def test_prepared_adapter_changes_only_candidate_inputs():
    prepared = {"candidate_arms": {"p1_c1": {"candidate_partition": "old", "membership_constraints": "old",
                    "seed_partition": "seed"}}, "fasta_inputs": ["fasta"]}
    old = deepcopy(prepared)
    arm = {"label": "norm_low", "partition": {"path": "new"}, "constraints": {"path": "newtrace"}}
    adapted = adapted_prepared(prepared, arm)
    assert prepared == old
    assert adapted["candidate_arms"]["p1_c1"] == {"candidate_partition": arm["partition"],
        "membership_constraints": arm["constraints"], "seed_partition": "seed"}
    assert adapted["fasta_inputs"] == ["fasta"]


def test_baseline_checks_native_evidence_not_validator_location():
    current = {key: key for key in ("cell", "status", "native_manifest", "native_metrics", "species_tree",
        "partition", "membership", "native_outputs_validated", "accuracy_evaluated")}
    prior = {**current, "verifier": "frozen/validator.py"}
    check_baseline({**current, "verifier": "main/validator.py"}, prior)
    for key in current:
        with pytest.raises(ValueError, match="source changed"):
            check_baseline({**current, key: "changed"}, prior)


@pytest.mark.parametrize("cpm", [False, True])
@pytest.mark.parametrize("problem", [None, "status", "command", "postflight", "raw_job", "scored", "provenance", "array_job"])
def test_execution_identity(problem, cpm):
    cell, arm = {"label": "candidate_norm_low", "argv": ["infer"]}, {"label": "norm_low"}
    pre = {"job_id": "21317", "array_job_id": "21316", "array_task_id": "0", "cell": cell,
           "arm": arm, "cwd": "/launcher"}
    if cpm:
        pre["array_job_id"] = "21324"
        cell["label"] = "candidate_cpm_low"
        arm["label"] = "cpm_low"
    post = {"status": "complete_pending_native_validation", "accuracy_evaluated": False,
            "baseline_source_unchanged": True, "cell": cell}
    status = {"status": "finished_pending_native_validation", "failed_methods": [], "accuracy_evaluated": False,
        "native_outputs_validated": False, "dataset": cell["label"], "methods": {cell["label"]: {}},
        "provenance": deepcopy(pre)}
    if problem == "status":
        status["failed_methods"] = [cell["label"]]
    elif problem == "command":
        pre["cell"] = {"label": "other"}
    elif problem == "postflight":
        post["baseline_source_unchanged"] = False
    elif problem == "raw_job":
        pre["job_id"] = "21316"
        status["provenance"] = deepcopy(pre)
    elif problem == "scored":
        status["accuracy_evaluated"] = True
    elif problem == "provenance":
        status["provenance"] = {}
    elif problem == "array_job":
        pre["array_job_id"] = "21316" if cpm else "21324"
        status["provenance"] = deepcopy(pre)
    args = (status, pre, post, cell, arm, {"JobIDRaw": "21317"}, 0, Path("/launcher"), cpm)
    if problem:
        with pytest.raises(ValueError):
            check_execution(*args)
    else:
        check_execution(*args)
