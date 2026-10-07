"""Admission remains solely with the unchanged full reviewer."""

from copy import deepcopy
import inspect

import pytest

from benchmark_tools import finalize_native10_review as finalization


def values():
    flags = dict(accuracy_evaluated=False, terminal_reviewed=False,
                 next_identity_authorized=False, automatic_retry=False)
    partition = dict(flags, index=10, job_id=23902, status="diagnosis_completed",
                     frozen_root_coverage_gate=dict(status="passed"),
                     all_output_partitions_identical=True, source_payload_matches_native_input=True,
                     parsers_identical=True, stages=dict(root=dict(complete_unique_input_partition=True)))
    context = dict(flags, index=10, job_id=23902, status="semantic_probe_passed",
                   semantic_result=dict(native_outputs_validated=True))
    failure = dict(flags, index=10, job_id=23902, status="terminal_factorial_review_failed")
    request = dict(index=10, job_id=23902)
    return partition, context, failure, request


def test_one_distinct_review_preserves_failed_history_and_never_retries_inference():
    result = finalization.decision(*values(), 99999, "2")
    assert result["status"] == "one_shot_fresh_review_permitted"
    assert result["native_retry"] is False
    assert result["automatic_retry"] is False
    assert result["original_failure_retained"] is True
    assert result["historical_failure_cause_established"] is False


@pytest.mark.parametrize("job,cpus", [(23902, "2"), (23910, "2"), (0, "2"),
                                     (True, "2"), (99999, "1"), (99999, None)])
def test_wrong_allocation_identity_rejected(job, cpus):
    with pytest.raises(ValueError):
        finalization.decision(*values(), job, cpus)


@pytest.mark.parametrize("which,key,value", [
    (0, "index", 11), (0, "job_id", 1000), (0, "status", "failed"),
    (0, "frozen_root_coverage_gate", {"status": "failed"}),
    (0, "all_output_partitions_identical", False), (0, "source_payload_matches_native_input", False),
    (0, "parsers_identical", False), (0, "stages", {}),
    (0, "stages", {"root": {"complete_unique_input_partition": False}}),
    (1, "status", "semantic_probe_failed"), (1, "semantic_result", {"native_outputs_validated": False}),
    (2, "status", "native_success"), (2, "job_id", 1000), (3, "index", 11),
])
def test_failed_or_wrong_diagnosis_cannot_trigger_fresh_review(which, key, value):
    data = deepcopy(values())
    data[which][key] = value
    with pytest.raises(ValueError):
        finalization.decision(*data, 99999, "2")


@pytest.mark.parametrize("which", [0, 1, 2])
@pytest.mark.parametrize("key", ["accuracy_evaluated", "terminal_reviewed", "next_identity_authorized", "automatic_retry"])
def test_no_diagnostic_admission_promotion(which, key):
    data = deepcopy(values())
    data[which][key] = True
    with pytest.raises(ValueError):
        finalization.decision(*data, 99999, "2")


def test_prospective_wrapper_does_not_patch_or_replace_the_frozen_reviewer():
    source = inspect.getsource(finalization.finalize)
    assert "from benchmark_tools.review_allocated_native_factorial_attempt import review" in source
    assert "result = review(request_ref, CONTROL / \"review\")" in source
    assert "CONTROL.mkdir(exist_ok=False)" in source
    assert "Finalization identity already used; no retry or overwrite" in source
    assert "mock" not in source
    assert "setattr" not in source


def test_consumed_namespace_rejected_before_frozen_reviewer(tmp_path, monkeypatch):
    monkeypatch.setattr(finalization, "CONTROL", tmp_path)
    with pytest.raises(ValueError, match="identity already used"):
        finalization.finalize(99999, "2")
