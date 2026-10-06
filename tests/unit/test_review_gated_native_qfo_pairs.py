"""Deferred conversion gates, with expensive native conversion explicitly stubbed."""

import json
from pathlib import Path
import shutil

import pytest

from benchmark_tools import run_review_gated_native_qfo_pairs as gate
from benchmark_tools.prepare_ob_candidate_neighborhood import record


@pytest.fixture
def fixture(tmp_path, monkeypatch):
    monkeypatch.setattr(gate, "ROOT", tmp_path)
    work = tmp_path / "benchmarks/work"
    work.mkdir(parents=True)
    tools = tmp_path / "benchmark_tools"
    tools.mkdir()
    reviewer = tools / "review_native_factorial_attempt.py"
    shutil.copyfile(Path(gate.conversion.__file__).with_name(reviewer.name), reviewer)
    batch = tmp_path / "batch.sh"
    batch.write_text("bound review batch\n")
    plan = tmp_path / "plan.json"
    plan.write_text("{}\n")
    request_path = tmp_path / "request.json"
    request_path.write_text(json.dumps(dict(job_id=101, index=8, plan=record(plan))))
    request_ref = record(request_path)
    run = dict(index=8, dataset="qfo_corrected", cell="p0_c1_r0", repeat=0)
    monkeypatch.setattr(gate.conversion, "validate_request", lambda req, ref, job: None)
    monkeypatch.setattr(gate.conversion, "validate_plan", lambda p: [None] * 8 + [run])
    review_dir = work / "native_factorial_terminal_review_101"
    review_dir.mkdir()
    review = dict(schema="native_factorial_terminal_review_v1", request=request_ref,
        plan=record(plan), job_id=101, index=8, dataset="qfo_corrected", cell="p0_c1_r0",
        repeat=0, status="native_success", scheduler_state="COMPLETED", scheduler_exit_code="0:0",
        terminal_reviewed=True, native_outputs_validated=True, primary_resources_replayed=True,
        shared_host_resources_reviewed=True, execution_scope=gate.conversion.SCOPE,
        resource_scopes=gate.conversion.SCOPES, uncontended_timing=False, automatic_retry=False,
        reviews=dict.fromkeys(("runtime", "resources", "environment", "outputs_or_failure")))
    review_path = review_dir / "review.json"
    review_path.write_text(json.dumps(review))
    submission_path = tmp_path / "submission.json"
    submission = dict(schema="native_factorial_review_held_submission_v1", job_id=102,
        native_job_id=101, index=8, automatic_retry=False, accuracy_evaluated=False,
        request=request_ref, destination=str(review_dir), review_source=record(reviewer), batch=record(batch))
    submission_path.write_text(json.dumps(submission))
    fields = dict(JobIDRaw="102", State="COMPLETED", ExitCode="0:0", NodeList="bizon", AllocCPUS="2", ReqMem="32Gn")
    monkeypatch.setattr(gate, "accounting", lambda job, include_memory: ("accounting", fields))
    monkeypatch.setattr(gate, "runtime_environment", lambda sub: {"fixture": True})
    context = dict(available_memory_bytes=64 * 2**30, available_disk_bytes=256 * 2**30)
    monkeypatch.setattr(gate, "capacity", lambda: context)
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
    monkeypatch.setenv("SLURM_JOB_ID", "103")
    calls = []

    def convert(request, review, destination):
        calls.append((request, review, destination))
        destination.mkdir()
        result = dict(status="full_native_factorial_qfo_pairs_prepared_unscored", native_job_id=101,
            native_index=8, cell="p0_c1_r0", job_id="103", source=record(gate.conversion.__file__),
            terminal_review=review, accuracy_evaluated=False, next_identity_authorized=False)
        output = destination / "results.json"
        output.write_text(json.dumps(result))
        return record(output)

    monkeypatch.setattr(gate.conversion, "prepare", convert)
    return dict(submission=submission, submission_path=submission_path, fields=fields, review=review,
                review_path=review_path, context=context, calls=calls, convert=convert,
                directory=work / "native_factorial_qfo_conversion_gate_101",
                destination=work / "native_factorial_qfo_pairs_101")


def execute(fixture):
    return gate.execute(record(fixture["submission_path"]), record(gate.__file__)["sha256"])


def result(fixture):
    return json.loads((fixture["directory"] / "results.json").read_text())


def test_success_invokes_only_the_original_converter_and_remains_unscored(fixture):
    ref = execute(fixture)
    report = result(fixture)
    assert ref == record(fixture["directory"] / "results.json")
    assert len(fixture["calls"]) == 1
    assert report["status"] == "review_gated_native_qfo_pairs_prepared_unscored"
    assert report["conversion_started"] is True and report["conversion_kind"] == "group"
    for flag in ("accuracy_evaluated", "native_inference_reexecuted", "next_identity_authorized",
                 "automatic_retry", "publication_ready"):
        assert report[flag] is False
    assert report["terminal_review"] == record(fixture["review_path"])
    assert report["conversion"] == record(fixture["destination"] / "results.json")


@pytest.mark.parametrize("key,value", (("JobIDRaw", "999"), ("State", "RUNNING"),
    ("State", "FAILED"), ("ExitCode", "1:0"), ("NodeList", "other"),
    ("AllocCPUS", "8"), ("ReqMem", "128Gn")))
def test_incomplete_or_wrong_reviewer_refuses_before_conversion(fixture, key, value):
    fixture["fields"][key] = value
    with pytest.raises(ValueError, match="Reviewer is not successfully completed"):
        execute(fixture)
    assert fixture["calls"] == [] and not fixture["destination"].exists()
    assert result(fixture)["conversion_started"] is False
    assert result(fixture)["status"] == "review_gated_native_qfo_conversion_failed_retained"


@pytest.mark.parametrize("key,value", (("status", "native_failure_retained"),
    ("native_outputs_validated", False), ("primary_resources_replayed", False),
    ("job_id", 999), ("cell", "p0_c0_r0")))
def test_completed_reviewer_does_not_bypass_native_admission(fixture, key, value):
    fixture["review"][key] = value
    fixture["review_path"].write_text(json.dumps(fixture["review"]))
    with pytest.raises(ValueError, match="bound successful"):
        execute(fixture)
    assert fixture["calls"] == [] and result(fixture)["conversion_started"] is False


@pytest.mark.parametrize("key", ("available_memory_bytes", "available_disk_bytes"))
def test_unsafe_future_capacity_does_not_launch_conversion(fixture, key):
    fixture["context"][key] = 1
    with pytest.raises(ValueError, match="Unsafe conversion"):
        execute(fixture)
    assert fixture["calls"] == [] and not fixture["destination"].exists()
    assert result(fixture)["conversion_started"] is False


@pytest.mark.parametrize("role", ("directory", "destination"))
def test_occupied_output_namespaces_are_not_overwritten(fixture, role):
    fixture[role].mkdir()
    marker = fixture[role] / "keep"
    marker.write_text("untouched")
    with pytest.raises(FileExistsError):
        execute(fixture)
    assert marker.read_text() == "untouched" and fixture["calls"] == []


def test_wrong_worker_hash_is_rejected_before_any_output(fixture):
    with pytest.raises(ValueError, match="worker source changed"):
        gate.execute(record(fixture["submission_path"]), "d" * 64)
    assert not fixture["directory"].exists() and fixture["calls"] == []


@pytest.mark.parametrize("key,value", (("index", True), ("index", 5), ("index", 13),
    ("automatic_retry", True), ("native_job_id", 102), ("destination", "/unbound")))
def test_wrong_submission_identity_never_creates_conversion_state(fixture, key, value):
    fixture["submission"][key] = value
    fixture["submission_path"].write_text(json.dumps(fixture["submission"]))
    with pytest.raises(ValueError):
        execute(fixture)
    assert not fixture["directory"].exists() and fixture["calls"] == []


def test_actual_runtime_checker_refuses_an_unrelated_python_invocation():
    with pytest.raises(ValueError, match="original Python3.10 venv"):
        gate.runtime_environment(dict(python_invocation_path="/not/the/retained/python"))


def test_same_native_or_review_job_cannot_be_used_for_conversion(fixture, monkeypatch):
    monkeypatch.setenv("SLURM_JOB_ID", "102")
    with pytest.raises(ValueError, match="separate scheduled"):
        execute(fixture)
    assert not fixture["directory"].exists()


def test_preparer_failure_is_retained_without_retry(fixture, monkeypatch):
    def fail(*args):
        fixture["calls"].append(args)
        raise RuntimeError("fixture converter failure")
    monkeypatch.setattr(gate.conversion, "prepare", fail)
    with pytest.raises(RuntimeError):
        execute(fixture)
    assert len(fixture["calls"]) == 1
    report = result(fixture)
    assert report["conversion_started"] is True and report["error_type"] == "RuntimeError"
    assert report["automatic_retry"] is False and report["accuracy_evaluated"] is False


def test_postflight_input_drift_retains_converted_output_but_refuses_gate(fixture, monkeypatch):
    def changed(*args):
        ref = fixture["convert"](*args)
        fixture["submission_path"].write_text("{}\n")
        return ref
    monkeypatch.setattr(gate.conversion, "prepare", changed)
    with pytest.raises(ValueError):
        execute(fixture)
    assert (fixture["destination"] / "results.json").exists()
    assert result(fixture)["status"] == "review_gated_native_qfo_conversion_failed_retained"
    assert result(fixture)["conversion_started"] is True


@pytest.mark.parametrize("key,value", (("native_job_id", 999), ("native_index", 9),
    ("cell", "p0_c1_r1"), ("job_id", "102"), ("source", {}),
    ("terminal_review", {}), ("accuracy_evaluated", True),
    ("next_identity_authorized", True), ("status", "failure")))
def test_unexpected_preparer_stage_is_retained_but_not_admitted(fixture, monkeypatch, key, value):
    def changed(*args):
        ref = fixture["convert"](*args)
        output = Path(ref["path"])
        stage = json.loads(output.read_text())
        stage[key] = value
        output.write_text(json.dumps(stage))
        return record(output)
    monkeypatch.setattr(gate.conversion, "prepare", changed)
    with pytest.raises(ValueError, match="Unexpected converted stage"):
        execute(fixture)
    assert len(fixture["calls"]) == 1
    assert result(fixture)["status"] == "review_gated_native_qfo_conversion_failed_retained"
    assert result(fixture)["accuracy_evaluated"] is False


@pytest.mark.parametrize("key,value", (("SLURM_CPUS_PER_TASK", "8"),
    ("SLURM_JOB_ID", ""), ("SLURM_JOB_ID", "101")))
def test_invalid_own_allocation_never_creates_output(fixture, monkeypatch, key, value):
    monkeypatch.setenv(key, value)
    with pytest.raises(ValueError, match="separate scheduled"):
        execute(fixture)
    assert not fixture["directory"].exists() and fixture["calls"] == []


@pytest.mark.parametrize("role", ("directory", "destination"))
def test_dangling_output_symlink_is_not_followed_or_removed(fixture, role):
    target = fixture[role].parent / "absent-target"
    fixture[role].symlink_to(target)
    with pytest.raises(ValueError, match="direct output namespaces"):
        execute(fixture)
    assert fixture[role].is_symlink() and not target.exists() and fixture["calls"] == []
