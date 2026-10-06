"""Deferred assessment gates; expensive endpoint execution is explicitly stubbed."""

import json
from pathlib import Path

import pytest

from benchmark_tools import run_conversion_gated_native_qfo_assessment as gate
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def put(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))
    return record(path)


@pytest.fixture
def fixture(tmp_path, monkeypatch):
    monkeypatch.setattr(gate, "ROOT", tmp_path)
    (tmp_path / "benchmarks/work").mkdir(parents=True)
    plan = put(tmp_path / "plan.json", {})
    request_ref = put(tmp_path / "request.json", dict(job_id=101, index=8, plan=plan))
    request = json.loads(Path(request_ref["path"]).read_text())
    run = dict(index=8, dataset="qfo_corrected", cell="p0_c1_r0", repeat=0)
    reviewer_ref = put(tmp_path / "reviewer_submission.json", {})
    reviewer = dict(job_id=102, native_job_id=101)
    monkeypatch.setattr(gate.conversion_gate, "submission_binding",
                        lambda ref: (reviewer, request_ref, request, run))
    batch = put(tmp_path / "conversion.sh", "fixture")
    conversion_paths = [tmp_path / "benchmarks/work" / f"{prefix}_101" for prefix in
        ("native_factorial_qfo_conversion_gate", "native_factorial_qfo_pairs")]
    submission = dict(schema="review_gated_native_qfo_conversion_held_submission_v1", job_id=103,
        held_inspection_passed=True, pair_conversion_completed=False, accuracy_evaluated=False,
        automatic_retry=False, next_identity_authorized=False, publication_ready=False,
        native_job_id=101, reviewer_job_id=102, index=8, cell="p0_c1_r0", request=request_ref,
        plan=plan, reviewer_submission=reviewer_ref, worker=record(gate.conversion_gate.__file__),
        converter=record(gate.conversion_gate.conversion.__file__), batch=batch,
        output_namespaces=[str(p) for p in conversion_paths])
    submission_path = tmp_path / "submission.json"
    put(submission_path, submission)
    review_ref = put(tmp_path / "review.json", {})
    pairs = put(conversion_paths[1] / "pairs.tsv", "A1 B1")
    stage = dict(schema="full_native_factorial_qfo_conversion_v1",
        status="full_native_factorial_qfo_pairs_prepared_unscored", native_index=8,
        cell="p0_c1_r0", native_job_id=101, job_id="103", request=request_ref, plan=plan,
        participant="ohmm_qfo_full_native_p0_c1_r0", conversion_kind="group",
        semantics="cross-species group-derived clique pairs", accuracy_evaluated=False,
        native_inference_reexecuted=False, automatic_retry=False, next_identity_authorized=False,
        publication_ready=False, total_pairs=2, retained_pairs=2, expected_pairs=2,
        removed_mapping_pairs=0, empty_predictions=False, pairs=pairs, filtered_pairs=pairs,
        pair_coverage=dict(pair_rows=2, input_accessions=3, accessions_in_any_pair=3,
                           fraction_inputs_in_any_pair=1.), source=submission["converter"],
        terminal_review=review_ref)
    pairs_ref = put(conversion_paths[1] / "results.json", stage)
    conversion = dict(schema="review_gated_native_qfo_conversion_v1",
        status="review_gated_native_qfo_pairs_prepared_unscored", source=submission["worker"],
        converter=submission["converter"], reviewer_submission=reviewer_ref, request=request_ref,
        job_id="103", native_job_id=101, reviewer_job_id=102, index=8, cell="p0_c1_r0",
        destination=submission["output_namespaces"][1], conversion=pairs_ref, conversion_started=True,
        native_inference_reexecuted=False, accuracy_evaluated=False, next_identity_authorized=False,
        automatic_retry=False, publication_ready=False, terminal_review=review_ref)
    conversion_path = conversion_paths[0] / "results.json"
    put(conversion_path, conversion)
    fields = dict(JobIDRaw="103", State="COMPLETED", ExitCode="0:0", NodeList="bizon",
                  AllocCPUS="2", ReqMem="32G")
    monkeypatch.setattr(gate.assessment, "accounting", lambda job, include_memory: ("fixture", fields))
    monkeypatch.setattr(gate.conversion_gate, "runtime_environment", lambda ref: {"fixture": True})
    context = dict(available_memory_bytes=128 * 2**30, available_disk_bytes=256 * 2**30)
    monkeypatch.setattr(gate.conversion_gate, "capacity", lambda: context)
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
    monkeypatch.setenv("SLURM_JOB_ID", "104")
    paths = gate.output_paths(run, 101)
    calls = []

    def assess(root, pairs_ref, conversion_job):
        calls.append((root, pairs_ref, conversion_job))
        execution = dict(status="process_succeeded_pending_independent_admission", job_id="104",
            source=record(gate.assessment.__file__), pairs_manifest=pairs_ref, exit_code=0,
            native_index=8, cell="p0_c1_r0", native_job_id=101, accuracy_admitted=False,
            **{k: str(paths[k]) for k in ("cwd", "work", "results")},
            publication_ready=False, automatic_retry=False, native_inference_reexecuted=False,
            next_identity_authorized=False)
        put(paths["cwd"] / "results.json", execution)
        return execution

    monkeypatch.setattr(gate.assessment, "run", assess)
    return dict(submission=submission, submission_path=submission_path, conversion=conversion,
        conversion_path=conversion_path, stage=stage, stage_path=Path(pairs_ref["path"]), fields=fields,
        context=context, paths=paths, calls=calls, assess=assess)


def execute(fixture):
    return gate.execute(record(fixture["submission_path"]), record(gate.__file__)["sha256"])


def report(fixture):
    return json.loads((fixture["paths"]["gate"] / "results.json").read_text())


def test_success_invokes_original_driver_once_and_still_requires_admission(fixture):
    ref = execute(fixture)
    assert ref == record(fixture["paths"]["gate"] / "results.json")
    result = report(fixture)
    assert result["status"] == "native_qfo_assessment_process_succeeded_pending_independent_admission"
    assert result["assessment_driver_invoked"] is result["endpoint_process_completed"] is True
    assert len(fixture["calls"]) == 1 and fixture["calls"][0][2] == "103"
    for key in ("accuracy_admitted", "publication_ready", "automatic_retry",
                "native_inference_reexecuted", "next_identity_authorized"):
        assert result[key] is False


@pytest.mark.parametrize("key,value", (("JobIDRaw", "999"), ("State", "RUNNING"),
    ("State", "FAILED"), ("ExitCode", "1:0"), ("NodeList", "other"),
    ("AllocCPUS", "8"), ("ReqMem", "128G")))
def test_unsuccessful_conversion_never_reads_unfinished_output(fixture, key, value):
    fixture["fields"][key] = value
    fixture["conversion_path"].unlink()
    fixture["stage_path"].unlink()
    with pytest.raises(ValueError, match="Conversion is not successfully completed"):
        execute(fixture)
    assert fixture["calls"] == [] and report(fixture)["assessment_driver_invoked"] is False


@pytest.mark.parametrize("key,value", (("status", "conversion_failed"), ("source", {}),
    ("converter", {}), ("reviewer_submission", {}), ("request", {}), ("job_id", "999"),
    ("native_job_id", 999), ("reviewer_job_id", 999), ("index", 9), ("cell", "p0_c1_r1"),
    ("destination", "/unbound"), ("conversion", {}), ("conversion_started", False),
    ("native_inference_reexecuted", True), ("accuracy_evaluated", True),
    ("next_identity_authorized", True), ("automatic_retry", True), ("publication_ready", True)))
def test_successful_scheduler_does_not_bypass_conversion_gate(fixture, key, value):
    fixture["conversion"][key] = value
    put(fixture["conversion_path"], fixture["conversion"])
    with pytest.raises(ValueError, match="successful bound conversion gate"):
        execute(fixture)
    assert fixture["calls"] == [] and report(fixture)["assessment_driver_invoked"] is False


@pytest.mark.parametrize("key,value", (("native_job_id", 999), ("request", {}),
    ("plan", {}), ("terminal_review", {}), ("native_index", 9), ("retained_pairs", 1)))
def test_changed_converted_stage_refuses_before_endpoint_execution(fixture, key, value):
    fixture["stage"][key] = value
    stage_ref = put(fixture["stage_path"], fixture["stage"])
    fixture["conversion"]["conversion"] = stage_ref
    put(fixture["conversion_path"], fixture["conversion"])
    with pytest.raises(ValueError):
        execute(fixture)
    assert fixture["calls"] == [] and report(fixture)["assessment_driver_invoked"] is False


@pytest.mark.parametrize("key", ("available_memory_bytes", "available_disk_bytes"))
def test_unsafe_future_capacity_retains_refusal(fixture, key):
    fixture["context"][key] = 1
    with pytest.raises(ValueError, match="Unsafe assessment"):
        execute(fixture)
    assert fixture["calls"] == [] and report(fixture)["assessment_driver_invoked"] is False


@pytest.mark.parametrize("role", ("gate", "cwd", "work", "results"))
def test_no_output_namespace_can_be_reused(fixture, role):
    fixture["paths"][role].mkdir(parents=True)
    marker = fixture["paths"][role] / "keep"
    marker.write_text("untouched")
    with pytest.raises(FileExistsError):
        execute(fixture)
    assert marker.read_text() == "untouched" and fixture["calls"] == []


@pytest.mark.parametrize("key,value", (("held_inspection_passed", False), ("job_id", True),
    ("job_id", 101), ("native_job_id", 999), ("index", 9), ("cell", "p0_c1_r1"),
    ("automatic_retry", True), ("pair_conversion_completed", True), ("worker", {}),
    ("converter", {}), ("output_namespaces", []), ("native_job_id", 101.0),
    ("reviewer_job_id", 102.0), ("index", 8.0), ("index", True)))
def test_wrong_conversion_submission_never_creates_assessment_state(fixture, key, value):
    fixture["submission"][key] = value
    put(fixture["submission_path"], fixture["submission"])
    with pytest.raises(ValueError):
        execute(fixture)
    assert not fixture["paths"]["gate"].exists() and fixture["calls"] == []


@pytest.mark.parametrize("key,value", (("SLURM_CPUS_PER_TASK", "2"), ("SLURM_JOB_ID", ""),
    ("SLURM_JOB_ID", "101"), ("SLURM_JOB_ID", "102"), ("SLURM_JOB_ID", "103")))
def test_wrong_own_allocation_refuses_before_any_output(fixture, monkeypatch, key, value):
    monkeypatch.setenv(key, value)
    with pytest.raises(ValueError, match="separate scheduled"):
        execute(fixture)
    assert not fixture["paths"]["gate"].exists()


def test_driver_exception_retained_without_retry(fixture, monkeypatch):
    def fail(*args):
        fixture["calls"].append(args)
        raise RuntimeError("fixture endpoint failure")
    monkeypatch.setattr(gate.assessment, "run", fail)
    with pytest.raises(RuntimeError):
        execute(fixture)
    assert len(fixture["calls"]) == 1
    assert report(fixture)["status"] == "conversion_gated_native_qfo_assessment_failed_retained"
    assert report(fixture)["endpoint_process_completed"] is False and report(fixture)["automatic_retry"] is False


@pytest.mark.parametrize("key,value", (("status", "failed"), ("job_id", "999"),
    ("source", {}), ("pairs_manifest", {}), ("exit_code", 1), ("native_index", 9),
    ("cell", "p0_c1_r1"), ("native_job_id", 999), ("accuracy_admitted", True),
    ("cwd", "/unbound"), ("work", "/unbound"), ("results", "/unbound")))
def test_wrong_driver_result_not_admitted(fixture, monkeypatch, key, value):
    def changed(*args):
        execution = fixture["assess"](*args)
        execution[key] = value
        put(fixture["paths"]["cwd"] / "results.json", execution)
        return execution
    monkeypatch.setattr(gate.assessment, "run", changed)
    with pytest.raises(ValueError, match="Unexpected original assessment result"):
        execute(fixture)
    assert (fixture["paths"]["cwd"] / "results.json").exists()
    assert report(fixture)["accuracy_admitted"] is report(fixture)["endpoint_process_completed"] is False


def test_postflight_input_drift_retains_scoring_but_refuses_gate(fixture, monkeypatch):
    def changed(*args):
        execution = fixture["assess"](*args)
        put(fixture["submission_path"], {})
        return execution
    monkeypatch.setattr(gate.assessment, "run", changed)
    with pytest.raises(ValueError):
        execute(fixture)
    assert (fixture["paths"]["cwd"] / "results.json").exists()
    assert report(fixture)["status"] == "conversion_gated_native_qfo_assessment_failed_retained"


def test_changed_worker_hash_never_creates_state(fixture):
    with pytest.raises(ValueError, match="Assessment gate worker changed"):
        gate.execute(record(fixture["submission_path"]), "f" * 64)
    assert not fixture["paths"]["gate"].exists()


@pytest.mark.parametrize("role", ("gate", "cwd", "work", "results"))
def test_dangling_assessment_namespace_symlink_remains_untouched(fixture, role):
    path = fixture["paths"][role]
    path.parent.mkdir(parents=True, exist_ok=True)
    target = path.parent / "absent-target"
    path.symlink_to(target)
    with pytest.raises(ValueError, match="direct assessment namespaces"):
        execute(fixture)
    assert path.is_symlink() and not target.exists() and fixture["calls"] == []
