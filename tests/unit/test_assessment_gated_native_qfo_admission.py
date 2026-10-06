"""Deferred admission gates; the expensive independent validator is stubbed."""

import json
from pathlib import Path

import pytest

from benchmark_tools import run_assessment_gated_native_qfo_admission as gate
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def put(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))
    return record(path)


@pytest.fixture
def fixture(tmp_path, monkeypatch):
    monkeypatch.setattr(gate, "ROOT", tmp_path)
    monkeypatch.setattr(gate.assessment_gate, "ROOT", tmp_path)
    (tmp_path / "benchmarks/work").mkdir(parents=True)
    plan = put(tmp_path / "plan.json", {})
    request = put(tmp_path / "request.json", dict(job_id=101, index=8, plan=plan))
    reviewer_ref = put(tmp_path / "reviewer_submission.json", {})
    reviewer = dict(job_id=102, native_job_id=101)
    run = dict(index=8, cell="p0_c1_r0", dataset="qfo_corrected", repeat=0)
    conversion_paths = [tmp_path / "benchmarks/work" / f"{prefix}_101" for prefix in
        ("native_factorial_qfo_conversion_gate", "native_factorial_qfo_pairs")]
    conversion = dict(job_id=103, native_job_id=101, reviewer_job_id=102, index=8,
        cell="p0_c1_r0", request=request, plan=plan, reviewer_submission=reviewer_ref,
        worker=record(gate.assessment_gate.conversion_gate.__file__),
        converter=record(gate.assessment_gate.conversion_gate.conversion.__file__),
        output_namespaces=[str(p) for p in conversion_paths])
    conversion_ref = put(tmp_path / "conversion_submission.json", conversion)
    monkeypatch.setattr(gate.assessment_gate, "submission_binding",
        lambda ref: (conversion, reviewer, {}, run, conversion_paths))
    paths = gate.assessment_gate.output_paths(run, 101)
    batch_ref = put(tmp_path / "assessment.sh", "fixture batch")
    sub = dict(schema="conversion_gated_native_qfo_assessment_held_submission_v1", job_id=104,
        native_job_id=101, reviewer_job_id=102, conversion_job_id=103, index=8, cell="p0_c1_r0",
        held_inspection_passed=True, assessment_completed=False, accuracy_admitted=False,
        automatic_retry=False, next_identity_authorized=False, publication_ready=False,
        conversion_submission=conversion_ref, reviewer_submission=reviewer_ref, request=request, plan=plan,
        worker=record(gate.assessment_gate.__file__), driver=record(gate.assessment_gate.assessment.__file__),
        batch=batch_ref, output_namespaces={k: str(p) for k, p in paths.items()})
    submission_path = tmp_path / "assessment_submission.json"
    put(submission_path, sub)
    pairs_ref = put(conversion_paths[1] / "results.json", dict(native_index=8, cell="p0_c1_r0"))
    conversion_gate = dict(schema="review_gated_native_qfo_conversion_v1",
        status="review_gated_native_qfo_pairs_prepared_unscored", source=conversion["worker"],
        converter=conversion["converter"], reviewer_submission=reviewer_ref, request=request,
        job_id="103", native_job_id=101, reviewer_job_id=102, index=8, cell="p0_c1_r0",
        destination=str(conversion_paths[1]), conversion=pairs_ref, conversion_started=True,
        native_inference_reexecuted=False, accuracy_evaluated=False, next_identity_authorized=False,
        automatic_retry=False, publication_ready=False)
    conversion_gate_path = conversion_paths[0] / "results.json"
    conversion_gate_ref = put(conversion_gate_path, conversion_gate)
    execution_ref = put(paths["cwd"] / "results.json", {"fixture": "original execution placeholder"})
    assessment_gate = dict(schema="conversion_gated_native_qfo_assessment_v1",
        status="native_qfo_assessment_process_succeeded_pending_independent_admission",
        source=sub["worker"], driver=sub["driver"], conversion_submission=conversion_ref,
        job_id="104", native_job_id=101, reviewer_job_id=102, conversion_job_id=103,
        index=8, cell="p0_c1_r0", output_namespaces=sub["output_namespaces"], execution=execution_ref,
        assessment_driver_invoked=True, endpoint_process_completed=True, accuracy_admitted=False,
        native_inference_reexecuted=False, automatic_retry=False, next_identity_authorized=False,
        publication_ready=False, conversion_gate=conversion_gate_ref, pairs=pairs_ref)
    assessment_gate_path = paths["gate"] / "results.json"
    put(assessment_gate_path, assessment_gate)
    fields = dict(JobIDRaw="104", State="COMPLETED", ExitCode="0:0", NodeList="bizon",
                  AllocCPUS="8", ReqMem="64G")
    monkeypatch.setattr(gate.admission, "accounting", lambda job, include_memory: ("fixture accounting", fields))
    monkeypatch.setattr(gate.assessment_gate.conversion_gate, "runtime_environment", lambda ref: {"fixture": True})
    context = dict(available_memory_bytes=64 * 2**30, available_disk_bytes=64 * 2**30)
    monkeypatch.setattr(gate.assessment_gate.conversion_gate, "capacity", lambda: context)
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
    monkeypatch.setenv("SLURM_JOB_ID", "105")
    calls = []
    output = gate.output_paths(run, 101)

    def admit(root, pairs_ref, conversion_job, assessment_job, destination):
        calls.append((root, pairs_ref, conversion_job, assessment_job, destination))
        report = dict(schema="full_native_factorial_qfo_admission_v1",
            status="full_native_factorial_qfo_assessment_admitted", source=record(gate.admission.__file__),
            pairs_manifest=pairs_ref, execution_report=execution_ref, accuracy_admitted=True,
            native_index=8, cell="p0_c1_r0", native_job_id=101,
            conversion=json.loads(Path(pairs_ref["path"]).read_text()),
            participant="ohmm_qfo_full_native_p0_c1_r0", publication_ready=False,
            automatic_retry=False, next_identity_authorized=False, scheduler=dict(fields))
        put(destination / "results.json", report)
        return report

    monkeypatch.setattr(gate.admission, "admit", admit)
    return dict(sub=sub, submission_path=submission_path, fields=fields, context=context, calls=calls,
        paths=output, assessment_gate=assessment_gate, assessment_gate_path=assessment_gate_path,
        conversion_gate=conversion_gate, conversion_gate_path=conversion_gate_path,
        execution_path=Path(execution_ref["path"]), pairs_path=Path(pairs_ref["path"]), admit=admit)


def execute(fixture):
    return gate.execute(record(fixture["submission_path"]), record(gate.__file__)["sha256"])


def result(fixture):
    return json.loads((fixture["paths"]["gate"] / "results.json").read_text())


def test_success_calls_original_validator_once_without_rerunning_endpoint_process(fixture):
    ref = execute(fixture)
    report = result(fixture)
    assert ref == record(fixture["paths"]["gate"] / "results.json")
    assert report["status"] == "assessment_gated_native_qfo_accuracy_admitted"
    assert report["validator_invoked"] is report["accuracy_admitted"] is True
    assert len(fixture["calls"]) == 1 and fixture["calls"][0][2:4] == ("103", "104")
    for key in ("native_inference_reexecuted", "automatic_retry", "next_identity_authorized", "publication_ready"):
        assert report[key] is False
    assert "successful terminal accounting" in " ".join(report["limitations"])


@pytest.mark.parametrize("key,value", (("JobIDRaw", "999"), ("State", "RUNNING"),
    ("State", "FAILED"), ("ExitCode", "1:0"), ("NodeList", "other"),
    ("AllocCPUS", "2"), ("ReqMem", "32G")))
def test_wrong_or_unfinished_assessment_never_reads_producer_outputs(fixture, key, value):
    fixture["fields"][key] = value
    fixture["assessment_gate_path"].unlink()
    fixture["execution_path"].unlink()
    with pytest.raises(ValueError, match="Assessment is not successfully completed"):
        execute(fixture)
    assert fixture["calls"] == [] and result(fixture)["validator_invoked"] is False


@pytest.mark.parametrize("key,value", (("status", "failed"), ("source", {}), ("driver", {}),
    ("conversion_submission", {}), ("job_id", "999"), ("native_job_id", 999),
    ("reviewer_job_id", 999), ("conversion_job_id", 999), ("index", 9), ("cell", "p0_c1_r1"),
    ("output_namespaces", {}), ("execution", {}), ("assessment_driver_invoked", False),
    ("endpoint_process_completed", False), ("accuracy_admitted", True),
    ("native_inference_reexecuted", True), ("automatic_retry", True),
    ("next_identity_authorized", True), ("publication_ready", True),
    ("conversion_gate", {}), ("pairs", {})))
def test_successful_scheduler_cannot_bypass_failed_or_mismatched_gate(fixture, key, value):
    fixture["assessment_gate"][key] = value
    put(fixture["assessment_gate_path"], fixture["assessment_gate"])
    with pytest.raises(ValueError):
        execute(fixture)
    assert fixture["calls"] == [] and result(fixture)["validator_invoked"] is False


def test_residual_conversion_stage_cannot_bypass_failed_conversion_gate(fixture):
    fixture["conversion_gate"]["status"] = "failed"
    put(fixture["conversion_gate_path"], fixture["conversion_gate"])
    with pytest.raises(ValueError, match="successful bound conversion gate"):
        execute(fixture)
    assert fixture["pairs_path"].exists() and fixture["calls"] == []


@pytest.mark.parametrize("key", ("available_memory_bytes", "available_disk_bytes"))
def test_unsafe_future_capacity_retains_refusal(fixture, key):
    fixture["context"][key] = 1
    with pytest.raises(ValueError, match="Unsafe admission"):
        execute(fixture)
    assert fixture["calls"] == [] and result(fixture)["accuracy_admitted"] is False


@pytest.mark.parametrize("key,value", (("job_id", True), ("job_id", 101),
    ("native_job_id", 101.0), ("reviewer_job_id", 102.0), ("conversion_job_id", 103.0),
    ("conversion_job_id", 999), ("index", 8.0), ("index", True), ("index", 5),
    ("cell", "p0_c1_r1"), ("request", {}), ("plan", {}), ("worker", {}), ("driver", {}),
    ("output_namespaces", {}), ("held_inspection_passed", False), ("accuracy_admitted", True),
    ("assessment_completed", True), ("automatic_retry", True)))
def test_invalid_submission_never_creates_admission_state(fixture, key, value):
    fixture["sub"][key] = value
    put(fixture["submission_path"], fixture["sub"])
    with pytest.raises(ValueError):
        execute(fixture)
    assert not fixture["paths"]["gate"].exists() and fixture["calls"] == []


@pytest.mark.parametrize("key,value", (("SLURM_CPUS_PER_TASK", "8"), ("SLURM_JOB_ID", ""),
    ("SLURM_JOB_ID", "101"), ("SLURM_JOB_ID", "102"), ("SLURM_JOB_ID", "103"), ("SLURM_JOB_ID", "104")))
def test_wrong_own_allocation_refuses_without_any_output(fixture, monkeypatch, key, value):
    monkeypatch.setenv(key, value)
    with pytest.raises(ValueError, match="separate scheduled"):
        execute(fixture)
    assert not fixture["paths"]["gate"].exists()


@pytest.mark.parametrize("role", ("gate", "admission"))
def test_occupied_namespaces_never_overwritten(fixture, role):
    fixture["paths"][role].mkdir(parents=True)
    marker = fixture["paths"][role] / "keep"
    marker.write_text("untouched")
    with pytest.raises(FileExistsError):
        execute(fixture)
    assert marker.read_text() == "untouched" and fixture["calls"] == []


@pytest.mark.parametrize("role", ("gate", "admission"))
def test_dangling_namespaces_never_followed(fixture, role):
    path = fixture["paths"][role]
    path.parent.mkdir(parents=True, exist_ok=True)
    target = path.parent / "absent"
    path.symlink_to(target)
    with pytest.raises(ValueError, match="direct admission namespaces"):
        execute(fixture)
    assert path.is_symlink() and not target.exists()


def test_original_validator_failure_is_retained_without_retry(fixture, monkeypatch):
    def fail(*args):
        fixture["calls"].append(args)
        raise ValueError("fixture independent validation failure")
    monkeypatch.setattr(gate.admission, "admit", fail)
    with pytest.raises(ValueError, match="fixture independent validation failure"):
        execute(fixture)
    assert len(fixture["calls"]) == 1
    assert result(fixture)["status"] == "assessment_gated_native_qfo_admission_failed_retained"
    assert result(fixture)["accuracy_admitted"] is result(fixture)["automatic_retry"] is False


@pytest.mark.parametrize("key,value", (("status", "failed"), ("source", {}), ("pairs_manifest", {}),
    ("execution_report", {}), ("accuracy_admitted", False), ("native_index", 9),
    ("native_job_id", 999), ("cell", "p0_c1_r1"), ("conversion", {}),
    ("participant", "other"), ("publication_ready", True), ("scheduler", {})))
def test_wrong_original_admission_is_retained_but_not_exportable(fixture, monkeypatch, key, value):
    def changed(*args):
        admitted = fixture["admit"](*args)
        admitted[key] = value
        put(fixture["paths"]["admission"] / "results.json", admitted)
        return admitted
    monkeypatch.setattr(gate.admission, "admit", changed)
    with pytest.raises(ValueError, match="Unexpected original independent admission"):
        execute(fixture)
    assert (fixture["paths"]["admission"] / "results.json").exists()
    assert result(fixture)["accuracy_admitted"] is False


def test_postflight_drift_preserves_original_admitted_report_but_refuses_wrapper(fixture, monkeypatch):
    def changed(*args):
        admitted = fixture["admit"](*args)
        put(fixture["submission_path"], {})
        return admitted
    monkeypatch.setattr(gate.admission, "admit", changed)
    with pytest.raises(ValueError):
        execute(fixture)
    original = json.loads((fixture["paths"]["admission"] / "results.json").read_text())
    assert original["accuracy_admitted"] is True and result(fixture)["accuracy_admitted"] is False
    assert "successful terminal accounting" in " ".join(result(fixture)["limitations"])


def test_changed_worker_refuses_before_creating_any_state(fixture):
    with pytest.raises(ValueError, match="Admission gate worker changed"):
        gate.execute(record(fixture["submission_path"]), "f" * 64)
    assert not fixture["paths"]["gate"].exists()


def test_changed_original_validator_refuses_before_independent_validation(fixture, monkeypatch):
    changed = fixture["submission_path"].parent / "changed_validator.py"
    changed.write_text("# fixture, never executed\n")
    monkeypatch.setattr(gate.admission, "__file__", str(changed))
    with pytest.raises(ValueError, match="Original independent validator changed"):
        execute(fixture)
    assert fixture["calls"] == [] and result(fixture)["validator_invoked"] is False


def test_wrong_runtime_refuses_before_reading_completed_producer_files(fixture, monkeypatch):
    def fail(ref):
        raise ValueError("fixture original runtime differs")
    monkeypatch.setattr(gate.assessment_gate.conversion_gate, "runtime_environment", fail)
    fixture["assessment_gate_path"].unlink()
    fixture["execution_path"].unlink()
    with pytest.raises(ValueError, match="fixture original runtime differs"):
        execute(fixture)
    assert fixture["calls"] == [] and result(fixture)["accuracy_admitted"] is False


@pytest.mark.parametrize("memory", ("64Gn", "65536M", "65536Mn"))
def test_equivalent_bound_scheduler_memory_units_are_supported(fixture, memory):
    fixture["fields"]["ReqMem"] = memory
    execute(fixture)
    assert result(fixture)["accuracy_admitted"] is True
