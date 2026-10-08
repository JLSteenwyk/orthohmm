import copy
import os
import pwd

import pytest

from benchmark_tools import review_native11_composed_attempt as current
from benchmark_tools.prepare_allocated_native_factorial_qfo_pairs import admit_conversion


def controller():
    fields = dict(JobId="25000", JobName="ohmm_native11_composed", JobState="RUNNING",
        Partition="gpu", NodeList="bizon", NumNodes="1", NumCPUs="2", NumTasks="1",
        MinMemoryNode="128G", TimeLimit="06:00:00", Requeue="0", Restarts="0",
        Command=str(current.BATCH), WorkDir=str(current.ROOT), Comment="digest",
        UserId=f"{pwd.getpwuid(os.getuid()).pw_name}({os.getuid()})")
    fields["CPUs/Task"] = "2"
    return fields


def test_exact_owned_scheduled_review_allocation_is_required():
    fields = controller()
    raw = " ".join(f"{k}={v}" for k, v in fields.items())
    assert current.allocation_gate(raw, 25000, "digest") == fields


@pytest.mark.parametrize("field", ["JobId", "JobName", "NumCPUs", "CPUs/Task", "MinMemoryNode",
    "TimeLimit", "Comment", "Requeue", "Restarts", "UserId", "Command", "WorkDir"])
def test_wrong_review_allocation_is_rejected(field):
    fields = controller()
    fields[field] = "wrong"
    with pytest.raises(ValueError):
        current.allocation_gate(" ".join(f"{k}={v}" for k, v in fields.items()), 25000, "digest")


@pytest.mark.parametrize("extra", ["JobId=25000", "ArrayJobId=25000", "HetJobId=25000"])
def test_duplicate_array_and_heterogeneous_allocations_are_rejected(extra):
    fields = controller()
    raw = " ".join(f"{k}={v}" for k, v in fields.items()) + " " + extra
    with pytest.raises(ValueError):
        current.allocation_gate(raw, 25000, "digest")


def runtime_fixture():
    request = dict(plan={"plan": 1}, amendment={"amendment": 1})
    report = dict(schema=current.runtime_v2.SCHEMA,
        status="historical_brackets_and_current_exact_additions_revalidated",
        source=current.record(current.runtime_v2.__file__),
        runtime_kernel_source=current.record(current.runtime_kernel.__file__),
        request={"request": 1}, plan=request["plan"], amendment=request["amendment"],
        job_id=23985, index=11, cell="p1_c1_r0", historical_first_terminal_rechecked=True,
        current_original_entries_unchanged=True, current_private_inventory_equality=True,
        current_original_inventory_equality=False, continuous_runtime_integrity_established=False,
        full_review_admitted=False, terminal_reviewed=False, next_identity_authorized=False,
        accuracy_evaluated=False, automatic_retry=False, publication_ready=False)
    return report, request


def test_truthful_nonadmitting_runtime_component_gate():
    report, request = runtime_fixture()
    current.runtime_gate(report, {"request": 1}, request)


@pytest.mark.parametrize("field", ["schema", "source", "request", "plan", "amendment", "job_id",
    "cell", "historical_first_terminal_rechecked", "current_original_entries_unchanged",
    "current_private_inventory_equality", "current_original_inventory_equality",
    "continuous_runtime_integrity_established", "terminal_reviewed", "next_identity_authorized"])
def test_runtime_gate_rejects_changed_identity_or_false_current_claim(field):
    report, request = runtime_fixture()
    report[field] = not report[field] if type(report[field]) is bool else "wrong"
    with pytest.raises(ValueError):
        current.runtime_gate(report, {"request": 1}, request)


@pytest.fixture
def output_fixture(monkeypatch, tmp_path):
    def record(path):
        path = str(path)
        digests = {str(current.OUTPUT): current.OUTPUT_SHA, str(current.READBACK): current.READBACK_SHA}
        return dict(path=path, bytes=1, sha256=digests.get(path, "fixture"))
    monkeypatch.setattr(current, "record", record)
    request_ref = dict(request=1)
    request = dict(plan=dict(plan=1), amendment=dict(amendment=1))
    execution = dict(historical_plan=request["plan"])
    run = dict(output_root=str(tmp_path))
    output = dict(schema="allocated_native_factorial_output_review_v1", status="native_outputs_validated",
        job_id=23985, index=11, cell="p1_c1_r0", request=request_ref,
        plan=request["plan"], amendment=request["amendment"],
        source=record(current.ROOT / "benchmark_tools/validate_allocated_native_factorial_outputs.py"),
        semantic_validator_source=record(current.ROOT / "benchmark_tools/validate_native_factorial_outputs.py"),
        native_outputs_validated=True, terminal_scheduler_confirmed=True, execution_scope=current.SCOPE,
        allocated_ready=record(tmp_path / "measurement/ready.json"), checked_files=[record("checked")],
        evidence=[record("evidence")], accuracy_evaluated=False, terminal_reviewed=False,
        resource_measurements_admitted=False, next_identity_authorized=False, uncontended_timing=False)
    readback = dict(status="standalone_semantic_diagnostic_completed_bound", output=record(current.OUTPUT),
        job_id=24031, standalone_semantics_passed=True, full_review_admitted=False,
        submission=record("submission"), release=record("release"),
        direct_invocation_files=dict(time=record("time")), independent_checked_records=[record("independent")])
    monkeypatch.setattr(current, "read", lambda ref: output if ref["path"] == str(current.OUTPUT) else readback)
    fields = dict(JobIDRaw="24031", State="COMPLETED", ExitCode="0:0", NodeList="bizon", AllocCPUS="2", ReqMem="32G")
    monkeypatch.setattr(current, "accounting", lambda *a, **k: ("raw", fields))
    checked = []
    monkeypatch.setattr(current, "check", lambda ref: checked.append(ref))
    class Evidence:
        def bind(self, path):
            pass
    return dict(request_ref=request_ref, request=request, execution=execution, run=run,
        evidence=Evidence()), output, readback, fields, checked


def test_output_component_reuses_original_producer_gate_and_rebinds_all_files(output_fixture):
    arguments, output, readback, fields, checked = output_fixture
    ref, observed, producer = current.output_binding(**arguments)
    assert observed == output and ref["sha256"] == current.OUTPUT_SHA
    assert producer["producer_job_id"] == 24031
    assert producer["semantic_validator_reexecuted"] is False
    assert {"checked", "evidence", "independent", "submission", "release", "time"} <= {r["path"] for r in checked}


@pytest.mark.parametrize("change", ["output_source", "semantic_source", "placement", "producer",
                                   "readback", "readback_admission", "output_admission"])
def test_invalid_standalone_component_cannot_be_adopted(output_fixture, change):
    arguments, output, readback, fields, checked = output_fixture
    if change == "output_source":
        output["source"] = {}
    elif change == "semantic_source":
        output["semantic_validator_source"] = {}
    elif change == "placement":
        output["allocated_ready"] = {}
    elif change == "producer":
        fields["State"] = "FAILED"
    elif change == "readback":
        readback["output"] = {}
    elif change == "readback_admission":
        readback["full_review_admitted"] = True
    else:
        output["terminal_reviewed"] = True
    with pytest.raises(ValueError):
        current.output_binding(**arguments)


def test_new_composed_type_is_not_an_original_ordinary_review():
    report = dict(schema=current.SCHEMA, terminal_reviewed=True, composed_full_review_complete=True,
        next_identity_authorized=False)
    with pytest.raises(ValueError):
        admit_conversion(report, {}, dict(plan={}, amendment={}, job_id=23985),
            dict(index=11, dataset="qfo_corrected", cell="p1_c1_r0", repeat=0))


def test_unscheduled_invocation_refuses_before_output_creation(monkeypatch, tmp_path):
    destination = tmp_path / "composed"
    monkeypatch.setattr(current, "DESTINATION", destination)
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
    ref = current.record(current.__file__)
    with pytest.raises(ValueError, match="scheduled"):
        current.review(ref["sha256"])
    assert not destination.exists()
