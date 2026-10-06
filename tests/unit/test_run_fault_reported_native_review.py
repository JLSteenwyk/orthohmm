"""Synthetic gates and child execution; no scheduler submission or raw inference."""

from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import run_fault_reported_native_review as module


def store(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))
    return module.record(path)


@pytest.fixture
def gate():
    request = dict(path="/fixture/request.json", bytes=1, sha256="a" * 64)
    plan = dict(path="/fixture/plan.json", bytes=1, sha256="b" * 64)
    output = dict(schema="native_factorial_output_review_v1", status="native_outputs_validated",
        job_id=22444, index=8, cell="p0_c1_r0", request=request, plan=plan, source=module.record(module.VALIDATOR),
        native_outputs_validated=True, terminal_scheduler_confirmed=True, accuracy_evaluated=False,
        terminal_reviewed=False, resource_measurements_admitted=False, next_identity_authorized=False,
        uncontended_timing=False)
    scheduler = dict(JobIDRaw="22734", State="COMPLETED", ExitCode="0:0", NodeList="bizon",
        AllocCPUS="2", ReqMem="32G")
    return output, scheduler, request, plan


def test_original_diagnostic_gate(gate):
    assert module.diagnostic_gate(*gate) is None


@pytest.mark.parametrize("side,key,value", [
    (0, "index", True), (0, "index", 9), (0, "cell", "p0_c1_r1"),
    (0, "job_id", 22445), (0, "accuracy_evaluated", True),
    (0, "terminal_reviewed", True), (0, "resource_measurements_admitted", True),
    (0, "next_identity_authorized", True), (0, "native_outputs_validated", False),
    (0, "source", {}), (0, "request", {}), (0, "plan", {}),
    (1, "JobIDRaw", "22445"), (1, "State", "RUNNING"),
    (1, "State", "FAILED"), (1, "ExitCode", "0:11"),
    (1, "AllocCPUS", "8"), (1, "ReqMem", "128G"), (1, "NodeList", "other"),
])
def test_diagnostic_cannot_substitute_for_full_admission(gate, side, key, value):
    changed = list(deepcopy(gate))
    changed[side][key] = value
    with pytest.raises(ValueError):
        module.diagnostic_gate(*changed)


def controller(job=22735, digest="d" * 64):
    user = module.pwd.getpwuid(module.os.getuid()).pw_name
    fields = dict(JobId=str(job), JobState="RUNNING", Partition="gpu", NodeList="bizon",
        NumNodes="1", NumCPUs="2", NumTasks="1", MinMemoryNode="128G", TimeLimit="06:00:00",
        Requeue="0", Restarts="0", Command=str(module.BATCH), WorkDir=str(module.ROOT), Comment=digest,
        UserId=f"{user}({module.os.getuid()})")
    fields["CPUs/Task"] = "2"
    return fields


def raw(fields):
    return " ".join(f"{k}={v}" for k, v in fields.items()) + "\n"


def test_declared_review_only_envelope():
    fields = controller()
    assert module.allocation_gate(raw(fields), 22735, "d" * 64) == fields


@pytest.mark.parametrize("key,value", [
    ("JobState", "PENDING"), ("MinMemoryNode", "32G"), ("NumCPUs", "64"),
    ("CPUs/Task", "8"), ("NodeList", "other"), ("TimeLimit", "1-02:00:00"),
    ("Requeue", "1"), ("Restarts", "1"), ("Command", "/wrong"),
    ("Comment", "e" * 64), ("ArrayJobId", "22735"), ("UserId", "other(2)"),
])
def test_allocation_refusals(key, value):
    fields = controller()
    fields[key] = value
    with pytest.raises(ValueError):
        module.allocation_gate(raw(fields), 22735, "d" * 64)


def test_duplicate_controller_field_refused():
    with pytest.raises(ValueError, match="Duplicate"):
        module.allocation_gate(raw(controller()).rstrip("\n") + " JobId=22735", 22735, "d" * 64)


def test_multiple_controller_records_refused():
    with pytest.raises(ValueError, match="one fresh"):
        module.allocation_gate(raw(controller()) * 2, 22735, "d" * 64)


@pytest.fixture
def attempt(tmp_path, monkeypatch):
    monkeypatch.setattr(module, "ROOT", tmp_path)
    paths = dict(REQUEST="request.json", REVIEWER="reviewer.py", VALIDATOR="validator.py",
        DIAGNOSTIC="diagnostic/outputs.json", PYTHON="python", BATCH="batch.sh",
        CONTROL="control", DESTINATION="review")
    for name, relative in paths.items():
        monkeypatch.setattr(module, name, tmp_path / relative)
    for name in ("REVIEWER", "VALIDATOR", "PYTHON", "BATCH"):
        getattr(module, name).write_text("synthetic " + name)
    monkeypatch.setattr(module, "REVIEWER_SHA", module.record(module.REVIEWER)["sha256"])
    monkeypatch.setattr(module, "VALIDATOR_SHA", module.record(module.VALIDATOR)["sha256"])
    plan = store(tmp_path / "plan.json", dict(helper_sources=[], evidence=[]))
    request = store(module.REQUEST, dict(plan=plan))
    monkeypatch.setattr(module, "REQUEST_SHA", request["sha256"])
    output = dict(schema="native_factorial_output_review_v1", status="native_outputs_validated",
        job_id=22444, index=8, cell="p0_c1_r0", request=request, plan=plan, source=module.record(module.VALIDATOR),
        native_outputs_validated=True, terminal_scheduler_confirmed=True, accuracy_evaluated=False,
        terminal_reviewed=False, resource_measurements_admitted=False, next_identity_authorized=False,
        uncontended_timing=False, checked_files=[], evidence=[])
    diagnostic = store(module.DIAGNOSTIC, output)
    producer = store(tmp_path / "producer.json", dict(schema="native_reviewer_signal11_diagnostic_submission_v1",
        job_id=22734, native_job_id=22444, original_failed_reviewer=22445, index=8,
        request=request, destination=str(module.DIAGNOSTIC.parent), automatic_retry=False,
        batch=module.record(module.BATCH)))
    monkeypatch.setattr(module, "DIAGNOSTIC_SUBMISSION", Path(producer["path"]))
    monkeypatch.setattr(module, "DIAGNOSTIC_SUBMISSION_SHA", producer["sha256"])
    scheduler = dict(JobIDRaw="22734", State="COMPLETED", ExitCode="0:0", NodeList="bizon",
        AllocCPUS="2", ReqMem="32G")
    monkeypatch.setattr(module, "accounting", lambda *a, **kw: ("synthetic", scheduler))
    monkeypatch.setattr(module, "validate_plan", lambda p: [{}] * 8 + [dict(cell="p0_c1_r0")])
    monkeypatch.setattr(module, "validate_request", lambda *a: None)
    monkeypatch.setattr(module, "runtime_environment", lambda *a: dict(synthetic=True))
    monkeypatch.setattr(module, "available_memory", lambda text: 2 * module.MEMORY)
    monkeypatch.setenv("SLURM_JOB_ID", "22735")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
    state = dict(mode="success", calls=[])

    def child(command, **kwargs):
        state["calls"].append(command)
        if command[0] == "sacct":
            return SimpleNamespace(stdout="22445|FAILED|0:11\n", stderr="", returncode=0)
        if command[0] == "scontrol":
            digest = module.record(module.__file__)["sha256"] if state["mode"] == "observed" else diagnostic["sha256"]
            return SimpleNamespace(stdout=raw(controller(digest=digest)), stderr="", returncode=0)
        assert command[:6] == [str(module.PYTHON), "-B", "-X", "faulthandler", "-m",
                              "benchmark_tools.review_native_factorial_attempt"]
        assert kwargs["cwd"] == tmp_path and kwargs["check"] is False
        if state["mode"] == "failed":
            return SimpleNamespace(returncode=-11)
        store(module.DESTINATION / "review.json", dict(schema="native_factorial_terminal_review_v1",
            status="native_success", job_id=22444, index=8, request=request,
            source=module.record(module.REVIEWER), native_outputs_validated=True,
            terminal_reviewed=True, next_identity_authorized=True))
        if state["mode"] == "mutated":
            module.VALIDATOR.write_text("mutated synthetic validator")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(module.subprocess, "run", child)
    return diagnostic, state, scheduler


def run(attempt):
    return module.execute(attempt[0]["sha256"], module.record(module.__file__)["sha256"])


def test_executes_original_full_reviewer_and_keeps_scope(attempt):
    report = module.read(run(attempt))
    assert len(attempt[1]["calls"]) == 3
    assert report["status"] == "original_full_review_returned_pending_independent_completion_check"
    assert report["original_full_review_reexecuted"] is True
    for key in ("native_inference_reexecuted", "accuracy_evaluated", "next_identity_authorized",
                "scientific_timings_admitted", "automatic_retry", "publication_ready"):
        assert report[key] is False
    assert (module.CONTROL / "preflight.json").exists()
    assert (module.DESTINATION / "review.json").exists()


def test_failed_child_retained_without_retry(attempt):
    attempt[1]["mode"] = "failed"
    with pytest.raises(ValueError, match="failed; retain"):
        run(attempt)
    result = json.loads((module.CONTROL / "results.json").read_text())
    assert result["status"] == "fault_reported_full_review_failed_retained"
    assert result["child_exit_code"] == -11 and result["automatic_retry"] is False
    assert len(attempt[1]["calls"]) == 3


def test_postflight_source_mutation_refused(attempt):
    attempt[1]["mode"] = "mutated"
    with pytest.raises(ValueError):
        run(attempt)
    assert json.loads((module.CONTROL / "results.json").read_text())["status"] == "fault_reported_full_review_failed_retained"


def test_pending_diagnostic_never_executes_review(attempt):
    attempt[2]["State"] = "RUNNING"
    with pytest.raises(ValueError, match="successful original diagnostic"):
        run(attempt)
    assert not attempt[1]["calls"] and not module.CONTROL.exists()


def test_existing_namespace_never_overwritten(attempt):
    module.DESTINATION.mkdir()
    with pytest.raises(ValueError, match="fresh direct"):
        run(attempt)
    assert len(attempt[1]["calls"]) == 2 and not module.CONTROL.exists()


def test_wrong_worker_digest_never_queries_scheduler(attempt):
    with pytest.raises(ValueError, match="wrapper source changed"):
        module.execute(attempt[0]["sha256"], "0" * 64)
    assert not attempt[1]["calls"]


def test_future_digest_observed_only_after_successful_producer(attempt):
    attempt[1]["mode"] = "observed"
    report = module.read(module.execute(None, module.record(module.__file__)["sha256"]))
    assert report["diagnostic"] == attempt[0]
    assert report["diagnostic_digest_mode"] == "observed_after_producer_completion"


def test_pending_producer_blocks_before_opening_output(attempt, monkeypatch):
    attempt[2]["State"] = "RUNNING"
    real_record = module.record
    opened = []

    def observed(path):
        opened.append(Path(path))
        return real_record(path)

    monkeypatch.setattr(module, "record", observed)
    with pytest.raises(ValueError, match="successful original diagnostic"):
        module.execute(None, real_record(module.__file__)["sha256"])
    assert module.DIAGNOSTIC not in opened and not attempt[1]["calls"]


def test_wrong_producer_binding_blocks_output_and_execution(attempt, monkeypatch):
    monkeypatch.setattr(module, "DIAGNOSTIC_SUBMISSION_SHA", "0" * 64)
    with pytest.raises(ValueError, match="submission changed"):
        module.execute(None, module.record(module.__file__)["sha256"])
    assert not attempt[1]["calls"] and not module.CONTROL.exists()
