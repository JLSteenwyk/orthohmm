"""Gates for one prospective child full review; no scientific execution in tests."""

import ast
import copy
import os
from pathlib import Path
import pwd
import subprocess

import pytest

from benchmark_tools import run_native11_fault_reported_review as worker


DIGEST = "a" * 64
REF = dict(path="/fixture/request.json", bytes=1, sha256=DIGEST)
PLAN = dict(path="/fixture/plan.json", bytes=1, sha256="b" * 64)
AMENDMENT = dict(path="/fixture/amendment.json", bytes=1, sha256="c" * 64)
VALIDATOR = dict(path="/fixture/validator.py", bytes=1, sha256="d" * 64)


def scheduler():
    return dict(JobIDRaw="24031", State="COMPLETED", ExitCode="0:0", NodeList="bizon",
                AllocCPUS="2", ReqMem="32G")


def diagnostic():
    return dict(schema="allocated_native_factorial_output_review_v1", status="native_outputs_validated",
        job_id=23985, index=11, cell="p1_c1_r0", request=REF, plan=PLAN, amendment=AMENDMENT,
        source=VALIDATOR, native_outputs_validated=True, terminal_scheduler_confirmed=True,
        accuracy_evaluated=False, terminal_reviewed=False, resource_measurements_admitted=False,
        next_identity_authorized=False, uncontended_timing=False)


def controller():
    return dict(JobId="25000", JobName="ohmm_native11_fullreview", JobState="RUNNING",
        Partition="gpu", NodeList="bizon", NumNodes="1", NumCPUs="2", NumTasks="1",
        MinMemoryNode="128G", TimeLimit="06:00:00", Requeue="0", Restarts="0",
        Command=str(worker.BATCH), WorkDir=str(worker.ROOT), Comment=DIGEST,
        UserId=f"{pwd.getpwuid(os.getuid()).pw_name}({os.getuid()})", **{"CPUs/Task": "2"})


def raw(fields):
    return " ".join(f"{k}={v}" for k, v in fields.items()) + "\n"


def test_bound_diagnostic_only_result_passes():
    worker.diagnostic_gate(diagnostic(), scheduler(), REF, PLAN, AMENDMENT, VALIDATOR)
    assert worker.allocation_gate(raw(controller()), 25000, DIGEST) == controller()


@pytest.mark.parametrize("key,value", [
    ("JobIDRaw", "23986"), ("State", "RUNNING"), ("State", "FAILED"),
    ("ExitCode", "0:11"), ("NodeList", "another-host"), ("AllocCPUS", "1"), ("ReqMem", "128G"),
])
def test_diagnostic_must_be_successfully_terminal_in_declared_envelope(key, value):
    fields = scheduler()
    fields[key] = value
    with pytest.raises(ValueError, match="successful diagnostic24031"):
        worker.diagnostic_gate(diagnostic(), fields, REF, PLAN, AMENDMENT, VALIDATOR)


@pytest.mark.parametrize("key,value", [
    ("schema", "native_factorial_output_review_v1"), ("status", "failed"),
    ("job_id", 22444), ("index", True), ("index", 10), ("cell", "p0_c1_r0"),
    ("request", PLAN), ("plan", REF), ("amendment", REF), ("source", REF),
    ("native_outputs_validated", False), ("terminal_scheduler_confirmed", False),
    ("accuracy_evaluated", True), ("terminal_reviewed", True),
    ("resource_measurements_admitted", True), ("next_identity_authorized", True), ("uncontended_timing", True),
])
def test_wrong_identity_or_diagnostic_promoted_to_admission_refused(key, value):
    output = copy.deepcopy(diagnostic())
    output[key] = value
    with pytest.raises(ValueError, match="bound standalone semantic result"):
        worker.diagnostic_gate(output, scheduler(), REF, PLAN, AMENDMENT, VALIDATOR)


@pytest.mark.parametrize("key,value", [
    ("JobState", "PENDING"), ("JobId", "24031"), ("JobName", "different"),
    ("NodeList", "other"), ("NumNodes", "2"), ("NumCPUs", "4"),
    ("CPUs/Task", "1"), ("MinMemoryNode", "32G"), ("TimeLimit", "01:00:00"),
    ("Requeue", "1"), ("Restarts", "1"), ("Command", "/another/script.sh"),
    ("WorkDir", "/tmp"), ("Comment", "0" * 64), ("UserId", "other(999)"), ("ArrayJobId", "25000"),
])
def test_owned_fresh_two_cpu_full_review_allocation_required(key, value):
    fields = controller()
    fields[key] = value
    with pytest.raises(ValueError, match="distinct owned"):
        worker.allocation_gate(raw(fields), 25000, DIGEST)


def test_duplicate_or_multiple_controller_records_refused():
    with pytest.raises(ValueError, match="Duplicate"):
        worker.allocation_gate(raw(controller()).strip() + " JobId=25000\n", 25000, DIGEST)
    with pytest.raises(ValueError, match="one fresh"):
        worker.allocation_gate(raw(controller()) * 2, 25000, DIGEST)


def test_failed_accounting_precedes_future_diagnostic_output_access(monkeypatch):
    seen = []
    pins = {worker.REQUEST: worker.REQUEST_SHA, worker.REVIEWER: worker.REVIEWER_SHA,
        worker.VALIDATOR: worker.VALIDATOR_SHA, worker.SEMANTIC: worker.SEMANTIC_SHA,
        worker.PREPARATION: worker.PREPARATION_SHA, worker.SUBMISSION: worker.SUBMISSION_SHA}

    def observe(path):
        path = Path(path)
        seen.append(path)
        assert path != worker.DIAGNOSTIC
        return dict(path=str(path), bytes=1, sha256=pins.get(path, DIGEST))

    def fail_accounting(*args, **kwargs):
        raise ValueError("producer has not completed")

    monkeypatch.setattr(worker, "record", observe)
    monkeypatch.setattr(worker, "read", lambda ref: dict(amendment=AMENDMENT, plan=PLAN))
    monkeypatch.setattr(worker, "amendment", lambda ref: ({}, {"runs": [{}] * 11 + [{"cell": "p1_c1_r0"}]}))
    monkeypatch.setattr(worker, "validate_request", lambda *args: None)
    monkeypatch.setattr(worker, "accounting", fail_accounting)
    with pytest.raises(ValueError, match="producer has not completed"):
        worker.execute(DIGEST, DIGEST)
    assert worker.DIAGNOSTIC not in seen


def test_batch_syntax_and_original_child_cli_without_alterations():
    result = subprocess.run(["bash", "-n", str(worker.BATCH)], capture_output=True, text=True, timeout=5)
    assert result.returncode == 0
    source = Path(worker.__file__).read_text()
    assert '"-m", "benchmark_tools.review_allocated_native_factorial_attempt"' in source
    calls = {ast.unparse(node.func) for node in ast.walk(ast.parse(source)) if isinstance(node, ast.Call)}
    assert not calls.intersection({"gc.disable", "gc.set_threshold", "setattr", "monkeypatch.setattr", "validate_semantics"})
    batch = worker.BATCH.read_text()
    for directive in ("--cpus-per-task=2", "--mem=128G", "--time=06:00:00", "--no-requeue"):
        assert "#SBATCH " + directive in batch
    assert "--diagnostic-sha256" in batch and "--worker-sha256" in batch
    assert "native_factorial_review_py310_20261004/bin/python" in batch


@pytest.mark.parametrize("args", [[], [DIGEST], ["x", DIGEST], [DIGEST, "F" * 64], [DIGEST, DIGEST, "extra"]])
def test_batch_refuses_invalid_digest_arguments_before_execution(args):
    env = {k: v for k, v in os.environ.items() if not k.startswith("SLURM_")}
    result = subprocess.run(["bash", str(worker.BATCH), *args], env=env, capture_output=True, text=True, timeout=5)
    assert result.returncode == 2 and not result.stdout
    assert "Require verified diagnostic" in result.stderr


@pytest.mark.parametrize("job", ["", "0", "23985", "23986", "24031", "024032"])
def test_batch_refuses_unscheduled_or_reused_job_before_execution(job):
    env = {k: v for k, v in os.environ.items() if not k.startswith("SLURM_")}
    env.update(SLURM_CPUS_PER_TASK="2", SLURM_JOB_ID=job)
    result = subprocess.run(["bash", str(worker.BATCH), DIGEST, DIGEST], env=env,
                            capture_output=True, text=True, timeout=5)
    assert result.returncode == 2 and not result.stdout
    assert "Require a distinct scheduled" in result.stderr
