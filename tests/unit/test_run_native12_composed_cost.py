"""Prospective controller fixtures, not actual inference or launch evidence."""

from contextlib import nullcontext
import json
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest

from benchmark_tools import run_native12_composed_cost as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def store(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding="ascii")
    return record(path)


def held_raw(job=25000, **overrides):
    fields = dict(JobId=str(job), JobName=module.JOB_NAME, JobState="PENDING", Reason="JobHeldUser",
        Partition="gpu", ReqNodeList="bizon", NumCPUs="64", NumTasks="1", MinMemoryNode="128G",
        Requeue="0", Restarts="0", TimeLimit="1-02:00:00", Command=str(module.SCRIPT),
        WorkDir=str(module.ROOT), NumNodes="1", UserId=f"fixture({os.getuid()})")
    fields["CPUs/Task"] = "64"
    fields.update(overrides)
    return " ".join(f"{key}={value}" for key, value in fields.items())


def test_held_scheduler_gate_binds_new_batch_and_matched_resources():
    assert module.held_gate(held_raw(), 25000)["Command"] == str(module.SCRIPT)


@pytest.mark.parametrize("field,value", [
    ("JobName", "other"), ("JobState", "RUNNING"), ("Reason", "Resources"),
    ("Partition", "other"), ("ReqNodeList", "other"), ("NumCPUs", "32"),
    ("CPUs/Task", "32"), ("MinMemoryNode", "32G"), ("NumTasks", "2"),
    ("Requeue", "1"), ("Restarts", "1"), ("TimeLimit", "06:00:00"),
    ("Command", "/old/batch.sh"), ("WorkDir", "/other"), ("NumNodes", "2"),
    ("UserId", "other(999999)"), ("ArrayJobId", "25000"), ("HetJobId", "25000"),
])
def test_wrong_held_identity_rejected(field, value):
    with pytest.raises(ValueError):
        module.held_gate(held_raw(**{field: value}), 25000)


def test_duplicate_fields_and_old_job_rejected():
    with pytest.raises(ValueError):
        module.held_gate(held_raw() + " JobId=25000", 25000)
    with pytest.raises(ValueError):
        module.held_gate(held_raw(job=24034), 24034)


@pytest.fixture
def fixture(tmp_path, monkeypatch):
    monkeypatch.setattr(module, "REQUEST", tmp_path / "request.json")
    core = tmp_path / "core"
    core.mkdir()
    baseline_ref = store(tmp_path / "baseline.json", dict(core_root=str(core), core_commit=module.CORE_COMMIT,
        environment_overrides={}, tool_entrypoints=dict(orthohmm_python=dict(absolute_path="/fixture/native/python"))))
    controller_ref = record(sys.executable)
    binding_ref = store(tmp_path / "lookup_binding.json", dict(controller_python=controller_ref,
        runtime_specs=[["/fixture/os.json", "os"], ["/fixture/private.json", "private"]]))
    lookup_ref = store(tmp_path / "lookup.json", dict(binding=binding_ref, baseline=baseline_ref))
    run = dict(index=12, cell="p1_c1_r1", dataset="qfo_corrected", inputs=[],
               output_root=str(tmp_path / "run_12"))
    plan = dict(runs=[{} for _ in range(12)] + [run], panel_root=str(tmp_path / "panel"),
                runtime_lookup=lookup_ref, baseline=baseline_ref)
    plan_ref = store(tmp_path / "plan.json", plan)
    policy_ref = store(tmp_path / "policy.json", dict(plan_sha256=plan_ref["sha256"],
        minimum_available_memory_bytes=module.MEMORY))
    amendment_ref = store(tmp_path / "amendment.json", dict(historical_plan=plan_ref))
    prior_ref = store(tmp_path / "prior.json", dict(prior=True))
    original = dict(plan=plan_ref, amendment=amendment_ref, policy=policy_ref, history=[prior_ref] * 11)
    original_ref = store(tmp_path / "original_request.json", original)
    runtime_ref = store(tmp_path / "runtime.json", dict(runtime=True))
    review = dict(reviews=dict(runtime=runtime_ref))
    review_ref = store(tmp_path / "composed_review.json", review)
    context = (original, dict(historical_plan=plan_ref), plan, {}, review, {}, "group", {}, [], {})
    history = dict(schema=module.HISTORY_SCHEMA, status="prefix_revalidated_for_native12_preparation",
        original_request=original_ref, composed_review=review_ref, historical_prefix=original["history"],
        next_unrun_index=12, original_review_translated=False, historical_failures_retained=[23986, 24033],
        next_identity_authorized=False, evidence=[])
    request = dict(schema=module.REQUEST_SCHEMA, job_id=25000, index=12, cell="p1_c1_r1",
        execution_authorized=True, plan=plan_ref, amendment=amendment_ref, policy=policy_ref,
        scheduler_command=str(module.SCRIPT), allocation_cwd=str(module.ROOT),
        automatic_retry=False, scientific_method_modified=False, original_review_translated=False,
        scientific_timings_admitted=False, publication_ready=False,
        source=record(module.__file__), history=original["history"] + [review_ref],
        original_request=original_ref, composed_review=review_ref,
        runtime_basis=runtime_ref, new_sources=module.sources(), history_basis=history)
    return SimpleNamespace(context=context, history=history, request=request, run=run,
        original_ref=original_ref, review_ref=review_ref, core=core, tmp=tmp_path)


def test_request_gate_preserves_new_schema_and_original_history(fixture):
    module.request_gate(fixture.request, 25000, fixture.context)
    assert fixture.request["schema"] != "allocated_native_factorial_request_v1"
    assert fixture.request["history"][:-1] == fixture.context[0]["history"]


@pytest.mark.parametrize("key,value", [
    ("schema", "allocated_native_factorial_request_v1"), ("job_id", 24034),
    ("index", 11), ("index", True), ("cell", "p1_c1_r0"), ("execution_authorized", False),
    ("plan", {}), ("amendment", {}), ("policy", {}), ("scheduler_command", "/old/batch.sh"),
    ("allocation_cwd", "/other"), ("automatic_retry", True), ("scientific_method_modified", True),
    ("original_review_translated", True), ("scientific_timings_admitted", True),
    ("publication_ready", True), ("source", {}), ("history", []), ("runtime_basis", {}),
    ("new_sources", []), ("history_basis", {}),
])
def test_wrong_request_rejected(fixture, key, value):
    fixture.request[key] = value
    with pytest.raises(ValueError):
        module.request_gate(fixture.request, 25000, fixture.context)


@pytest.mark.parametrize("field,value", [
    ("schema", "old"), ("status", "unreviewed"), ("original_request", {}),
    ("composed_review", {}), ("historical_prefix", []), ("next_unrun_index", 11),
    ("original_review_translated", True), ("historical_failures_retained", []),
    ("next_identity_authorized", True),
])
def test_wrong_history_basis_rejected(fixture, field, value):
    fixture.request["history_basis"][field] = value
    with pytest.raises(ValueError):
        module.request_gate(fixture.request, 25000, fixture.context)


def prepare_stubs(monkeypatch, fixture):
    calls = []
    monkeypatch.setattr(module, "history_binding",
        lambda a, b: calls.append((a, b)) or (fixture.context, fixture.history))
    monkeypatch.setattr(module, "available_memory", lambda _: module.MEMORY)
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k:
        SimpleNamespace(stdout=held_raw(), stderr="", returncode=0))
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: "fixture_commit\n")
    return calls


def test_prepare_creates_one_bound_request_not_native_output(fixture, monkeypatch):
    calls = prepare_stubs(monkeypatch, fixture)
    ref = module.prepare(25000, fixture.original_ref, fixture.review_ref, record(module.__file__)["sha256"])
    request = module.read(ref)
    module.request_gate(request, 25000, fixture.context)
    assert calls == [(fixture.original_ref, fixture.review_ref)]
    assert request["capacity_precheck"]["background_cpu_used_for_eligibility"] is False
    assert not Path(fixture.run["output_root"]).exists()
    assert not (fixture.tmp / "panel/sessions/run_12").exists()
    with pytest.raises(ValueError):
        module.prepare(25000, fixture.original_ref, fixture.review_ref, record(module.__file__)["sha256"])


@pytest.mark.parametrize("failure", ["source", "memory", "held", "prior_output"])
def test_prepare_rejects_unsafe_or_repeated_attempt_without_request(fixture, monkeypatch, failure):
    prepare_stubs(monkeypatch, fixture)
    digest = record(module.__file__)["sha256"]
    if failure == "source":
        digest = "0" * 64
    elif failure == "memory":
        monkeypatch.setattr(module, "available_memory", lambda _: module.MEMORY - 1)
    elif failure == "held":
        monkeypatch.setattr(module.subprocess, "run", lambda *a, **k:
            SimpleNamespace(stdout=held_raw(Command="/old/batch.sh"), stderr="", returncode=0))
    else:
        Path(fixture.run["output_root"]).mkdir()
    with pytest.raises(ValueError):
        module.prepare(25000, fixture.original_ref, fixture.review_ref, digest)
    assert not module.REQUEST.exists()


def test_execution_binding_requires_actual_consumers_and_source_pins(fixture, monkeypatch):
    ref = store(module.REQUEST, fixture.request)
    calls = []
    monkeypatch.setattr(module, "history_binding", lambda a, b:
        calls.append((a, b)) or (fixture.context, fixture.history))
    request, context, _ = module.execution_binding(ref, 25000)
    assert request == fixture.request and context is fixture.context
    assert calls == [(fixture.original_ref, fixture.review_ref)]
    wrong_ref = store(fixture.tmp / "other_request.json", fixture.request)
    with pytest.raises(ValueError):
        module.execution_binding(wrong_ref, 25000)


def execution_stubs(fixture, monkeypatch, phase_failure=None):
    from benchmark_tools import isolated_numba_cache, isolated_native_tmp
    from benchmark_tools import measure_threadripper_run, measure_allocated_threadripper_scaling
    from benchmark_tools import review_threadripper_process_stream, run_simulation_methods

    ref = store(module.REQUEST, fixture.request)
    calls = []
    monkeypatch.setattr(module, "execution_binding", lambda a, b:
        (fixture.request, fixture.context, fixture.history))
    monkeypatch.setenv("SLURM_JOB_ID", "25000")
    monkeypatch.setenv("PYTHONHASHSEED", "0")
    monkeypatch.setenv("PYTHONNOUSERSITE", "1")
    monkeypatch.setattr(module.sys, "dont_write_bytecode", True)
    monkeypatch.setattr(module.sys, "pycache_prefix", str(fixture.tmp / "absent_python_cache"))
    monkeypatch.setattr(module.sys, "executable", str(Path(sys.executable).resolve()))
    monkeypatch.chdir(module.ROOT)
    fields = dict(Comment=ref["sha256"], JobName=module.JOB_NAME, UserId=f"fixture({os.getuid()})")
    monkeypatch.setattr(module, "controller", lambda *a, **k: dict(fields=fields))
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k:
        SimpleNamespace(stdout="stub controller", stderr="", returncode=0))
    monkeypatch.setattr(module, "amendment", lambda _: (fixture.context[1], fixture.context[2]))
    monkeypatch.setattr(module, "original_inputs", lambda _: calls.append("original_inputs"))
    monkeypatch.setattr(module, "prepared_inputs", lambda *a: dict(copied=True))
    monkeypatch.setattr(measure_threadripper_run, "assert_environment", lambda *a: None)
    monkeypatch.setattr(run_simulation_methods, "verify_environment", lambda *a: None)
    monkeypatch.setattr(run_simulation_methods, "execution_environment", lambda *a: ({}, {}))
    monkeypatch.setattr(isolated_numba_cache, "fresh_cache", lambda *a: nullcontext())
    monkeypatch.setattr(isolated_native_tmp, "fresh_tmp", lambda *a: nullcontext())
    monkeypatch.setattr(review_threadripper_process_stream, "shared_environment", lambda _: True)
    monkeypatch.setattr(review_threadripper_process_stream, "audit", lambda *a, **k:
        (store(Path(a[0]) / "process_stream_review.json", dict(stub=True)),
         dict(sampled_environment_policy_satisfied=True)))

    class RuntimeStub:
        def __init__(self, plan, runtime_ref, output):
            assert plan is fixture.context[2] and runtime_ref == fixture.request["runtime_basis"]
            output.mkdir()
            self.count = 0

        def __call__(self, specs):
            self.count += 1
            calls.append(f"runtime_{self.count}")
            if phase_failure == self.count:
                raise ValueError("stub runtime drift")
            return dict(stub_runtime=True)

    monkeypatch.setattr(module, "ProspectiveRuntimeChecker", RuntimeStub)

    def prepare_inputs(run, baseline):
        calls.append("prepare")
        (Path(run["output_root"]) / "input").mkdir()

    def measure(command, directory, job, cpus, memory, timeout, interval, **kwargs):
        calls.append(("measure", command, job, cpus, memory, timeout, interval, kwargs))
        directory.mkdir()
        store(directory / "done.json", dict(stub_done=True))
        store(directory / "environment_preflight.json", dict(stub_preflight=True))
        return dict(status="command_exited_zero", native=dict(exit_code=0))

    monkeypatch.setattr(module, "prepare_inputs", prepare_inputs)
    monkeypatch.setattr(measure_allocated_threadripper_scaling, "measure", measure)
    return ref, calls, fields


def test_controller_reuses_original_native_command_and_measurement_envelope(fixture, monkeypatch):
    ref, calls, _ = execution_stubs(fixture, monkeypatch)
    before = Path.cwd()
    result = module.execute(ref)
    assert Path.cwd() == before
    assert result["schema"] == module.SESSION_SCHEMA
    assert result["status"] == "measurement_returned_pending_independent_review"
    assert result["next_identity_authorized"] is False
    assert result["sampled_environment_evidence_valid"] is True
    measured = next(call for call in calls if isinstance(call, tuple))
    assert measured[1] == module.native_command(fixture.request["amendment"], fixture.run,
        module.read(fixture.context[2]["baseline"]))
    assert measured[2:7] == (25000, 32, module.MEMORY, 85800, 1.)
    assert measured[7]["monitor_host"] is True and measured[7]["host_interval_s"] == 30.
    assert isinstance(measured[7]["release_guard"], module.ReleaseGuard)
    assert calls.index("runtime_1") < calls.index("prepare") < calls.index(measured) < calls.index("runtime_2")
    retained = json.loads((fixture.tmp / "panel/sessions/run_12/result.json").read_text())
    assert retained == result
    assert (Path(fixture.run["output_root"]) / "verification.json").is_file()


@pytest.mark.parametrize("phase", [1, 2])
def test_runtime_failure_prevents_or_preserves_native_measurement(fixture, monkeypatch, phase):
    ref, calls, _ = execution_stubs(fixture, monkeypatch, phase_failure=phase)
    result = module.execute(ref)
    assert result["status"] == "factorial_attempt_failed_retained"
    measured = [call for call in calls if isinstance(call, tuple)]
    assert len(measured) == (0 if phase == 1 else 1)
    assert result["wrapper"]["status"] == ("verified_wrapper_failed" if phase == 1 else
                                           "runtime_changed_or_unverifiable")
    assert result["wrapper"]["scientific_results_admitted"] is False
    assert (fixture.tmp / "panel/sessions/run_12/result.json").is_file()


@pytest.mark.parametrize("field,value", [("Comment", "wrong"), ("JobName", "other"),
                                         ("UserId", "other(999999)")])
def test_running_scheduler_identity_rechecked_before_attempt(fixture, monkeypatch, field, value):
    ref, calls, fields = execution_stubs(fixture, monkeypatch)
    fields[field] = value
    with pytest.raises(ValueError):
        module.execute(ref)
    assert not (fixture.tmp / "panel/sessions/run_12").exists()
    assert not Path(fixture.run["output_root"]).exists()
    assert not calls


def test_controller_does_not_reuse_output(fixture, monkeypatch):
    ref, calls, _ = execution_stubs(fixture, monkeypatch)
    Path(fixture.run["output_root"]).mkdir()
    with pytest.raises(ValueError):
        module.execute(ref)
    assert not calls


def release_fixture(fixture, monkeypatch, *, capacity=None, process_valid=True, comment=None):
    from benchmark_tools import probe_host_counters, observe_threadripper_process_identity
    from benchmark_tools import review_threadripper_process_policy, slurm_resource_snapshot

    process_ref = store(fixture.tmp / "process_policy.json", dict(fixture_policy=True))
    policy = module.read(fixture.request["policy"])
    policy["process_policy"] = process_ref
    fixture.request["policy"] = store(fixture.tmp / "release_policy.json", policy)
    ref = store(module.REQUEST, fixture.request)
    directory = fixture.tmp / "release"
    directory.mkdir()
    first = dict(snapshot=dict(first=True), observer_pid=88)
    (directory / "host_processes.jsonl").write_text(json.dumps(first) + "\n", encoding="ascii")
    allowed = list(range(32))
    store(directory / "ready.json", dict(allocated_placement=dict(stub=True),
        placement=dict(affinity=allowed), cgroup="fixture_cgroup"))
    calls = []
    monkeypatch.setattr(module, "amendment", lambda _: calls.append("amendment"))
    monkeypatch.setattr(module, "prepared_inputs", lambda *a: calls.append("prepared_inputs"))
    monkeypatch.setattr(module, "validate_placement", lambda *a: allowed)
    monkeypatch.setattr(module, "available_memory", lambda _: module.MEMORY if capacity is None else capacity)
    monkeypatch.setattr(probe_host_counters, "snapshot", lambda: dict(counters=True))
    monkeypatch.setattr(observe_threadripper_process_identity, "enriched_snapshot", lambda: dict(second=True))
    monkeypatch.setattr(slurm_resource_snapshot, "scoped_path", lambda *a: Path("/slice/job_25000/user/task"))

    def attribution(*args, **kwargs):
        calls.append(("process_attribution", args, kwargs))
        return dict(process_policy_matched=process_valid, observed_foreign_average_cores=150)

    monkeypatch.setattr(review_threadripper_process_policy, "review", attribution)

    class BudgetStub:
        def __init__(self, job, **kwargs):
            calls.append(("budget", job, kwargs))

        def __call__(self, destination):
            budget = dict(allocation=dict(fields=dict(Comment=ref["sha256"] if comment is None else comment,
                JobName=module.JOB_NAME)))
            store(destination / "release_budget.json", dict(stub=True, budget=budget))
            return budget

    monkeypatch.setattr(module, "ReleaseBudgetGuard", BudgetStub)
    return module.ReleaseGuard(ref, fixture.request, fixture.run, {}), directory, calls


def test_release_gate_accepts_recorded_competition_not_an_isolation_claim(fixture, monkeypatch):
    guard, directory, calls = release_fixture(fixture, monkeypatch)
    guard(directory)
    preflight = json.loads((directory / "environment_preflight.json").read_text())
    observation = json.loads((directory / "launch_environment_observation.json").read_text())
    assert preflight["controller_schema"] == module.REQUEST_SCHEMA
    assert preflight["background_competition_recorded"] is True
    assert preflight["uncontended_timing"] is False
    assert observation["background_cpu_used_for_eligibility"] is False
    assert observation["process_review"]["observed_foreign_average_cores"] == 150
    assert observation["capacity_guaranteed_through_run"] is False
    budget_call = next(call for call in calls if isinstance(call, tuple) and call[0] == "budget")
    assert budget_call == ("budget", 25000, dict(command=str(module.SCRIPT), cwd=str(module.ROOT),
                                               allocation_mode="shared"))
    attribution = next(call for call in calls if isinstance(call, tuple) and call[0] == "process_attribution")
    assert attribution[1][1:] == (dict(first=True), dict(second=True))
    assert attribution[2]["observer_pid"] == 88


@pytest.mark.parametrize("failure", ["memory", "attribution", "comment", "request_changed"])
def test_release_gate_rejects_unsafe_or_unbound_launch(fixture, monkeypatch, failure):
    guard, directory, calls = release_fixture(fixture, monkeypatch,
        capacity=module.MEMORY - 1 if failure == "memory" else None,
        process_valid=failure != "attribution", comment="wrong" if failure == "comment" else None)
    if failure == "request_changed":
        request = module.read(guard.request_ref)
        request["cell"] = "p0_c0_r0"
        store(module.REQUEST, request)
    with pytest.raises(ValueError):
        guard(directory)
    assert not (directory / "environment_preflight.json").exists()
    if failure in {"memory", "attribution"}:
        assert (directory / "launch_environment_observation.json").is_file()
        assert not any(isinstance(call, tuple) and call[0] == "budget" for call in calls)


def test_batch_fixed_envelope_and_unscheduled_guard(tmp_path):
    raw = module.SCRIPT.read_text()
    for value in ("--cpus-per-task=64", "--mem=128G", "--time=1-02:00:00", "--no-requeue",
                  "PYTHONHASHSEED=0", "PYTHONDONTWRITEBYTECODE=1", "scheduler-comment"):
        assert value in raw
    assert str(module.SCRIPT) not in raw  # imported exact batch identity, not historical batch
    result = subprocess.run(["bash", "-n", str(module.SCRIPT)], capture_output=True, text=True)
    assert result.returncode == 0
    env = dict(os.environ)
    env.pop("SLURM_JOB_ID", None)
    result = subprocess.run(["bash", str(module.SCRIPT)], env=env, cwd=tmp_path,
                            capture_output=True, text=True)
    assert result.returncode != 0
    assert "SLURM_JOB_ID" in result.stderr
