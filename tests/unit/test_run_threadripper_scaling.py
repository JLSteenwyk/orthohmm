import json
from pathlib import Path
import subprocess
import sys
import time
from types import SimpleNamespace

import pytest

from benchmark_tools import run_threadripper_scaling as driver
from benchmark_tools.prepare_scaling_inputs import planned_runs


def put(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))
    return driver.record(path)


@pytest.fixture
def setup(tmp_path, monkeypatch, synthetic_linux_boot_id):
    tools = tmp_path / "benchmark_tools"
    tools.mkdir()
    (tools / "run_threadripper_scaling.py").write_text("# synthetic executor\n")
    (tools / "run_threadripper_scaling.sh").write_text("# synthetic submission\n")
    baseline = put(tmp_path / "baseline.json", {})
    binding = put(tmp_path / "binding.json", {"runtime_specs": [["fixture", "digest"]]})
    results = tools / "results"
    runs = planned_runs()
    for row in runs:
        row.update(cwd=str(tmp_path), measurement_directory=str(tmp_path / "runs" / f"run_{row['index']:02d}" / "measurement"))
    plan = put(results / "threadripper_scaling_commands_20260928.json", {"runs": runs})
    lookup = put(results / "threadripper_python_lookup_v3_20260928.json", dict(baseline=baseline, binding=binding))
    monkeypatch.setattr(driver, "PLAN_SHA", plan["sha256"])
    monkeypatch.setattr(driver, "LOOKUP_SHA", lookup["sha256"])
    recipe = put(tmp_path / "recipe.json", dict(schema="threadripper_executor_recipe_v1", root=str(tmp_path),
        sources=[driver.record(p) for p in sorted(tools.iterdir()) if p.is_file()]))
    support = put(tmp_path / "support.json", {"synthetic_test_only": True})
    policy = put(tmp_path / "environment_policy.json", dict(schema="threadripper_environment_policy_v1",
        decision="reviewed", host="bizon", plan_sha256=plan["sha256"],
        maximum_foreign_average_cores=.1, maximum_sample_period_s=35.,
        maximum_pressure_sample_period_s=3., maximum_pressure_percent=dict(cpu=10., io=10., memory=10.),
        review_reference="synthetic test only", evidence=[support]))
    ready = put(tmp_path / "ready.json", dict(schema="threadripper_readiness_review_v1", decision="passed",
        plan_sha256=plan["sha256"], lookup_sha256=lookup["sha256"], recipe_sha256=recipe["sha256"],
        full_scale_observer_validated=True, environment_policy_frozen=True,
        review_reference="synthetic test, not real readiness", evidence=[support], environment_policy=policy))
    preflight = put(tmp_path / "preflight.json", dict(schema="threadripper_environment_preflight_v1",
        decision="passed", job_id=42, index=0, plan_sha256=plan["sha256"], recipe_sha256=recipe["sha256"],
        readiness_review_sha256=ready["sha256"], whole_run_observer_ready=True,
        unrelated_scientific_work_present=False, observation_started_unix_ns=10**12,
        observation_finished_unix_ns=10**12+10**9,
        boot_id=synthetic_linux_boot_id,
        review_reference="synthetic test only", evidence=[support]))
    request = dict(schema="threadripper_execution_request_v1", job_id=42, execution_authorized=True,
        plan_sha256=plan["sha256"], lookup_sha256=lookup["sha256"], index=0,
        allocation_cwd=str(tmp_path), scheduler_command=str(tools / "run_threadripper_scaling.sh"),
        recipe=recipe, history=[], readiness_review=ready,
        environment_preflight_path=str(tmp_path / "runs/sessions/run_00/environment_preflight.json"))
    return tmp_path, request, runs


def test_select_preserves_exact_run_and_has_no_side_effects(setup):
    root, request, runs = setup
    run, lookup, history, sources = driver.select(request, root, 42)
    assert run == runs[0]
    assert history["progress"]["index"] == 0
    assert not (root / "runs").exists()
    assert len(sources) == 3


@pytest.mark.parametrize("key,value", [("job_id", 43), ("index", True), ("index", -1),
    ("index", 27), ("index", 1), ("execution_authorized", 1), ("execution_authorized", False),
    ("plan_sha256", "wrong"), ("lookup_sha256", "wrong"),
    ("scheduler_command", "/other.sh"), ("allocation_cwd", "/elsewhere")])
def test_wrong_request_refused(setup, key, value):
    root, request, _ = setup
    request[key] = value
    with pytest.raises(ValueError):
        driver.select(request, root, 42)


def test_uninventoried_helper_refused(setup):
    root, request, _ = setup
    (root / "benchmark_tools/new.py").write_text("pass\n")
    with pytest.raises(ValueError, match="inventory"):
        driver.select(request, root, 42)


def test_changed_helper_refused(setup):
    root, request, _ = setup
    (root / "benchmark_tools/run_threadripper_scaling.py").write_text("changed\n")
    with pytest.raises(ValueError):
        driver.select(request, root, 42)


@pytest.mark.parametrize("key,value", [("full_scale_observer_validated", False),
    ("environment_policy_frozen", False), ("review_reference", ""), ("evidence", []),
    ("recipe_sha256", "other")])
def test_unready_panel_refused(setup, key, value):
    root, request, _ = setup
    path = Path(request["readiness_review"]["path"])
    data = json.loads(path.read_text())
    data[key] = value
    request["readiness_review"] = put(path, data)
    with pytest.raises(ValueError):
        driver.select(request, root, 42)


def guard(setup, clock=None, budget=None):
    root, request, _ = setup
    ref = put(root / "request.json", request)
    directory = root / "measurement"
    directory.mkdir()
    calls = []
    budget = budget or (lambda p: calls.append(p) or {"budget": "checked"})
    times = iter([10**12, 10**12+2*10**9, 10**12+2*10**9])
    def waiter(path, seconds):
        assert seconds == 20
        put(path, json.loads((root / "preflight.json").read_text()))
    return driver.EnvironmentalReleaseGuard(request, ref, budget,
        clock=clock or (lambda: next(times)), waiter=waiter), directory, calls


def test_release_calls_budget_only_after_bound_review(setup):
    gate, directory, calls = guard(setup)
    assert gate(directory) == {"budget": "checked"}
    assert calls == [directory]
    receipt = json.loads((directory / "environment_release.json").read_text())
    assert receipt["status"] == "environment_review_bound"
    assert receipt["observational_validity_independently_established"] is False


@pytest.mark.parametrize("key,value", [("job_id", 41), ("index", 1),
    ("boot_id", "different"), ("unrelated_scientific_work_present", True),
    ("whole_run_observer_ready", False), ("observation_started_unix_ns", True),
    ("observation_finished_unix_ns", 10**15), ("observation_started_unix_ns", 1),
    ("readiness_review_sha256", "different")])
def test_release_refuses_bad_or_stale_preflight(setup, key, value):
    root, request, _ = setup
    path = root / "preflight.json"
    data = json.loads(path.read_text())
    data[key] = value
    put(path, data)
    gate, directory, calls = guard(setup)
    with pytest.raises(ValueError):
        gate(directory)
    assert not calls
    assert json.loads((directory / "environment_release.json").read_text())["status"] == "environment_release_failed"


def test_review_expiring_during_budget_is_refused(setup):
    times = iter([10**12, 10**12+2*10**9, 10**12+121*10**9])
    gate, directory, calls = guard(setup, clock=lambda: next(times))
    with pytest.raises(ValueError, match="expired"):
        gate(directory)
    assert len(calls) == 1


def test_changed_support_refused_at_release(setup):
    root, _, _ = setup
    gate, directory, calls = guard(setup)
    (root / "support.json").write_text("{}")
    with pytest.raises(ValueError):
        gate(directory)
    assert not calls


def test_missing_review_times_out_without_budget_check(setup):
    gate, directory, calls = guard(setup)
    def absent(path, seconds):
        assert seconds == 20
        raise TimeoutError("no external review")
    gate.waiter = absent
    with pytest.raises(TimeoutError):
        gate(directory)
    assert not calls
    assert (directory / "environment_review_requested.json").is_file()


def test_preexisting_review_is_not_reused(setup):
    root, request, _ = setup
    gate, directory, calls = guard(setup)
    put(Path(request["environment_preflight_path"]), {})
    with pytest.raises(FileExistsError):
        gate(directory)
    assert not calls


def test_source_change_at_gate_is_refused(setup):
    root, _, _ = setup
    gate, directory, calls = guard(setup)
    source = root / "benchmark_tools/run_threadripper_scaling.py"
    gate.evidence = [driver.record(source)]
    source.write_text("changed")
    with pytest.raises(ValueError):
        gate(directory)
    assert not calls


def test_real_file_handshake_binds_new_synthetic_review(setup):
    root, request, _ = setup
    gate, directory, calls = guard(setup, clock=time.time_ns)
    gate.waiter = driver.wait_file
    code = """
import json, pathlib, sys, time
root, directory = map(pathlib.Path, sys.argv[1:])
signal = directory / 'environment_review_requested.json'
deadline = time.monotonic() + 10
while not signal.exists():
    if time.monotonic() > deadline:
        raise TimeoutError('gate not observed')
    time.sleep(.01)
request = json.loads(signal.read_text())
data = json.loads((root / 'preflight.json').read_text())
data['observation_started_unix_ns'] = time.time_ns()
data['observation_finished_unix_ns'] = time.time_ns()
target = pathlib.Path(request['review_path'])
target.parent.mkdir(parents=True, exist_ok=True)
temporary = target.with_suffix('.tmp')
temporary.write_text(json.dumps(data))
temporary.rename(target)
"""
    process = subprocess.Popen([sys.executable, "-I", "-B", "-c", code, str(root), str(directory)],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        assert gate(directory) == {"budget": "checked"}
        out, err = process.communicate(timeout=15)
        assert process.returncode == 0, (out, err)
    finally:
        if process.poll() is None:
            process.kill()
            process.communicate()
    assert calls == [directory]
    assert json.loads((directory / "environment_release.json").read_text())["review"] == driver.record(
        request["environment_preflight_path"])


@pytest.mark.parametrize("raises", [False, True])
@pytest.mark.parametrize("worker_fails", [False, True])
@pytest.mark.parametrize("stream_fails", [False, True])
def test_execute_one_preserves_attempt_and_never_submits(setup, monkeypatch, raises, worker_fails, stream_fails):
    root, request, runs = setup
    ref = put(root / "request.json", request)
    monkeypatch.setattr(driver, "__file__", str(root / "benchmark_tools/run_threadripper_scaling.py"))
    monkeypatch.setattr(driver.sys, "dont_write_bytecode", True)
    monkeypatch.setattr(driver.sys, "pycache_prefix", "/dev/shm/orthohmm_scaling_driver_42")
    monkeypatch.setattr(driver.os, "uname", lambda: SimpleNamespace(nodename="bizon"))
    monkeypatch.setenv("SLURM_JOB_ID", "42")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "64")
    monkeypatch.setenv("SLURM_MEM_PER_NODE", "131072")
    monkeypatch.setenv("PYTHONHASHSEED", "0")
    monkeypatch.setenv("PYTHONDONTWRITEBYTECODE", "1")
    monkeypatch.setenv("PYTHONPYCACHEPREFIX", "/dev/shm/orthohmm_scaling_driver_42")
    monkeypatch.chdir(root)
    monkeypatch.setattr(driver, "execution_environment", lambda baseline: ({}, None))
    monkeypatch.setattr(driver, "RuntimeChecker", lambda *args: "checker")
    calls = []
    worker_events = []
    class Worker:
        def __init__(self, session, *args): self.session = session
        def __enter__(self):
            worker_events.append("started")
            return self
        def finish(self):
            worker_events.append("joined")
            if worker_fails: raise RuntimeError("synthetic worker failure")
        def wait_response(self, path, seconds): raise AssertionError("Synthetic measurement does not request preflight")
        def __exit__(self, kind, error, traceback):
            worker_events.append("cleaned")
            put(self.session / "environment_worker_lifecycle.json", dict(synthetic_test_only=True))
    monkeypatch.setattr(driver, "EnvironmentWorker", Worker)
    monkeypatch.setattr(driver, "ReleaseBudgetGuard", lambda *args, **kwargs:
                        lambda directory: worker_events.append("budget"))
    def measured(run, baseline, specs, collector, job, **kwargs):
        calls.append(run)
        assert run == runs[0] and job == 42 and kwargs["runtime_checker"] == "checker"
        assert isinstance(kwargs["release_guard"], driver.EnvironmentalReleaseGuard)
        assert worker_events == ["started"]
        if raises:
            raise RuntimeError("synthetic infrastructure failure")
        kwargs["release_guard"].budget(root)
        put(Path(request["environment_preflight_path"]), dict(synthetic_test_only=True))
        return {"status": "command_exited_zero"}
    monkeypatch.setattr(driver, "measure_run", measured)
    def audit_stream(directory, policy_ref, preflight_ref, **kwargs):
        assert directory == runs[0]["measurement_directory"]
        assert kwargs == dict(job_id=42, index=0)
        assert worker_events[-1] == "cleaned"
        worker_events.append("stream_review")
        return put(root / "stream_review.json", {}), dict(sampled_environment_policy_satisfied=not stream_fails)
    monkeypatch.setattr(driver, "audit_process_stream", audit_stream)
    if raises or worker_fails:
        with pytest.raises(RuntimeError, match="synthetic"):
            driver.execute(Path(ref["path"]), ref["sha256"])
    elif stream_fails:
        with pytest.raises(ValueError, match="sampled environment policy"):
            driver.execute(Path(ref["path"]), ref["sha256"])
    else:
        result = driver.execute(Path(ref["path"]), ref["sha256"])
        assert result["status"] == "measurement_returned_pending_independent_review"
    saved = json.loads((root / "runs/sessions/run_00/result.json").read_text())
    assert not saved["scientific_timings_admitted"] and not saved["next_submission_authorized"]
    with pytest.raises(FileExistsError):
        driver.execute(Path(ref["path"]), ref["sha256"])
    assert len(calls) == 1
    expected = ["started"] + ([] if raises else ["joined"] + ([] if worker_fails else ["budget"])) + ["cleaned"]
    if not raises and not worker_fails:
        expected.append("stream_review")
        assert saved["process_stream_review"]["path"].endswith("stream_review.json")
    assert worker_events == expected
    assert saved["environment_worker"]["path"].endswith("environment_worker_lifecycle.json")
