import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
from types import SimpleNamespace

import pytest

from benchmark_tools import run_threadripper_scaling as driver
from benchmark_tools.prepare_scaling_inputs import planned_runs
from benchmark_tools.prepare_threadripper_overhead import build
from tests.unit.test_verify_threadripper_controller import RAW


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


def private_setup(setup, monkeypatch):
    root, request, runs = setup
    tools = root / "benchmark_tools"
    results = tools / "results"
    script = tools / "run_threadripper_private_scaling.sh"
    script.write_text("# synthetic private submission\n")
    plan = put(results / "threadripper_private_commands_20260928.json", {"runs": runs, "deployment": "synthetic private interpreter"})
    controller_path = root / "synthetic_controller_python"
    controller_path.write_bytes(b"synthetic controller identity; never executed\n")
    monkeypatch.setattr(driver.sys, "executable", str(controller_path))
    controller = driver.record(controller_path)
    binding = put(root / "binding.json", dict(runtime_specs=[["fixture", "digest"]],
        controller_python=controller, baseline=driver.record(root / "baseline.json"),
        command_plan=plan, baseline_paths=[], retired_roots=[]))
    lookup = put(results / "threadripper_private_lookup_v2_20260928.json", dict(
        baseline=driver.record(root / "baseline.json"), binding=binding))
    monkeypatch.setattr(driver, "PRIVATE_PLAN_SHA", plan["sha256"])
    monkeypatch.setattr(driver, "PRIVATE_LOOKUP_SHA", lookup["sha256"])
    recipe = put(root / "recipe.json", dict(schema="threadripper_executor_recipe_v1", root=str(root),
        sources=[driver.record(tools / "run_threadripper_scaling.py"), driver.record(script)]))
    request.update(deployment="private_v2_20260928", plan_sha256=plan["sha256"], lookup_sha256=lookup["sha256"],
                   scheduler_command=str(script), recipe=recipe)
    policy = json.loads((root / "environment_policy.json").read_text())
    policy["plan_sha256"] = plan["sha256"]
    policy.update(schema="threadripper_environment_policy_v2", native_pressure_role="diagnostic_only")
    policy_ref = put(root / "environment_policy.json", policy)
    ready = json.loads((root / "ready.json").read_text())
    ready.update(plan_sha256=plan["sha256"], lookup_sha256=lookup["sha256"], recipe_sha256=recipe["sha256"],
                 environment_policy=policy_ref)
    request["readiness_review"] = put(root / "ready.json", ready)
    preflight = json.loads((root / "preflight.json").read_text())
    preflight.update(plan_sha256=plan["sha256"], recipe_sha256=recipe["sha256"], readiness_review_sha256=request["readiness_review"]["sha256"])
    put(root / "preflight.json", preflight)
    monkeypatch.setenv("ORTHOHMM_THREADRIPPER_DEPLOYMENT", "private_v2_20260928")
    return root, request, runs


def test_explicit_private_selection_and_gate_preserve_identity(setup, monkeypatch):
    root, request, runs = private_setup(setup, monkeypatch)
    run, lookup, history, sources = driver.select(request, root, 42)
    assert run == runs[0] and history["progress"]["index"] == 0
    assert lookup["binding"] == driver.record(root / "binding.json")
    assert not (root / "runs").exists() and len(sources) == 3
    gate, directory, calls = guard((root, request, runs))
    assert gate(directory) == {"budget": "checked"} and calls == [directory]


@pytest.mark.parametrize("fault", ["old_plan", "old_lookup", "old_script", "no_deployment", "wrong_ready"])
def test_private_request_never_mixes_shared_bindings(setup, monkeypatch, fault):
    root, request, _ = private_setup(setup, monkeypatch)
    if fault == "old_plan": request["plan_sha256"] = driver.PLAN_SHA
    elif fault == "old_lookup": request["lookup_sha256"] = driver.LOOKUP_SHA
    elif fault == "old_script": request["scheduler_command"] = str(root / "benchmark_tools/run_threadripper_scaling.sh")
    elif fault == "no_deployment": request.pop("deployment")
    else:
        ready = json.loads((root / "ready.json").read_text())
        ready["plan_sha256"] = driver.PLAN_SHA
        request["readiness_review"] = put(root / "ready.json", ready)
    with pytest.raises(ValueError):
        driver.select(request, root, 42)
    assert not (root / "runs").exists()


@pytest.mark.parametrize("name", [None, True, "private_latest", "", []])
def test_unknown_deployment_refused(name):
    with pytest.raises(ValueError, match="Unknown retained"):
        driver.deployment(dict(deployment=name))


def explicit_lookup_setup(setup, monkeypatch):
    root, request, runs = private_setup(setup, monkeypatch)
    retained_path = root / "benchmark_tools/results/threadripper_private_lookup_v2_20260928.json"
    retained = json.loads(retained_path.read_text())
    refreshed = dict(retained, status="native_lookup_repeated_identity_match",
                     scientific_execution_authorized=False, refresh="synthetic only")
    ref = put(root / "refreshed/lookup.json", refreshed)
    request.update(runtime_lookup=ref, lookup_sha256=ref["sha256"])
    ready_path = root / "ready.json"
    ready = json.loads(ready_path.read_text())
    ready["lookup_sha256"] = ref["sha256"]
    request["readiness_review"] = put(ready_path, ready)
    return root, request, runs


def overhead_setup(setup, monkeypatch):
    root, request, _ = explicit_lookup_setup(setup, monkeypatch)
    parent_path = root / "benchmark_tools/results/threadripper_private_commands_20260928.json"
    real_parent = Path(__file__).resolve().parents[2] / "benchmark_tools/results/threadripper_private_commands_20260928.json"
    parent = json.loads(real_parent.read_text())
    for row in parent["runs"]:
        row["cwd"] = str(root)
    parent_ref = put(parent_path, parent)
    monkeypatch.setattr(driver, "PRIVATE_PLAN_SHA", parent_ref["sha256"])
    binding_path = root / "binding.json"
    binding = json.loads(binding_path.read_text())
    binding["command_plan"] = parent_ref
    binding_ref = put(binding_path, binding)
    for lookup_path in (root / "benchmark_tools/results/threadripper_private_lookup_v2_20260928.json",
                        Path(request["runtime_lookup"]["path"])):
        lookup = json.loads(lookup_path.read_text())
        lookup["binding"] = binding_ref
        ref = put(lookup_path, lookup)
        if lookup_path == Path(request["runtime_lookup"]["path"]):
            request.update(runtime_lookup=ref, lookup_sha256=ref["sha256"])
        else:
            monkeypatch.setattr(driver, "PRIVATE_LOOKUP_SHA", ref["sha256"])
    plan = build(parent, root / "engineering", Path("/dev/shm") / ("engineering-" + root.name))
    plan.update(sources=[parent_ref, driver.record(root / "support.json")], helpers=[])
    plan_ref = put(root / "overhead_plan.json", plan)
    request.update(overhead_plan=plan_ref, plan_sha256=plan_ref["sha256"],
        environment_preflight_path=str(root / "engineering/sessions/run_00/environment_preflight.json"))
    policy_path = root / "environment_policy.json"
    policy = json.loads(policy_path.read_text())
    policy["plan_sha256"] = plan_ref["sha256"]
    policy_ref = put(policy_path, policy)
    ready_path = root / "ready.json"
    ready = json.loads(ready_path.read_text())
    ready.update(schema="threadripper_overhead_readiness_review_v1", plan_sha256=plan_ref["sha256"],
        lookup_sha256=request["lookup_sha256"], environment_policy=policy_ref,
        boundary_control_validated=True, incremental_overhead_measurement_pending=True)
    request["readiness_review"] = put(ready_path, ready)
    return root, request, [task["run"] for task in plan["runs"]]


def overhead_prior(root, request, index):
    plan = json.loads(Path(request["overhead_plan"]["path"]).read_text())
    history = []
    for task in plan["runs"][:index]:
        n, job = task["index"], 1000 + task["index"]
        directory = root / f"synthetic_history_{n:02d}"
        raw = RAW.replace("JobId=42", f"JobId={job}").replace("RUNNING", "COMPLETED")
        raw = raw.replace("Command=/recipe/run.sh", "Command=" + request["scheduler_command"])
        raw = raw.replace("WorkDir=/recipe", "WorkDir=" + str(root)).replace("1-00:00:00", driver.TIME_LIMIT)
        controller = put(directory / "controller.json", dict(
            command=["scontrol", "show", "job", str(job), "--oneliner"], returncode=0, stdout=raw))
        reviews = {k: put(directory / (k + ".json"), dict(schema="threadripper_overhead_review_v1",
            index=n, pair=task["pair"], arm=task["arm"], job_id=job,
            plan_sha256=request["plan_sha256"], category=k, decision="passed",
            evidence=[driver.record(root / "support.json")]))
            for k in ("runtime", "environment", "resources", "outputs_or_failure")}
        native = put(directory / "native.json", dict(job_id=job, native_outcome="exited_zero",
            status="native_success_outputs_verified"))
        history.append(put(directory / "session.json", dict(schema="threadripper_overhead_session_v1",
            index=n, pair=task["pair"], arm=task["arm"], job_id=job,
            plan_sha256=request["plan_sha256"], controller=controller, reviews=reviews,
            native_audit=native, native_outcome="exited_zero", phase="terminal")))
    request.update(index=index, history=history,
        environment_preflight_path=str(root / f"engineering/sessions/run_{index:02d}/environment_preflight.json"))


@pytest.mark.parametrize("index", [1, 26, 27, 53])
def test_overhead_later_identity_uses_its_arm_and_reviewed_prefix(setup, monkeypatch, index):
    root, request, runs = overhead_setup(setup, monkeypatch)
    overhead_prior(root, request, index)
    run, _, history, _ = driver.select(request, root, 42)
    assert run == runs[index] and history['progress']['index'] == index
    assert history['engineering_task']['pair'] == index // 2
    assert len(history['attempts']) == index
    assert not (root / 'engineering').exists()


@pytest.mark.parametrize('fault', ['session_schema', 'session_pair', 'review_schema', 'review_arm', 'native_failure', 'missing_review'])
def test_overhead_history_never_borrows_production_or_skips_failure(setup, monkeypatch, fault):
    root, request, _ = overhead_setup(setup, monkeypatch)
    overhead_prior(root, request, 1)
    session_path = Path(request['history'][0]['path'])
    session = json.loads(session_path.read_text())
    if fault == 'session_schema': session['schema'] = 'threadripper_panel_session_v1'
    elif fault == 'session_pair': session['pair'] = True
    elif fault == 'missing_review': session['reviews'] = None
    elif fault == 'native_failure':
        session['native_outcome'] = 'exited_nonzero'
        session['native_audit'] = put(Path(session['native_audit']['path']), dict(job_id=1000,
            native_outcome='exited_nonzero', status='native_failure_requires_review'))
    else:
        path = Path(session['reviews']['runtime']['path'])
        review = json.loads(path.read_text())
        if fault == 'review_schema': review['schema'] = 'threadripper_panel_review_v1'
        else: review['arm'] = 'periodic' if review['arm'] == 'boundary' else 'boundary'
        session['reviews']['runtime'] = put(path, review)
    request['history'][0] = put(session_path, session)
    with pytest.raises(ValueError): driver.select(request, root, 42)
    assert not (root / 'engineering').exists()


def test_overhead_selection_preserves_native_work_and_separate_identity(setup, monkeypatch):
    root, request, runs = overhead_setup(setup, monkeypatch)
    run, _, history, evidence = driver.select(request, root, 42)
    assert run == runs[0]
    assert history["engineering_task"]["pair"] == 0
    assert history["engineering_task"]["arm"] in {"boundary", "periodic"}
    assert request["overhead_plan"] in evidence
    assert not (root / "engineering").exists()
    gate, directory, calls = guard((root, request, runs))
    preflight_path = root / "preflight.json"
    preflight = json.loads(preflight_path.read_text())
    preflight.update(plan_sha256=request["plan_sha256"], readiness_review_sha256=request["readiness_review"]["sha256"])
    put(preflight_path, preflight)
    assert gate(directory) == {"budget": "checked"} and calls == [directory]


def test_retained_overhead_plan_matches_private_parent_without_native_work():
    results = Path(__file__).resolve().parents[2] / "benchmark_tools/results"
    parent_path = results / "threadripper_private_commands_20260928.json"
    plan_path = results / "threadripper_native_overhead_plan_20260930.json"
    parent = driver.read_frozen(parent_path, driver.PRIVATE_PLAN_SHA)
    plan = driver.read_frozen(plan_path, "84affc2274bf3594f679661b95b49705fa36e2577c9b12594b4643e388dd3182")
    driver.overhead_design(plan, parent)
    assert plan['sources'][0]['sha256'] == driver.PRIVATE_PLAN_SHA
    assert plan['runs'][0]['arm'] == 'boundary' and plan['runs'][1]['arm'] == 'periodic'
    assert driver.overhead_position(plan['runs'], [])['index'] == 0


@pytest.mark.parametrize("fault", ["shared", "no_lookup", "size", "hash", "symlink", "missing_task",
    "arm", "command", "parent", "production_ready", "boundary_unvalidated", "already_measured", "bool_index", "index_54"])
def test_overhead_selection_rejects_mixed_or_drifted_work(setup, monkeypatch, fault):
    root, request, _ = overhead_setup(setup, monkeypatch)
    if fault == "shared": request["deployment"] = "shared_v3_20260928"
    elif fault == "no_lookup": request.pop("runtime_lookup")
    elif fault == "size": request["overhead_plan"]["bytes"] += 1
    elif fault == "hash": request["overhead_plan"]["sha256"] = "0" * 64
    elif fault == "symlink":
        link = root / "plan-link.json"
        link.symlink_to(request["overhead_plan"]["path"])
        request["overhead_plan"]["path"] = str(link)
    elif fault in {"bool_index", "index_54"}:
        request["index"] = False if fault == "bool_index" else 54
    elif fault in {"production_ready", "boundary_unvalidated", "already_measured"}:
        path = root / "ready.json"
        ready = json.loads(path.read_text())
        if fault == "production_ready": ready["schema"] = "threadripper_readiness_review_v1"
        elif fault == "boundary_unvalidated": ready["boundary_control_validated"] = False
        else: ready["incremental_overhead_measurement_pending"] = False
        request["readiness_review"] = put(path, ready)
    else:
        path = Path(request["overhead_plan"]["path"])
        plan = json.loads(path.read_text())
        if fault == "missing_task": plan["runs"].pop()
        elif fault == "arm": plan["runs"][0]["arm"] = "periodic" if plan["runs"][0]["arm"] == "boundary" else "boundary"
        elif fault == "command": plan["runs"][0]["run"]["native_argv"].append("--different")
        else: plan["sources"][0] = driver.record(root / "support.json")
        request["overhead_plan"] = put(path, plan)
        request["plan_sha256"] = request["overhead_plan"]["sha256"]
    with pytest.raises((ValueError, KeyError)):
        driver.select(request, root, 42)
    assert not (root / "engineering").exists()


def test_explicit_runtime_lookup_preserves_plan_controller_and_old_pin(setup, monkeypatch):
    root, request, runs = explicit_lookup_setup(setup, monkeypatch)
    retained_path = root / "benchmark_tools/results/threadripper_private_lookup_v2_20260928.json"
    retained = driver.record(retained_path)
    selected = driver.deployment(request)
    assert selected["lookup"] == request["runtime_lookup"]["path"]
    assert selected["lookup_sha"] == request["runtime_lookup"]["sha256"]
    assert selected["plan_sha"] == driver.PRIVATE_PLAN_SHA
    run, lookup, _, evidence = driver.select(request, root, 42)
    assert run == runs[0] and lookup["baseline"] == driver.record(root / "baseline.json")
    assert request["runtime_lookup"] in evidence and retained in evidence
    assert driver.record(retained_path) == retained
    assert not (root / "runs").exists()


@pytest.mark.parametrize("fault", ["null", "extra", "relative", "symlink", "size", "digest",
                                  "uppercase", "wrong_bytes", "wrong_hash", "baseline", "controller", "readiness",
                                  "lookup_status", "lookup_permission", "binding_baseline", "command_plan",
                                  "baseline_paths", "retired_roots", "private_manifest"])
def test_explicit_runtime_lookup_rejects_invalid_or_mixed_bindings(setup, monkeypatch, fault):
    root, request, _ = explicit_lookup_setup(setup, monkeypatch)
    ref = request["runtime_lookup"]
    if fault == "null": request["runtime_lookup"] = None
    elif fault == "extra": ref["discovered"] = True
    elif fault == "relative": ref["path"] = "refreshed/lookup.json"
    elif fault == "symlink":
        link = root / "lookup-link.json"
        link.symlink_to(ref["path"])
        ref["path"] = str(link)
    elif fault == "size": ref["bytes"] = True
    elif fault == "digest": ref["sha256"] = "invalid"
    elif fault == "uppercase": ref["sha256"] = "A" * 64
    elif fault == "wrong_bytes": ref["bytes"] += 1
    elif fault == "wrong_hash": ref["sha256"] = "0" * 64
    elif fault == "readiness":
        ready = json.loads((root / "ready.json").read_text())
        ready["lookup_sha256"] = driver.PRIVATE_LOOKUP_SHA
        request["readiness_review"] = put(root / "ready.json", ready)
    else:
        data = json.loads(Path(ref["path"]).read_text())
        if fault == "baseline": data["baseline"] = put(root / "other-baseline.json", {})
        elif fault == "lookup_status": data["status"] = "native_lookup_inspected"
        elif fault == "lookup_permission": data["scientific_execution_authorized"] = True
        else:
            binding = json.loads((root / "binding.json").read_text())
            if fault == "controller":
                binding["controller_python"] = dict(binding["controller_python"], path=str(root / "other-python"))
            elif fault == "binding_baseline": binding["baseline"] = put(root / "other-baseline.json", {})
            elif fault == "command_plan": binding["command_plan"] = dict(binding["command_plan"], sha256="0" * 64)
            elif fault == "private_manifest": binding["runtime_specs"].append(["other-private", "digest"])
            else: binding[fault] = ["changed"]
            data["binding"] = put(root / "other-binding.json", binding)
        request["runtime_lookup"] = put(Path(ref["path"]), data)
        request["lookup_sha256"] = request["runtime_lookup"]["sha256"]
    with pytest.raises((ValueError, KeyError)):
        driver.select(request, root, 42)
    assert not (root / "runs").exists()


@pytest.mark.parametrize("deployment", [None, "shared_v3_20260928", "unknown"])
def test_explicit_lookup_cannot_change_shared_or_unknown_deployment(deployment, tmp_path):
    request = dict(runtime_lookup=put(tmp_path / "lookup.json", {}))
    if deployment is not None:
        request["deployment"] = deployment
    with pytest.raises(ValueError, match="retained private"):
        driver.deployment(request)


def test_private_submission_bootstrap_does_not_fall_back_to_shared(setup, monkeypatch):
    root, request, _ = private_setup(setup, monkeypatch)
    ref = put(root / "request.json", request)
    monkeypatch.delenv("ORTHOHMM_THREADRIPPER_DEPLOYMENT")
    with pytest.raises(ValueError, match="bootstrap"):
        driver.execute(Path(ref["path"]), ref["sha256"])
    assert not (root / "runs").exists()


@pytest.mark.parametrize("relocated", [False, True])
def test_private_pin_constants_match_retained_artifacts_and_bootstrap(tmp_path, relocated):
    original = Path(__file__).resolve().parents[2]
    profile = driver.deployment(dict(deployment="private_v2_20260928"))
    root = tmp_path / "relocated checkout" if relocated else original
    if relocated:
        for name in ("results/" + profile["plan"], "results/" + profile["lookup"],
                     "results/threadripper_private_controller_20260928.json", profile["script"]):
            target = root / "benchmark_tools" / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(original / "benchmark_tools" / name, target)
    results = root / "benchmark_tools/results"
    assert driver.record(results / profile["plan"])["sha256"] == profile["plan_sha"]
    assert driver.record(results / profile["lookup"])["sha256"] == profile["lookup_sha"]
    document = json.loads((results / "threadripper_private_controller_20260928.json").read_text())
    controller = document["interpreter"]
    script = root / "benchmark_tools" / profile["script"]
    text = script.read_text()
    assert controller["sha256"] in text
    retained_root = Path(document["source"]["path"]).parents[1]
    interpreter_suffix = Path(controller["path"]).relative_to(retained_root).as_posix()
    assert f"ROOT={retained_root}" in text
    assert f'PYTHON="$ROOT/{interpreter_suffix}"' in text
    assert "/home/bizon/anaconda3" not in text
    assert "ORTHOHMM_THREADRIPPER_DEPLOYMENT=private_v2_20260928" in text
    for token in ("--exclusive", "--cpus-per-task=64", "--mem=128G", "--time=1-02:00:00", "--no-requeue"):
        assert token in text
    syntax = subprocess.run(["bash", "-n", str(script)], capture_output=True, text=True, timeout=10)
    assert syntax.returncode == 0, syntax.stderr


@pytest.mark.parametrize("fault", ["missing", "wrong_path", "wrong_hash", "legacy_policy"])
def test_private_controller_failure_precedes_session_or_native_work(setup, monkeypatch, fault):
    root, request, _ = private_setup(setup, monkeypatch)
    binding = json.loads((root / "binding.json").read_text())
    if fault == "legacy_policy":
        policy_path = root / "environment_policy.json"
        policy = json.loads(policy_path.read_bytes())
        policy.update(schema="threadripper_environment_policy_v1")
        policy.pop("native_pressure_role")
        policy_ref = put(policy_path, policy)
        ready_path = root / "ready.json"
        ready = json.loads(ready_path.read_bytes())
        ready["environment_policy"] = policy_ref
        request["readiness_review"] = put(ready_path, ready)
    elif fault == "missing": binding.pop("controller_python")
    elif fault == "wrong_path": binding["controller_python"]["path"] = str(root / "other_controller")
    else: binding["controller_python"]["sha256"] = "0" * 64
    binding_ref = put(root / "binding.json", binding)
    lookup_path = root / "benchmark_tools/results/threadripper_private_lookup_v2_20260928.json"
    lookup = json.loads(lookup_path.read_text())
    lookup["binding"] = binding_ref
    lookup_ref = put(lookup_path, lookup)
    request["lookup_sha256"] = lookup_ref["sha256"]
    monkeypatch.setattr(driver, "PRIVATE_LOOKUP_SHA", lookup_ref["sha256"])
    ready = json.loads((root / "ready.json").read_text())
    ready["lookup_sha256"] = lookup_ref["sha256"]
    request["readiness_review"] = put(root / "ready.json", ready)
    ref = put(root / "request.json", request)
    monkeypatch.setattr(driver, "__file__", str(root / "benchmark_tools/run_threadripper_scaling.py"))
    monkeypatch.setattr(driver.sys, "dont_write_bytecode", True)
    monkeypatch.setattr(driver.sys, "pycache_prefix", "/dev/shm/orthohmm_scaling_driver_42")
    monkeypatch.setattr(driver.os, "uname", lambda: SimpleNamespace(nodename="bizon"))
    for key, value in dict(SLURM_JOB_ID="42", SLURM_CPUS_PER_TASK="64", SLURM_MEM_PER_NODE="131072", PYTHONHASHSEED="0").items():
        monkeypatch.setenv(key, value)
    monkeypatch.chdir(root)
    with pytest.raises(ValueError, match="controller interpreter|Frozen input/source|Private timing requires"):
        driver.execute(Path(ref["path"]), ref["sha256"])
    assert not (root / "runs").exists()


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
@pytest.mark.parametrize("private,native_v2", [(False, False), (False, True), (True, True), ("explicit", True), ("overhead", True), ("overhead_second", True)])
def test_execute_one_preserves_attempt_and_never_submits(setup, monkeypatch, raises, worker_fails, stream_fails, private, native_v2):
    if private in {"overhead", "overhead_second"}:
        setup = overhead_setup(setup, monkeypatch)
        if private == "overhead_second": overhead_prior(setup[0], setup[1], 1)
    elif private == "explicit":
        setup = explicit_lookup_setup(setup, monkeypatch)
    elif private:
        setup = private_setup(setup, monkeypatch)
    root, request, runs = setup
    active = runs[request["index"]]
    if native_v2:
        policy_path = root / "environment_policy.json"
        policy = json.loads(policy_path.read_bytes())
        policy.update(schema="threadripper_environment_policy_v2", native_pressure_role="diagnostic_only")
        policy_ref = put(policy_path, policy)
        ready_path = root / "ready.json"
        ready = json.loads(ready_path.read_bytes())
        ready["environment_policy"] = policy_ref
        request["readiness_review"] = put(ready_path, ready)
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
    runtime_bindings = []
    monkeypatch.setattr(driver, "RuntimeChecker", lambda *args: runtime_bindings.append(args) or "checker")
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
        assert run == active and job == 42 and kwargs["runtime_checker"] == "checker"
        selected = driver.deployment(request)
        assert runtime_bindings[0][:2] == (root / "benchmark_tools/results" / selected["lookup"], selected["lookup_sha"])
        assert isinstance(kwargs["release_guard"], driver.EnvironmentalReleaseGuard)
        if private in {"overhead", "overhead_second"}:
            task = json.loads(Path(request["overhead_plan"]["path"]).read_text())["runs"][request["index"]]
            assert collector is (driver.measure_boundary if task["arm"] == "boundary" else driver.measure)
        else:
            assert collector is driver.measure
        assert worker_events == ["started"]
        if raises:
            raise RuntimeError("synthetic infrastructure failure")
        kwargs["release_guard"].budget(root)
        put(Path(request["environment_preflight_path"]), dict(synthetic_test_only=True))
        return {"status": "command_exited_zero"}
    monkeypatch.setattr(driver, "measure_run", measured)
    def audit_stream(directory, policy_ref, preflight_ref, **kwargs):
        assert directory == active["measurement_directory"]
        expected = dict(job_id=42, index=request["index"])
        if private in {"overhead", "overhead_second"}:
            task = json.loads(Path(request["overhead_plan"]["path"]).read_text())["runs"][request["index"]]
            expected["collector_arm"] = task["arm"]
        assert kwargs == expected
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
    saved = json.loads(Path(request["environment_preflight_path"]).with_name("result.json").read_text())
    if private in {"overhead", "overhead_second"}:
        assert saved["purpose"] == "native_overhead" and not saved["production_identity"]
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
