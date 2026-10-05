"""Terminal joins use synthetic evidence; no job submission or inference."""

from copy import deepcopy
import json
from pathlib import Path
import tempfile

import pytest

from benchmark_tools import review_native_factorial_attempt as module
from benchmark_tools.snapshot_runtime_trees import inventory


def store(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))
    return module.record(path)


@pytest.fixture
def runtime(tmp_path):
    source = tmp_path / "native.so"
    source.write_bytes(b"synthetic scientific runtime")
    manifest = store(tmp_path / "runtime.json", inventory([source]))
    specs = [[manifest["path"], manifest["sha256"]]]
    expected = [dict(path=manifest["path"], sha256=manifest["sha256"],
        status="runtime_tree_identity_matches", records=1, scientific_execution_authorized=False)]
    report = dict(python="fixture version", executable="/private/python", cwd="/frozen/core",
        requested=["orthohmm.orthohmm"], paths=[], modules={}, mapped_files=[],
        meta_path=[], path_hooks=[], editable={}, files=[], dont_write_bytecode=True,
        coverage=dict(missing=[], changed=[], all_covered=True), scientific_origin=dict(root="/frozen/core"))
    prior = store(tmp_path / "prior.json", report)
    baseline = store(tmp_path / "baseline.json", {})
    binding = store(tmp_path / "binding.json", dict(runtime_specs=specs))
    lookup = store(tmp_path / "lookup.json", dict(baseline=baseline, binding=binding,
        source=module.record(module.ROOT / "benchmark_tools/inspect_native_python_lookup.py"),
        interpreters={name: dict(reports=[prior]) for name in ("orthohmm", "orthofinder")}))
    original = tmp_path / "source.fa"
    original.write_text(">gene\nACDE\n")
    root = tmp_path / "run"
    (root / "input").mkdir(parents=True)
    copy = root / "input/source.fa"
    copy.write_bytes(original.read_bytes())
    run = dict(inputs=[module.record(original)], output_root=str(root), native_order=["source.fa"])
    session = tmp_path / "session"
    verification = dict(before_check_wall_s=.1, after_check_wall_s=.2)
    for number, side in enumerate(("before", "after"), start=1):
        comparisons = {}
        for name in ("orthohmm", "orthofinder"):
            observed = store(session / f"lookup_checks/check_{number:02d}/{name}.json", report)
            comparisons[name] = dict(module.compare_lookup(report, report), report=observed)
        checked = dict(status="runtime_and_lookup_checked", scientific_execution_authorized=False,
                       runtime=expected, lookup=comparisons)
        store(session / f"lookup_checks/checked_{number:02d}.json", checked)
        prepared = None if side == "before" else dict(datasets=[dict(native_order=["source.fa"],
            inputs_in_native_order=[module.record(copy)])])
        verification[side] = dict(runtime=checked, original_inputs=run["inputs"], prepared_inputs=prepared)
    return dict(plan=dict(runtime_lookup=lookup, baseline=baseline), run=run, session=session,
                verification=verification, source=source)


def run_runtime(data):
    evidence = module.Evidence()
    result = module.runtime_review(data["plan"], data["run"], data["session"], data["verification"], evidence)
    evidence.finish()
    return result


def test_runtime_brackets_lookup_and_fresh_inventory(runtime):
    result = run_runtime(runtime)
    assert result["status"] == "runtime_brackets_and_lookup_replayed"
    assert result["fresh_terminal_inventory"][0]["records"] == 1
    assert result["continuous_runtime_integrity_established"] is False


@pytest.mark.parametrize("change", ["record_count", "input", "precopy", "postcopy", "order", "lookup", "current_tree"])
def test_runtime_evidence_tampering(runtime, change):
    data = runtime
    verification = data["verification"]
    if change == "record_count":
        verification["before"]["runtime"]["runtime"][0]["records"] = 2
    elif change == "input":
        verification["after"]["original_inputs"] = []
    elif change == "precopy":
        verification["before"]["prepared_inputs"] = {}
    elif change == "postcopy":
        verification["after"]["prepared_inputs"]["datasets"][0]["inputs_in_native_order"] = []
    elif change == "order":
        verification["after"]["prepared_inputs"]["datasets"][0]["native_order"] = ["wrong.fa"]
    elif change == "lookup":
        path = data["session"] / "lookup_checks/check_02/orthohmm.json"
        value = json.loads(path.read_text())
        value["modules"] = {"orthohmm": "/untrusted/orthohmm.py"}
        path.write_text(json.dumps(value))
    else:
        data["source"].write_bytes(b"changed runtime")
    with pytest.raises(ValueError):
        run_runtime(data)


def outcome(kind="exited_zero"):
    successful = kind == "exited_zero"
    terminal = dict(verified=dict(fields=dict(JobState="COMPLETED" if successful else "FAILED",
                                             ExitCode="0:0" if successful else "1:0")))
    result = dict(status="measurement_returned_pending_independent_review" if successful else "factorial_attempt_failed_retained")
    verification = dict(status="command_exited_zero" if successful else "command_failed")
    return terminal, result, verification, dict(native_outcome=kind)


@pytest.mark.parametrize("kind,expected", [("exited_zero", "native_success"),
    ("exited_nonzero", "native_failure_retained"), ("timed_out", "native_timeout_retained")])
def test_native_failure_classes_preserved(kind, expected):
    assert module.classify(*outcome(kind)) == expected


@pytest.mark.parametrize("change", ["scheduler", "exit_code", "executor", "postflight", "timeout_signal"])
def test_infrastructure_failure_not_mislabeled_native(change):
    terminal, result, verification, replayed = outcome()
    if change == "scheduler": terminal["verified"]["fields"]["JobState"] = "FAILED"
    elif change == "exit_code": terminal["verified"]["fields"]["ExitCode"] = "1:0"
    elif change == "executor": result["status"] = "factorial_executor_failed"
    elif change == "postflight": verification["status"] = "runtime_changed_or_unverifiable"
    else: replayed["native_outcome"] = "timed_out"
    with pytest.raises(ValueError):
        module.classify(terminal, result, verification, replayed)


def test_accounting_terminal_format():
    args = outcome()
    args[0]["verified"] = dict(State="COMPLETED", ExitCode="0:0")
    assert module.classify(*args) == "native_success"


@pytest.fixture
def environment(tmp_path):
    directory = tmp_path / "measurement"
    directory.mkdir()
    job, observer, boot = 42, 20, "fixture-boot"
    scope = "/slurm/job_42"
    process_policy = dict(schema="threadripper_process_policy_v3", boot_id=boot,
        execution_scope=module.SCOPE, review_reference="synthetic fixture, not host evidence")
    process_ref = store(tmp_path / "process_policy.json", process_policy)
    plan_ref = dict(path="/synthetic/plan.json", bytes=1, sha256="plan-digest")
    policy = dict(schema="threadripper_environment_policy_v3", execution_scope=module.SCOPE,
        native_pressure_role="diagnostic_only", preflight_pressure_role="diagnostic_only",
        foreign_cpu_role="diagnostic_only", plan_sha256=plan_ref["sha256"], process_policy=process_ref,
        maximum_foreign_average_cores=0., maximum_sample_period_s=3.,
        maximum_pressure_percent=dict(cpu=0., memory=0., io=0.), maximum_pressure_sample_period_s=3.,
        minimum_available_memory_bytes=module.MEMORY, configuration_files=[], evidence=[])
    policy_ref = store(tmp_path / "policy.json", policy)
    rows, points = [], []
    for index, t in enumerate((10., 12., 14.)):
        processes = []
        for pid, group in ((10, "/other_analysis"), (observer, scope + "/step_batch")):
            processes.append(dict(pid=pid, created=float(pid), cgroup=group, name="fixture_process",
                user_s=100. * index if pid == 10 else 0., system_s=0., observed_monotonic_s=t + .1,
                kernel_identity=dict(pid=pid, tgid=pid, kthread=0, started_monotonic_s=t + .15,
                                     finished_monotonic_s=t + .2)))
        rows.append(dict(index=index, observer_pid=observer, snapshot=dict(
            schema="threadripper_typed_process_snapshot_v1", boot_id=boot,
            started_monotonic_s=t, finished_monotonic_s=t + .3, errors=[], processes=processes)))
        samples = [dict(started_monotonic_ns=int((t + offset) * 1e9),
            finished_monotonic_ns=int((t + offset + .1) * 1e9), errors=[],
            raw=dict(boot_id=boot + "\n", online_cpus="0-31\n", cgroup_membership="0::" + scope + "/step_0\n"),
            optional={"host_" + resource + "_pressure":
                "some avg10=0.00 avg60=0.00 avg300=0.00 total=%d\nfull avg10=0.00 avg60=0.00 avg300=0.00 total=0\n" % (index * 300000)
                for resource in ("cpu", "memory", "io")}) for offset in (0., .2)]
        points.append(dict(host=samples))
        store(directory / f"point_{index:06d}.json", points[-1])
    stream = directory / "host_processes.jsonl"
    stream.write_text("".join(json.dumps(row) + "\n" for row in rows))
    done = dict(started_ns=11_000_000_000, finished_ns=13_000_000_000)
    store(directory / "ready.json", dict(cgroup="0::" + scope + "/step_0/user/task_0\n"))
    request_ref = dict(path="/synthetic/request", bytes=1, sha256="request-digest")
    raw = (f"JobId=42 JobState=RUNNING Partition=gpu NodeList=bizon NumNodes=1 NumCPUs=64 "
        f"NumTasks=1 CPUs/Task=64 OverSubscribe=OK MinMemoryNode=128G Requeue=0 Restarts=0 "
        f"Command={module.SCRIPT} WorkDir={module.ROOT} TimeLimit=1-02:00:00 ExitCode=0:0 "
        f"RunTime=00:00:10 Comment={request_ref['sha256']}")
    budget = module.remaining_budget(raw, job, command=str(module.SCRIPT), cwd=str(module.ROOT),
                                    query_elapsed_s=.1, allocation_mode="shared")
    budget_ref = store(directory / "release_budget.json", dict(status="release_budget_check_passed",
        returncode=0, stdout=raw, query_elapsed_s=.1, budget=budget))
    launch_verdict = module.review_process_pair(process_policy, rows[0]["snapshot"], rows[1]["snapshot"],
        boot_id=boot, job_scope=scope, observer_pid=observer)
    raw_memory = "MemAvailable: 200000000 kB\n"
    launch_ref = store(directory / "launch_environment_observation.json", dict(
        process_snapshots=[rows[0]["snapshot"], rows[1]["snapshot"]], process_review=launch_verdict,
        raw_meminfo=raw_memory, available_memory_bytes=module.available_memory(raw_memory),
        background_cpu_used_for_eligibility=False))
    preflight = dict(schema="threadripper_environment_preflight_v1", decision="passed",
        job_id=job, index=0, environment_policy=policy_ref, plan_sha256=plan_ref["sha256"],
        execution_scope=module.SCOPE, background_competition_recorded=True, boot_id=boot,
        evidence=[launch_ref, budget_ref])
    store(directory / "environment_preflight.json", preflight)
    processes = module.process_evaluate((json.dumps(row) for row in rows), process_policy,
        boot_id=boot, job_scope=scope, observer_pid=observer, launch=11., end=13.,
        maximum_foreign_average_cores=0., maximum_sample_period_s=3.)
    pressure = module.pressure_evaluate(iter(points), boot_id=boot, job_scope=scope,
        launch_ns=done["started_ns"], end_ns=done["finished_ns"], limits=policy["maximum_pressure_percent"],
        maximum_period_s=3., pressure_role="diagnostic_only")
    retained = dict(processes, job_id=job, index=0, uncontended_timing=False,
        background_cpu_used_for_eligibility=False, pressure_thresholds_used_for_eligibility=False,
        pressure_review=pressure, sampled_environment_policy_satisfied=True, evidence=[],
        source=module.record(module.ROOT / "benchmark_tools/review_threadripper_process_stream.py"))
    retained_ref = store(directory / "process_stream_review.json", retained)
    return dict(request_ref=request_ref, request=dict(job_id=job, policy=policy_ref), plan_ref=plan_ref,
        run=dict(index=0), result=dict(environment_review=retained_ref, sampled_environment_evidence_valid=True),
        directory=directory, done=done)


def run_environment(data):
    evidence = module.Evidence()
    result = module.environment_review(**data, evidence=evidence)
    evidence.finish()
    return result


def test_shared_contention_magnitudes_not_rejection(environment):
    result = run_environment(environment)
    assert result["sampled_environment_evidence_valid"] is True
    assert result["processes"]["maximum_observed_foreign_average_cores"] == 50.
    assert result["pressure"]["diagnostic_thresholds_satisfied"] is False
    assert result["background_cpu_used_for_eligibility"] is False
    assert result["pressure_thresholds_used_for_eligibility"] is False
    assert result["uncontended_timing"] is False


@pytest.mark.parametrize("change", ["launch_capacity", "launch_counter", "budget_comment", "retained_cores", "stream", "point"])
def test_environment_raw_evidence_tampering(environment, change):
    data = environment
    directory = data["directory"]
    if change == "stream":
        path = directory / "host_processes.jsonl"
        rows = [json.loads(line) for line in path.read_text().splitlines()]
        rows[1]["observer_pid"] = 999
        path.write_text("".join(json.dumps(row) + "\n" for row in rows))
    else:
        name = {"launch_capacity": "launch_environment_observation.json", "launch_counter": "launch_environment_observation.json",
                "budget_comment": "release_budget.json", "retained_cores": "process_stream_review.json", "point": "point_000001.json"}[change]
        path = directory / name
        value = json.loads(path.read_text())
        if change == "launch_capacity":
            value.update(raw_meminfo="MemAvailable: 1 kB\n", available_memory_bytes=1024)
        elif change == "launch_counter": value["process_review"]["process_policy_matched"] = False
        elif change == "budget_comment": value["stdout"] = value["stdout"].replace("request-digest", "wrong")
        elif change == "retained_cores": value["maximum_observed_foreign_average_cores"] = 0.
        else: value["host"][0]["raw"]["boot_id"] = "other"
        path.write_text(json.dumps(value))
        # Re-seal the specific referenced objects: semantic replay, not merely
        # a checksum mismatch, must still reject a forged consistent summary.
        if change.startswith("launch") or change == "budget_comment":
            preflight = json.loads((directory / "environment_preflight.json").read_text())
            preflight["evidence"] = [module.record(r["path"]) for r in preflight["evidence"]]
            store(directory / "environment_preflight.json", preflight)
        elif change == "retained_cores":
            data["result"]["environment_review"] = module.record(path)
    with pytest.raises(ValueError):
        run_environment(data)


def test_resource_scope_includes_wrappers_without_subtraction(monkeypatch):
    derived = dict(primary=dict(wall_seconds=10., cpu_seconds=100., peak_memory_bytes=1000), primary_scopes=module.SCOPES)
    monkeypatch.setattr(module, "endpoints", lambda *args: deepcopy(derived))
    replayed = dict(affinity_observation_statuses=["observed_within_affinity", "incomplete"], narrow_flagged_intervals=[1])
    result = module.resource_review(replayed, {}, 42)
    assert result["primary"] == derived["primary"]
    assert result["narrow_flagged_intervals"] == [1]
    assert result["whole_job_cpu_seconds"] is None
    assert result["uncontended_timing"] is False
    replayed["affinity_observation_statuses"].append("violation")
    with pytest.raises(ValueError, match="affinity violation"):
        module.resource_review(replayed, {}, 42)


def test_native_command_is_bound_to_wrapper_not_cli():
    assert module.native_command(dict(path="/plan", sha256="digest"), dict(index=2),
        dict(tool_entrypoints=dict(orthohmm_python=dict(absolute_path="/python")))) == [
        "/python", "-B", str(module.ROOT / "benchmark_tools/run_native_factorial_cost.py"),
        "--native", "--plan", "/plan", "--plan-sha256", "digest", "--index", "2"]


def test_save_produces_an_immutable_reference(tmp_path):
    path = tmp_path / "result.json"
    ref = module.save(path, dict(status="synthetic"))
    assert ref == module.record(path)
    with pytest.raises(FileExistsError):
        module.save(path, {})


def test_live_poll_failure_does_not_create_review_or_restart(tmp_path, monkeypatch):
    request_ref = store(tmp_path / "request.json", dict(plan={"path": "/plan"}, job_id=42))
    monkeypatch.setattr(module, "read", lambda ref: dict(plan={"path": "/plan"}, job_id=42))
    monkeypatch.setattr(module, "validate_plan", lambda plan: [])
    monkeypatch.setattr(module, "validate_request", lambda *args: None)
    def running(job):
        raise ValueError("Require terminal controller evidence")
    monkeypatch.setattr(module, "verify_terminal", running)
    destination = tmp_path / "review"
    with pytest.raises(ValueError, match="terminal"):
        module.review(request_ref, destination)
    assert not destination.exists()


@pytest.fixture
def joined(monkeypatch):
    with tempfile.TemporaryDirectory(prefix="factorial_terminal_test_", dir=module.ROOT / "benchmarks/work") as temporary:
        base = Path(temporary)
        panel = base / "panel"
        root = panel / "run_00"
        session = panel / "sessions/run_00"
        run = dict(index=0, dataset="orthobench", cell="p0_c0_r0", repeat=0, output_root=str(root))
        baseline = store(base / "baseline.json", dict(tool_entrypoints=dict(
            orthohmm_python=dict(absolute_path="/synthetic/private/python"))))
        plan_ref = store(base / "plan.json", dict(panel_root=str(panel), baseline=baseline,
            helper_sources=[], evidence=[]))
        request_ref = store(base / "request.json", dict(schema="native_factorial_cost_request_v1",
            execution_authorized=True, job_id=42, index=0, plan=plan_ref, history=[],
            scheduler_command=str(module.SCRIPT), allocation_cwd=str(module.ROOT)))
        done = dict(exit_code=0, timed_out=False, started_ns=10_000_000_000, finished_ns=12_000_000_000)
        store(root / "measurement/done.json", done)
        measured = dict(native=done)
        verification = dict(status="command_exited_zero", scientific_results_admitted=False,
            source_sha256=module.record(module.ROOT / "benchmark_tools/run_verified_slurm_measurement.py")["sha256"],
            measurement=measured)
        result = dict(status="measurement_returned_pending_independent_review", job_id=42, index=0,
            cell=run["cell"], plan=plan_ref, request=request_ref, execution_scope=module.SCOPE,
            automatic_retry=False, next_identity_authorized=False, uncontended_timing=False, wrapper=verification)
        terminal = dict(source="live_controller", verified=dict(fields=dict(JobState="COMPLETED",
            ExitCode="0:0", Comment=request_ref["sha256"])))
        replayed = dict(native_outcome="exited_zero", native_exit_code=0, measured=measured,
            evidence=[module.record(root / "measurement/done.json")],
            affinity_observation_statuses=["observed_within_affinity"], narrow_flagged_intervals=[])
        derived = dict(primary=dict(wall_seconds=2., cpu_seconds=40., peak_memory_bytes=1000000), primary_scopes=module.SCOPES)
        environment = dict(sampled_environment_evidence_valid=True,
                           processes=dict(maximum_observed_foreign_average_cores=99.))
        outputs = dict(job_id=42, index=0, plan=plan_ref, native_outputs_validated=True,
                       evidence=[], checked_files=[], accuracy_evaluated=False)
        monkeypatch.setattr(module, "validate_plan", lambda plan: [run])
        monkeypatch.setattr(module, "verify_terminal", lambda job: deepcopy(terminal))
        monkeypatch.setattr(module, "runtime_review", lambda *args: dict(status="synthetic_kernel_result"))
        commands = []
        def replay(directory, job, command):
            commands.append(command)
            return deepcopy(replayed)
        monkeypatch.setattr(module, "replay", replay)
        monkeypatch.setattr(module, "endpoints", lambda *args: deepcopy(derived))
        monkeypatch.setattr(module, "environment_review", lambda *args: deepcopy(environment))
        monkeypatch.setattr(module, "validate_outputs", lambda ref: deepcopy(outputs))
        yield dict(base=base, root=root, session=session, request=request_ref, verification=verification,
            result=result, terminal=terminal, replayed=replayed, environment=environment,
            outputs=outputs, commands=commands)


def run_joined(data):
    store(data["root"] / "verification.json", data["verification"])
    store(data["session"] / "result.json", data["result"])
    return module.review(data["request"], data["base"] / "review")


@pytest.mark.parametrize("kind", ["exited_zero", "exited_nonzero", "timed_out"])
def test_joined_review_history_contract_and_retained_failures(joined, kind):
    data = joined
    if kind != "exited_zero":
        data["terminal"]["verified"]["fields"].update(JobState="FAILED", ExitCode="1:0")
        data["result"]["status"] = "factorial_attempt_failed_retained"
        data["verification"]["status"] = "command_failed"
        data["replayed"].update(native_outcome=kind, native_exit_code=124 if kind == "timed_out" else 1)
    ref = run_joined(data)
    report = module.read(ref)
    assert report["terminal_reviewed"] is True
    assert report["next_identity_authorized"] is True
    assert report["index"] == 0 and report["job_id"] == 42
    assert report["native_outputs_validated"] is (kind == "exited_zero")
    assert report["accuracy_evaluated"] is False
    assert report["automatic_retry"] is False
    assert report["uncontended_timing"] is False
    assert report["resources"] == dict(wall_seconds=2., cpu_seconds=40., peak_memory_bytes=1000000)
    assert report["whole_run_maximum_foreign_average_cores"] == 99.
    assert data["commands"][0][1] == "-B"
    for category in ("runtime", "resources", "environment", "outputs_or_failure"):
        module.read(report["reviews"][category])


def test_invalid_environment_blocks_next_not_by_contention_magnitude(joined):
    joined["environment"]["sampled_environment_evidence_valid"] = False
    report = module.read(run_joined(joined))
    assert report["terminal_reviewed"] is True
    assert report["next_identity_authorized"] is False
    assert report["shared_host_resources_reviewed"] is False
    assert report["primary_resources_replayed"] is True


@pytest.mark.parametrize("change", ["session", "collector", "outputs", "runtime"])
def test_join_failure_retained_without_next_or_retry(joined, monkeypatch, change):
    if change == "session": joined["result"]["job_id"] = 999
    elif change == "collector": joined["replayed"]["measured"] = {"wrong": True}
    elif change == "outputs": joined["outputs"]["plan"] = {}
    else:
        def broken(*args):
            raise ValueError("synthetic runtime drift")
        monkeypatch.setattr(module, "runtime_review", broken)
    with pytest.raises(ValueError):
        run_joined(joined)
    failure = json.loads((joined["base"] / "review/failure.json").read_text())
    assert failure["terminal_reviewed"] is False
    assert failure["next_identity_authorized"] is False
    assert failure["automatic_retry"] is False
    assert not (joined["base"] / "review/review.json").exists()
