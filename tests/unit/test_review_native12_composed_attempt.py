"""Small real runtime/policy kernels and stub full joins; no production reviews."""

from copy import deepcopy
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from types import SimpleNamespace

import pytest

from benchmark_tools import review_native12_composed_attempt as module
from tests.unit import test_review_native_factorial_attempt as historical
from tests.unit import test_native12_composed_execution as compatible
from tests.unit.test_run_native12_composed_cost import held_raw
from tests.unit import test_native_factorial_outputs as semantic_tests
from tests.unit.test_run_allocated_native_factorial_cost import ready


def store(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding="ascii")
    return module.record(path)


@pytest.fixture
def runtime(tmp_path, monkeypatch):
    roots, specs, basis, _ = compatible.real_trees.__wrapped__(tmp_path, monkeypatch)
    fresh = module.current_trees(specs, basis)
    report = dict(python="fixture version", executable="/private/python", cwd="/frozen/core",
        requested=[], paths=[], modules={}, mapped_files=[], meta_path=[], path_hooks=[], editable={},
        files=[], dont_write_bytecode=True, coverage=dict(missing=[], changed=[], all_covered=True),
        scientific_origin=dict(root="/frozen/core"))
    prior = store(tmp_path / "prior.json", report)
    baseline = store(tmp_path / "baseline.json", {})
    binding = store(tmp_path / "binding.json", dict(runtime_specs=specs))
    lookup = store(tmp_path / "lookup.json", dict(baseline=baseline, binding=binding,
        source=module.record(module.ROOT / "benchmark_tools/inspect_native_python_lookup.py"),
        interpreters={name: dict(reports=[prior]) for name in ("orthohmm", "orthofinder")}))
    original = tmp_path / "source.fa"
    original.write_text(">gene\nACDE\n", encoding="ascii")
    root = tmp_path / "run"
    (root / "input").mkdir(parents=True)
    copied = root / "input/source.fa"
    copied.write_bytes(original.read_bytes())
    run = dict(inputs=[module.record(original)], output_root=str(root), native_order=["source.fa"])
    session = tmp_path / "session"
    verification = dict(before_check_wall_s=.1, after_check_wall_s=.2)
    for number, side in enumerate(("before", "after"), start=1):
        comparisons = {}
        for name in ("orthohmm", "orthofinder"):
            observed = store(session / f"lookup_checks/check_{number:02d}/{name}.json", report)
            comparisons[name] = dict(module.compare_lookup(report, report), report=observed)
        checked = dict(status="runtime_and_lookup_checked", scientific_execution_authorized=False,
                       runtime=fresh, lookup=comparisons)
        store(session / f"lookup_checks/checked_{number:02d}.json", checked)
        store(session / f"lookup_checks/tree_check_{number:02d}.json", dict(
            schema="native12_current_runtime_check_v1", source=module.record(
                module.ROOT / "benchmark_tools/native12_composed_execution.py"), runtime_basis=basis,
            inventories=fresh, original_os_inventory_equality=False, current_identity_revalidated=True,
            continuous_runtime_integrity_established=False, check_wall_s=.1))
        prepared = None if side == "before" else dict(datasets=[dict(native_order=["source.fa"],
            inputs_in_native_order=[module.record(copied)])])
        verification[side] = dict(runtime=deepcopy(checked), original_inputs=run["inputs"], prepared_inputs=prepared)
    return dict(plan=dict(runtime_lookup=lookup, baseline=baseline), run=run, session=session,
        verification=verification, request=dict(runtime_basis=basis), roots=roots)


def runtime_review(data):
    evidence = module.Evidence()
    args = {key: data[key] for key in ("plan", "run", "session", "verification", "request")}
    result = module.runtime_review(**args, evidence=evidence)
    evidence.finish()
    return result


def test_runtime_review_replays_real_trees_and_lookup_with_original_inequality(runtime):
    result = runtime_review(runtime)
    assert result["schema"] == "native12_composed_runtime_review_v1"
    assert result["current_original_os_inventory_equality"] is False
    assert result["prospective_current_inventory_equality"] is True
    assert result["continuous_runtime_integrity_established"] is False
    assert result["fresh_terminal_check_wall_s"] > 0


@pytest.mark.parametrize("mutation", ["precopy", "input", "postcopy", "order", "lookup",
    "tree_receipt", "changed_tree", "source", "count", "duration", "basis"])
def test_runtime_review_rejects_semantically_resealed_tampering(runtime, mutation):
    verification = runtime["verification"]
    if mutation == "precopy":
        verification["before"]["prepared_inputs"] = {}
    elif mutation == "input":
        verification["after"]["original_inputs"] = []
    elif mutation == "postcopy":
        verification["after"]["prepared_inputs"]["datasets"][0]["inputs_in_native_order"] = []
    elif mutation == "order":
        verification["after"]["prepared_inputs"]["datasets"][0]["native_order"] = ["wrong.fa"]
    elif mutation == "lookup":
        path = runtime["session"] / "lookup_checks/check_02/orthohmm.json"
        report = json.loads(path.read_text())
        report["modules"] = dict(untrusted="/outside/runtime")
        store(path, report)
    elif mutation == "changed_tree":
        (runtime["roots"][0] / "lftp").write_text("changed", encoding="ascii")
    elif mutation in {"source", "tree_receipt", "basis"}:
        path = runtime["session"] / "lookup_checks/tree_check_02.json"
        receipt = json.loads(path.read_text())
        receipt[{"source": "source", "tree_receipt": "original_os_inventory_equality", "basis": "runtime_basis"}[mutation]] = True
        store(path, receipt)
    elif mutation == "count":
        verification["after"]["runtime"]["runtime"][0]["comparison"]["observed_records"] += 1
    else:
        verification["before_check_wall_s"] = 0
    with pytest.raises(ValueError):
        runtime_review(runtime)


@pytest.fixture
def environment(tmp_path):
    data = historical.environment.__wrapped__(tmp_path)
    request = data["request"]
    request["plan"] = data.pop("plan_ref")
    request["amendment"] = dict(path="/fixture/amendment.json", bytes=1, sha256="fixture")
    data["run"]["index"] = 12
    directory = data["directory"]
    budget_path = directory / "release_budget.json"
    budget = json.loads(budget_path.read_text())
    budget["stdout"] = budget["stdout"].replace(str(historical.module.SCRIPT), str(module.executor.SCRIPT))
    budget["stdout"] += " JobName=" + module.executor.JOB_NAME
    budget["budget"] = module.remaining_budget(budget["stdout"], request["job_id"],
        command=str(module.executor.SCRIPT), cwd=str(module.ROOT), query_elapsed_s=.1, allocation_mode="shared")
    store(budget_path, budget)
    path = directory / "environment_preflight.json"
    preflight = json.loads(path.read_text())
    preflight.update(index=12, amendment=request["amendment"], controller_schema=module.executor.REQUEST_SCHEMA,
        evidence=[module.record(ref["path"]) for ref in preflight["evidence"]])
    store(path, preflight)
    path = directory / "process_stream_review.json"
    retained = json.loads(path.read_text())
    retained["index"] = 12
    data["result"]["environment_review"] = store(path, retained)
    return data


def environment_review(data):
    evidence = module.Evidence()
    result = module.environment_review(**data, evidence=evidence)
    evidence.finish()
    return result


def test_environment_replays_actual_policy_kernels_and_retains_competition(environment):
    result = environment_review(environment)
    assert result["sampled_environment_evidence_valid"] is True
    assert result["processes"]["maximum_observed_foreign_average_cores"] == 50.
    assert result["pressure"]["diagnostic_thresholds_satisfied"] is False
    assert result["background_cpu_used_for_eligibility"] is False
    assert result["uncontended_timing"] is False


@pytest.mark.parametrize("mutation", ["batch", "comment", "job_name", "controller", "amendment",
    "capacity", "summary", "preflight_index", "retained_index", "pressure", "snapshot"])
def test_environment_rejects_tampered_contract_and_accounting(environment, mutation):
    directory = environment["directory"]
    preflight_path = directory / "environment_preflight.json"
    preflight = json.loads(preflight_path.read_text())
    if mutation in {"batch", "comment", "job_name"}:
        path = directory / "release_budget.json"
        budget = json.loads(path.read_text())
        if mutation == "batch":
            budget["stdout"] = budget["stdout"].replace(str(module.executor.SCRIPT), str(historical.module.SCRIPT))
        elif mutation == "comment":
            budget["stdout"] = budget["stdout"].replace("request-digest", "wrong")
        else:
            budget["stdout"] = budget["stdout"].replace(module.executor.JOB_NAME, "other")
        store(path, budget)
    elif mutation in {"controller", "amendment", "preflight_index"}:
        preflight[{"controller": "controller_schema", "amendment": "amendment", "preflight_index": "index"}[mutation]] = "wrong"
    elif mutation in {"capacity", "snapshot"}:
        path = directory / "launch_environment_observation.json"
        launch = json.loads(path.read_text())
        if mutation == "capacity":
            launch.update(raw_meminfo="MemAvailable: 1 kB\n", available_memory_bytes=1024)
        else:
            launch["process_snapshots"][0] = {}
        store(path, launch)
    else:
        path = directory / "process_stream_review.json"
        retained = json.loads(path.read_text())
        if mutation == "summary":
            retained["maximum_observed_foreign_average_cores"] = 0
        elif mutation == "retained_index":
            retained["index"] = 11
        else:
            retained["pressure_review"]["diagnostic_thresholds_satisfied"] = True
        environment["result"]["environment_review"] = store(path, retained)
    preflight["evidence"] = [module.record(ref["path"]) for ref in preflight["evidence"]]
    store(preflight_path, preflight)
    with pytest.raises(ValueError):
        environment_review(environment)


def test_fresh_terminal_accounting_uses_original_envelope_parser(monkeypatch):
    raw = "25000|COMPLETED|0:0|64|128G|bizon|gpu|1560|2026-10-08T12:00:00|2026-10-08T12:01:00|2026-10-08T13:00:00|orthohmm_allocated_factorial\n"
    calls = []

    def query(command, **kwargs):
        calls.append(command)
        if command[0] == "scontrol":
            return SimpleNamespace(returncode=1, stdout="", stderr="Invalid job id specified")
        return SimpleNamespace(returncode=0, stdout=raw, stderr="")

    monkeypatch.setattr(module.subprocess, "run", query)
    result = module.verify_terminal(25000)
    assert result["source"] == "fresh_accounting_after_controller_expiry"
    assert result["verified"]["AllocCPUS"] == "64"
    assert [call[0] for call in calls] == ["scontrol", "sacct"]


def test_controller_query_failure_is_not_assumed_terminal(monkeypatch):
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k:
        SimpleNamespace(returncode=1, stdout="", stderr="Transient timeout"))
    with pytest.raises(ValueError):
        module.verify_terminal(25000)
    with pytest.raises(ValueError):
        module.verify_terminal(24034)


def test_live_terminal_controller_binds_the_new_batch_path(monkeypatch):
    raw = (f"JobId=25000 JobState=COMPLETED Partition=gpu NodeList=bizon NumNodes=1 NumCPUs=64 "
        f"NumTasks=1 CPUs/Task=64 OverSubscribe=OK MinMemoryNode=128G Requeue=0 Restarts=0 "
        f"Command={module.executor.SCRIPT} WorkDir={module.ROOT} TimeLimit=1-02:00:00 ExitCode=0:0 "
        f"JobName={module.executor.JOB_NAME}")
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k:
        SimpleNamespace(returncode=0, stdout=raw, stderr=""))
    result = module.verify_terminal(25000)
    assert result["verified"]["fields"]["Command"] == str(module.executor.SCRIPT)
    wrong = raw.replace(str(module.executor.SCRIPT), str(historical.module.SCRIPT))
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k:
        SimpleNamespace(returncode=0, stdout=wrong, stderr=""))
    with pytest.raises(ValueError):
        module.verify_terminal(25000)


@pytest.fixture
def scientific(tmp_path):
    # CPU/command metadata are resealed ONLY in this synthetic relocated copy.
    context = semantic_tests.copy_context(tmp_path, "p1_c1_r1")
    original = module.validate_semantics(context)
    root = Path(context["output_root"])
    shutil.copytree(context["input_directory"], root / "input")
    run = dict(context, index=12, dataset="qfo_corrected", repeat=0)
    baseline = dict(core_root=context["cwd"], tool_entrypoints=dict(
        orthohmm_python=dict(absolute_path=str(Path(sys.executable).absolute())),
        mafft=dict(absolute_path=context["aligner"]), FastTree=dict(absolute_path=context["tree_builder"])))
    baseline_ref = store(tmp_path / "baseline.json", baseline)
    plan_ref = store(tmp_path / "plan.json", dict(fixture=True))
    amendment_ref = store(tmp_path / "amendment.json", dict(fixture=True))
    request = dict(job_id=42, index=12, plan=plan_ref, amendment=amendment_ref, new_sources=[])
    request_ref = store(tmp_path / "request.json", request)
    ready_ref = store(root / "measurement/ready.json", ready())
    current = dict(deepcopy(ready()["placement"]), pid=2)
    allowed = ready()["placement"]["affinity"]
    native = dict(schema="allocated_native_factorial_execution_v1",
        status="native_factorial_completed_pending_output_review", plan=plan_ref, amendment=amendment_ref,
        index=12, cell="p1_c1_r1", placement=current, parent_pid=1, allocated_ready=ready_ref,
        native_cpu_ids=allowed, native_order=context["native_order"], automatic_retry=False,
        source=module.record(module.ROOT / "benchmark_tools/run_allocated_native_factorial_cost.py"),
        factors=original["factors"])
    store(root / "native_execution.json", native)
    store(root / "preparation.json", dict(status="fresh_factorial_inputs_prepared", genes=26,
        gene_ownership_sha256=original["gene_ownership_sha256"], per_species_counts=original["per_species_counts"]))
    metrics_path = root / "metrics.json"
    metrics = json.loads(metrics_path.read_text())
    metrics["metadata"].update(cpu_budget=32, search_total_threads=32, search_threads_per_worker=4,
        search_workers=8, fasta_directory=str(root / "input"))
    metrics["command"] = module.native_command(amendment_ref, run, baseline, metrics=True)
    store(metrics_path, metrics)
    provenance_path = root / "native/orthohmm_phylogeny/provenance_manifest.json"
    provenance = json.loads(provenance_path.read_text())
    provenance["cpu_budget"] = 32
    store(provenance_path, provenance)
    terminal = dict(verified=dict(fields=dict(JobState="COMPLETED", ExitCode="0:0")))
    return dict(request_ref=request_ref, request=request, execution=dict(new_sources=[]),
        plan=dict(baseline=baseline_ref, helper_sources=[]), run=run, baseline=baseline, terminal=terminal)


def test_scientific_output_uses_real_semantic_and_placement_kernels(scientific):
    result = module.scientific_output(**scientific)
    assert result["schema"] == "native12_composed_output_review_v1"
    assert result["native_outputs_validated"] is True
    assert result["input_genes"] == 26 and result["checkpoint"]["hits"] == 116
    assert result["phylogeny"]["native_pair_rows"] == 43
    assert result["semantic_validator_source"] == module.record(
        module.ROOT / "benchmark_tools/validate_native_factorial_outputs.py")
    assert result["accuracy_evaluated"] is False


@pytest.mark.parametrize("mutation", ["native_source", "native_index", "cpu", "factor",
    "ownership", "parameter", "terminal", "checkpoint"])
def test_scientific_output_rejects_provenance_and_real_semantic_tampering(scientific, mutation):
    root = Path(scientific["run"]["output_root"])
    if mutation in {"native_source", "native_index", "cpu", "factor"}:
        path = root / "native_execution.json"
        native = json.loads(path.read_text())
        if mutation == "native_source":
            native["source"] = {}
        elif mutation == "native_index":
            native["index"] = 11
        elif mutation == "cpu":
            native["native_cpu_ids"] = list(range(32))
        else:
            native["factors"] = {}
        store(path, native)
    elif mutation == "ownership":
        path = root / "preparation.json"
        preparation = json.loads(path.read_text())
        preparation["gene_ownership_sha256"] = "wrong"
        store(path, preparation)
    elif mutation == "parameter":
        path = root / "metrics.json"
        metrics = json.loads(path.read_text())
        metrics["metadata"]["cpm_resolution"] = .2
        store(path, metrics)
    elif mutation == "checkpoint":
        path = root / "native/orthohmm_working_res/high_sensitivity_checkpoint/gene_names.txt"
        path.write_text("tampered", encoding="ascii")
    else:
        scientific["terminal"]["verified"]["fields"]["JobState"] = "FAILED"
    with pytest.raises(ValueError):
        module.scientific_output(**scientific)


@pytest.fixture
def joined(tmp_path, monkeypatch):
    root, session, destination = tmp_path / "panel/run_12", tmp_path / "panel/sessions/run_12", tmp_path / "review"
    monkeypatch.setattr(module, "DESTINATION", destination)
    baseline = store(tmp_path / "baseline.json", dict(core_root=str(tmp_path / "core"),
        tool_entrypoints=dict(orthohmm_python=dict(absolute_path="/fixture/native/python"))))
    run = dict(index=12, cell="p1_c1_r1", dataset="qfo_corrected", repeat=0, output_root=str(root))
    plan = dict(panel_root=str(tmp_path / "panel"), baseline=baseline, runs=[{}] * 12 + [run],
        helper_sources=[], evidence=[])
    plan_ref = store(tmp_path / "plan.json", plan)
    amendment_ref = store(tmp_path / "amendment.json", {})
    original_ref = store(tmp_path / "original_request.json", {})
    composed_ref = store(tmp_path / "composed_review.json", {})
    history = dict(schema=module.executor.HISTORY_SCHEMA, original_request=original_ref, composed_review=composed_ref,
        historical_prefix=[], next_unrun_index=12, original_review_translated=False, next_identity_authorized=False, evidence=[])
    request = dict(schema=module.executor.REQUEST_SCHEMA, job_id=25000, index=12, cell=run["cell"], plan=plan_ref,
        amendment=amendment_ref, original_request=original_ref, composed_review=composed_ref, history=[composed_ref],
        history_basis=history, new_sources=[], held_scheduler=dict(stdout=held_raw()))
    request_ref = store(tmp_path / "request.json", request)
    done = dict(exit_code=0, timed_out=False, started_ns=10_000_000_000, finished_ns=12_000_000_000)
    store(root / "measurement/done.json", done)
    measured = dict(native=done)
    verification = dict(status="command_exited_zero", scientific_results_admitted=False,
        source_sha256=module.record(module.ROOT / "benchmark_tools/run_verified_slurm_measurement.py")["sha256"], measurement=measured)
    result = dict(schema=module.executor.SESSION_SCHEMA, status="measurement_returned_pending_independent_review",
        job_id=25000, index=12, cell=run["cell"], plan=plan_ref, amendment=amendment_ref, request=request_ref,
        source=module.record(module.executor.__file__), execution_scope=module.SCOPE, automatic_retry=False,
        next_identity_authorized=False, uncontended_timing=False, wrapper=verification, history=history)
    terminal = dict(source="live_controller", verified=dict(fields=dict(JobState="COMPLETED", ExitCode="0:0",
        Comment=request_ref["sha256"])))
    replayed = dict(native_outcome="exited_zero", native_exit_code=0, measured=measured,
        evidence=[module.record(root / "measurement/done.json")], native_cpu_ids=list(range(52, 84)),
        allocated_placement=dict(stub=True), affinity_observation_statuses=["observed_within_affinity"],
        narrow_flagged_intervals=[])
    resources = dict(primary=dict(wall_seconds=2., cpu_seconds=40., peak_memory_bytes=1000000), primary_scopes=module.SCOPES)
    environment = dict(sampled_environment_evidence_valid=True)
    outputs = dict(schema="native12_composed_output_review_v1", job_id=25000, index=12, request=request_ref,
        plan=plan_ref, amendment=amendment_ref, native_outputs_validated=True, evidence=[], checked_files=[])
    execution = dict(new_sources=[])
    context = ({}, execution, plan, {}, {})
    monkeypatch.setattr(module.executor, "execution_binding", lambda *a: (request, context, history))
    monkeypatch.setattr(module, "verify_terminal", lambda job: deepcopy(terminal))
    monkeypatch.setattr(module, "runtime_review", lambda *a: dict(stub_runtime=True))
    monkeypatch.setattr(module, "replay", lambda *a: deepcopy(replayed))
    monkeypatch.setattr(module, "resource_review", lambda *a: deepcopy(resources))
    monkeypatch.setattr(module, "environment_review", lambda *a: deepcopy(environment))
    monkeypatch.setattr(module, "scientific_output", lambda *a: deepcopy(outputs))
    return dict(root=root, session=session, destination=destination, request_ref=request_ref,
        result=result, verification=verification, terminal=terminal, replayed=replayed, outputs=outputs,
        environment=environment)


def run_join(data):
    store(data["root"] / "verification.json", data["verification"])
    store(data["session"] / "result.json", data["result"])
    return module.read(module.review(data["request_ref"], data["destination"], module.record(module.__file__)["sha256"]))


@pytest.mark.parametrize("kind", ["exited_zero", "exited_nonzero", "timed_out"])
def test_join_retains_actual_outcome_without_scoring_or_retry(joined, kind):
    if kind != "exited_zero":
        joined["terminal"]["verified"]["fields"].update(JobState="FAILED", ExitCode="1:0")
        joined["result"]["status"] = "factorial_attempt_failed_retained"
        joined["verification"]["status"] = "command_failed"
        joined["replayed"].update(native_outcome=kind, native_exit_code=124 if kind == "timed_out" else 1)
    result = run_join(joined)
    assert result["schema"] == module.SCHEMA
    assert result["native_outputs_validated"] is (kind == "exited_zero")
    assert result["terminal_reviewed"] is True
    assert result["primary_resources_replayed"] is True
    assert result["current_original_os_inventory_equality"] is False
    assert all(result[key] is False for key in ("next_identity_authorized", "automatic_retry", "accuracy_evaluated",
        "scientific_timings_admitted", "uncontended_timing", "publication_ready", "original_review_translated"))
    assert not (joined["destination"] / "resource_replay.json").exists()
    summary = module.read(result["resource_replay"])
    assert summary["full_replay_executed"] is True


@pytest.mark.parametrize("mutation", ["session", "source", "amendment", "history", "comment", "collector", "outputs", "runtime"])
def test_failed_join_preserves_failure_without_admission(joined, monkeypatch, mutation):
    if mutation == "session":
        joined["result"]["schema"] = "allocated_native_factorial_session_v1"
    elif mutation == "source":
        joined["result"]["source"] = {}
    elif mutation == "amendment":
        joined["result"]["amendment"] = {}
    elif mutation == "history":
        joined["result"]["history"]["original_review_translated"] = True
    elif mutation == "comment":
        joined["terminal"]["verified"]["fields"]["Comment"] = "wrong"
    elif mutation == "collector":
        joined["replayed"]["measured"] = {}
    elif mutation == "outputs":
        joined["outputs"]["amendment"] = {}
    else:
        def failure(*args):
            raise ValueError("stub fresh runtime refusal")
        monkeypatch.setattr(module, "runtime_review", failure)
    with pytest.raises(ValueError):
        run_join(joined)
    assert not (joined["destination"] / "review.json").exists()
    if mutation != "comment":
        retained = json.loads((joined["destination"] / "failure.json").read_text())
        assert retained["terminal_reviewed"] is False and retained["automatic_retry"] is False


def test_review_batch_envelope_syntax_and_unscheduled_guard(tmp_path):
    batch = module.ROOT / "benchmark_tools/results/native12_composed_terminal_review_20261008_v1.sh"
    raw = batch.read_text()
    for value in ("--cpus-per-task=2", "--mem=128G", "--time=06:00:00", "--no-requeue",
                  "-X faulthandler", "--source-sha256", "PYTHONDONTWRITEBYTECODE=1"):
        assert value in raw
    assert subprocess.run(["bash", "-n", str(batch)], capture_output=True).returncode == 0
    env = dict(os.environ)
    env.pop("SLURM_JOB_ID", None)
    result = subprocess.run(["bash", str(batch)], env=env, cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode != 0 and "SLURM_JOB_ID" in result.stderr
