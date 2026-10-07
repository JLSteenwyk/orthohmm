from copy import deepcopy
import json
from pathlib import Path
import tempfile

import pytest

from benchmark_tools import review_allocated_native_factorial_attempt as reviewer
from tests.unit import test_review_native_factorial_attempt as historical_tests


@pytest.fixture
def environment(tmp_path):
    data = historical_tests.environment.__wrapped__(tmp_path)
    data["request"]["plan"] = data.pop("plan_ref")
    data["amendment_ref"] = dict(path="/synthetic/amendment", bytes=1, sha256="amendment-digest")
    directory = data["directory"]
    budget_path = directory / "release_budget.json"
    budget = json.loads(budget_path.read_text())
    budget["stdout"] = budget["stdout"].replace(str(historical_tests.module.SCRIPT), str(reviewer.SCRIPT))
    budget["stdout"] += " JobName=orthohmm_allocated_factorial"
    budget["budget"] = reviewer.remaining_budget(budget["stdout"], 42, command=str(reviewer.SCRIPT),
        cwd=str(reviewer.ROOT), query_elapsed_s=.1, allocation_mode="shared")
    budget_path.write_text(json.dumps(budget))
    preflight_path = directory / "environment_preflight.json"
    preflight = json.loads(preflight_path.read_text())
    preflight["amendment"] = data["amendment_ref"]
    preflight["evidence"] = [reviewer.record(r["path"]) for r in preflight["evidence"]]
    preflight_path.write_text(json.dumps(preflight))
    return data


def test_new_environment_preserves_shared_host_acceptance(environment):
    evidence = reviewer.Evidence()
    result = reviewer.environment_review(**environment, evidence=evidence)
    evidence.finish()
    assert result["sampled_environment_evidence_valid"] is True
    assert result["processes"]["maximum_observed_foreign_average_cores"] == 50.
    assert result["pressure"]["diagnostic_thresholds_satisfied"] is False
    assert result["background_cpu_used_for_eligibility"] is False
    assert result["uncontended_timing"] is False
    assert result["amendment"] == environment["amendment_ref"]


@pytest.mark.parametrize("change", ["old_script", "job_name", "comment", "amendment", "capacity", "summary"])
def test_new_environment_rejects_semantically_resealed_tampering(environment, change):
    directory = environment["directory"]
    preflight_path = directory / "environment_preflight.json"
    preflight = json.loads(preflight_path.read_text())
    if change in {"old_script", "job_name", "comment"}:
        path = directory / "release_budget.json"
        budget = json.loads(path.read_text())
        if change == "old_script":
            budget["stdout"] = budget["stdout"].replace(str(reviewer.SCRIPT), str(historical_tests.module.SCRIPT))
        elif change == "job_name":
            budget["stdout"] = budget["stdout"].replace("orthohmm_allocated_factorial", "orthohmm_factorial_cost")
        else:
            budget["stdout"] = budget["stdout"].replace("request-digest", "wrong")
        path.write_text(json.dumps(budget))
    elif change == "amendment":
        preflight["amendment"] = {}
    elif change == "capacity":
        path = directory / "launch_environment_observation.json"
        launch = json.loads(path.read_text())
        launch.update(raw_meminfo="MemAvailable: 1 kB\n", available_memory_bytes=1024)
        path.write_text(json.dumps(launch))
    else:
        path = directory / "process_stream_review.json"
        summary = json.loads(path.read_text())
        summary["maximum_observed_foreign_average_cores"] = 0.
        path.write_text(json.dumps(summary))
        environment["result"]["environment_review"] = reviewer.record(path)
    preflight["evidence"] = [reviewer.record(r["path"]) for r in preflight["evidence"]]
    preflight_path.write_text(json.dumps(preflight))
    with pytest.raises(ValueError):
        reviewer.environment_review(**environment, evidence=reviewer.Evidence())


def store(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))
    return reviewer.record(path)


@pytest.fixture
def joined(monkeypatch):
    with tempfile.TemporaryDirectory(prefix="allocated_native_review_test_", dir=reviewer.ROOT / "benchmarks/work") as tmp:
        base = Path(tmp)
        root, session = base / "panel/run_10", base / "panel/sessions/run_10"
        baseline_ref = store(base / "baseline.json", dict(core_root=str(base / "private_core"),
            tool_entrypoints=dict(orthohmm_python=dict(absolute_path="/synthetic/private/python"))))
        run = dict(index=10, cell="p1_c0_r1", dataset="qfo_corrected", repeat=0, output_root=str(root))
        plan = dict(panel_root=str(base / "panel"), baseline=baseline_ref, runs=[{}]*10+[run],
                    helper_sources=[], evidence=[])
        plan_ref = store(base / "plan.json", plan)
        amendment_ref = store(base / "amendment.json", {})
        execution = dict(historical_plan=plan_ref, new_sources=[])
        request_ref = store(base / "request.json", dict(schema="allocated_native_factorial_request_v1",
            job_id=42, index=10, cell=run["cell"], plan=plan_ref, amendment=amendment_ref))
        done = dict(exit_code=0, timed_out=False, started_ns=10_000_000_000, finished_ns=12_000_000_000)
        store(root / "measurement/done.json", done)
        measured = dict(native=done)
        verification = dict(status="command_exited_zero", scientific_results_admitted=False,
            source_sha256=reviewer.record(reviewer.ROOT / "benchmark_tools/run_verified_slurm_measurement.py")["sha256"],
            measurement=measured)
        result = dict(schema="allocated_native_factorial_session_v1", status="measurement_returned_pending_independent_review",
            job_id=42, index=10, cell=run["cell"], plan=plan_ref, amendment=amendment_ref, request=request_ref,
            source=reviewer.record(reviewer.ROOT / "benchmark_tools/run_allocated_native_factorial_cost.py"),
            execution_scope=reviewer.SCOPE, automatic_retry=False, next_identity_authorized=False,
            uncontended_timing=False, wrapper=verification)
        terminal = dict(source="live_controller", verified=dict(fields=dict(JobState="COMPLETED", ExitCode="0:0",
            Comment=request_ref["sha256"])))
        replayed = dict(native_outcome="exited_zero", native_exit_code=0, measured=measured,
            evidence=[reviewer.record(root / "measurement/done.json")], native_cpu_ids=list(range(52,84)),
            allocated_placement={"fixture": "synthetic"}, affinity_observation_statuses=["observed_within_affinity"],
            narrow_flagged_intervals=[])
        resources = dict(primary=dict(wall_seconds=2.,cpu_seconds=40.,peak_memory_bytes=1000000),primary_scopes=reviewer.SCOPES)
        env = dict(sampled_environment_evidence_valid=True, processes=dict(maximum_observed_foreign_average_cores=99.))
        outputs = dict(schema="allocated_native_factorial_output_review_v1", job_id=42, index=10, plan=plan_ref,
            amendment=amendment_ref, native_outputs_validated=True, evidence=[], checked_files=[])
        monkeypatch.setattr(reviewer, "amendment", lambda ref: (execution, plan))
        monkeypatch.setattr(reviewer, "validate_request", lambda *a: None)
        monkeypatch.setattr(reviewer, "verify_terminal", lambda job: deepcopy(terminal))
        monkeypatch.setattr(reviewer, "runtime_review", lambda *a: dict(status="synthetic_runtime_kernel"))
        monkeypatch.setattr(reviewer, "replay", lambda *a: deepcopy(replayed))
        monkeypatch.setattr(reviewer, "resource_review", lambda *a: deepcopy(resources))
        monkeypatch.setattr(reviewer, "environment_review", lambda *a: deepcopy(env))
        monkeypatch.setattr(reviewer, "validate_outputs", lambda *a: deepcopy(outputs))
        yield dict(base=base, root=root, session=session, request=request_ref, verification=verification,
            result=result, terminal=terminal, replayed=replayed, environment=env, outputs=outputs)


def run_review(data):
    store(data["root"] / "verification.json", data["verification"])
    store(data["session"] / "result.json", data["result"])
    return reviewer.read(reviewer.review(data["request"], data["base"] / "review"))


@pytest.mark.parametrize("kind", ["exited_zero", "exited_nonzero", "timed_out"])
def test_review_preserves_outcome_resources_and_no_accuracy(joined, kind):
    if kind != "exited_zero":
        joined["terminal"]["verified"]["fields"].update(JobState="FAILED", ExitCode="1:0")
        joined["result"]["status"] = "factorial_attempt_failed_retained"
        joined["verification"]["status"] = "command_failed"
        joined["replayed"].update(native_outcome=kind,native_exit_code=124 if kind=="timed_out" else 1)
    result = run_review(joined)
    assert result["schema"] == "allocated_native_factorial_terminal_review_v1"
    assert result["index"] == 10 and result["terminal_reviewed"] is True
    assert result["native_cpu_ids"] == list(range(52,84))
    assert result["native_outputs_validated"] is (kind == "exited_zero")
    assert result["next_identity_authorized"] is True
    assert all(result[k] is False for k in ("accuracy_evaluated", "automatic_retry", "uncontended_timing",
                                           "scientific_timings_admitted", "publication_ready"))
    assert result["whole_run_maximum_foreign_average_cores"] == 99.


def test_invalid_environment_defers_only_next_launch(joined):
    joined["environment"]["sampled_environment_evidence_valid"] = False
    result = run_review(joined)
    assert result["terminal_reviewed"] is True and result["primary_resources_replayed"] is True
    assert result["next_identity_authorized"] is False


@pytest.mark.parametrize("change", ["old_session", "amendment", "source", "collector", "outputs", "runtime"])
def test_failed_review_retains_failure_no_success_or_retry(joined, monkeypatch, change):
    if change == "old_session": joined["result"]["schema"] = "historical"
    elif change == "amendment": joined["result"]["amendment"] = {}
    elif change == "source": joined["result"]["source"] = {}
    elif change == "collector": joined["replayed"]["measured"] = {"wrong": True}
    elif change == "outputs": joined["outputs"]["amendment"] = {}
    else:
        def fail(*a): raise ValueError("synthetic drift")
        monkeypatch.setattr(reviewer, "runtime_review", fail)
    with pytest.raises(ValueError): run_review(joined)
    failure = json.loads((joined["base"] / "review/failure.json").read_text())
    assert all(failure[k] is False for k in ("terminal_reviewed","next_identity_authorized","automatic_retry"))
    assert not (joined["base"] / "review/review.json").exists()
