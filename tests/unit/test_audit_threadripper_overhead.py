from copy import deepcopy
import json
from pathlib import Path
import sys

import pytest

from benchmark_tools import audit_threadripper_overhead as module
from benchmark_tools.prepare_threadripper_overhead import build
from tests.unit.test_prepare_threadripper_overhead import parent, roots
from tests.unit.test_replay_threadripper_boundary import archive, evidence
from tests.unit.test_verify_threadripper_controller import RAW


@pytest.fixture
def plan(tmp_path):
    return build(parent(), *roots(tmp_path))


def rows(plan, *, status="native_outputs_replayed", periodic_ns=103 * 10**9):
    result = []
    for task in plan["runs"]:
        row = {k: task[k] for k in module.IDENTITY}
        duration = periodic_ns if task["arm"] == "periodic" else 100 * 10**9
        row.update(status=status, native_wall_s=duration / 1e9, native_duration_ns=duration,
            work_identity=dict(method=task["method"], proteomes=task["proteomes"]),
            scientific_timings_admitted=False)
        result.append(row)
    return result


def test_complete_design_and_exact_numerical_budget_are_not_full_admission(plan):
    module.design(plan, parent())
    result = module.summarize(plan, rows(plan), [])
    assert len(result["pairs"]) == 27 and len(result["cells"]) == 9
    assert result["raw_output_panel_complete"] is True
    assert result["complete_panel_numerical_budget_passed"] is True
    assert all(cell["median_signed_ratio"] == .03 for cell in result["cells"])
    assert result["engineering_budget_passed"] is None
    assert result["runtime_environment_admission_complete"] is False
    assert result["scientific_timings_admitted"] is False


@pytest.mark.parametrize("periodic_ns,passed", [(105 * 10**9, True), (105 * 10**9 + 1, False)])
def test_exact_median_boundary_not_float_rounding_tolerance(plan, periodic_ns, passed):
    result = module.summarize(plan, rows(plan, periodic_ns=periodic_ns), [])
    assert result["complete_panel_numerical_budget_passed"] is passed


@pytest.mark.parametrize("extra,passed", [(0, True), (1, False)])
def test_exact_per_pair_budget_boundary(plan, extra, passed):
    outcomes = rows(plan, periodic_ns=100 * 10**9)
    target = next(r for r in outcomes[:2] if r["arm"] == "periodic")
    target.update(native_duration_ns=110 * 10**9 + extra, native_wall_s=110 + extra / 1e9)
    result = module.summarize(plan, outcomes, [])
    assert result["pairs"][0]["numerical_pair_budget_passed"] is passed
    assert result["complete_panel_numerical_budget_passed"] is passed


def test_negative_ratios_are_retained(plan):
    result = module.summarize(plan, rows(plan, periodic_ns=90 * 10**9), [])
    assert all(p["signed_ratio"] == -.1 for p in result["pairs"])
    assert result["complete_panel_numerical_budget_passed"] is True


@pytest.mark.parametrize("status", sorted(module.STATUSES - {"native_outputs_replayed"}))
def test_any_failed_unrun_or_nonterminal_task_prevents_complete_median(plan, status):
    outcomes = rows(plan)
    outcomes[0]["status"] = status
    result = module.summarize(plan, outcomes, [])
    assert result["pairs"][0]["signed_ratio"] is None
    affected = next(c for c in result["cells"] if 0 in c["pair_indices"])
    assert affected["consistent_pairs"] == 2 and affected["median_signed_ratio"] is None
    assert result["complete_panel_numerical_budget_passed"] is None


@pytest.mark.parametrize("change", ["missing", "order", "bool_index", "admitted", "unknown", "identity"])
def test_outcome_contract_rejects_drift(plan, change):
    outcomes = rows(plan)
    if change == "missing": outcomes.pop()
    elif change == "order": outcomes.reverse()
    elif change == "bool_index": outcomes[0]["index"] = False
    elif change == "admitted": outcomes[0]["scientific_timings_admitted"] = True
    elif change == "unknown": outcomes[0]["status"] = "validated"
    else: outcomes[0]["pair"] = 27
    with pytest.raises(ValueError):
        module.summarize(plan, outcomes, [])


@pytest.mark.parametrize("value", [0, -1, True, 10**15])
def test_invalid_integer_clock_duration_cannot_form_a_pair(plan, value):
    outcomes = rows(plan)
    outcomes[0]["native_duration_ns"] = value
    result = module.summarize(plan, outcomes, [])
    assert result["pairs"][0]["signed_ratio"] is None
    assert result["complete_panel_numerical_budget_passed"] is None


def test_mismatched_outputs_and_panel_issues_prevent_panel_conclusion(plan):
    outcomes = rows(plan)
    outcomes[0]["work_identity"] = {"changed": True}
    result = module.summarize(plan, outcomes, [])
    assert result["pairs"][0]["status"] == "output_mismatch"
    assert result["complete_panel_numerical_budget_passed"] is None
    result = module.summarize(plan, rows(plan), [{"stage": "overlap"}])
    assert result["raw_output_panel_complete"] is False
    assert result["complete_panel_numerical_budget_passed"] is None


@pytest.mark.parametrize("change", ["seed", "budget", "resources", "status", "authorized",
    "statistic", "host", "missing", "order", "native", "input", "environment"])
def test_frozen_design_drift_is_rejected(plan, change):
    if change == "seed": plan["ordering_seed"] = "other"
    elif change == "budget": plan["engineering_budget"]["every_pair_max"] = .2
    elif change == "resources": plan["resources"]["native_workers"] = 20
    elif change == "status": plan["status"] = "completed"
    elif change == "authorized": plan["execution_authorized"] = True
    elif change == "statistic": plan["statistic"] = "absolute difference"
    elif change == "host": plan["common_host_observation"]["required_in_both_arms"] = False
    elif change == "missing": plan["runs"].pop()
    elif change == "order": plan["runs"].reverse()
    elif change == "native": plan["runs"][0]["run"]["native_argv"].append("--different")
    elif change == "input": plan["runs"][0]["run"]["dataset"]["inputs"][0]["bytes"] += 1
    else: plan["environment_overrides"]["NEW_VARIABLE"] = "changed"
    with pytest.raises(ValueError):
        module.design(plan, parent())


@pytest.mark.parametrize("code", [0, 7])
def test_boundary_task_composes_actual_raw_replay_and_existing_cpu_reader(archive, monkeypatch, code):
    import benchmark_tools.audit_threadripper_native_outcome as native_module
    (archive / "native.log").write_text("")
    report = archive / "boundary_report.json"
    measured = json.loads(report.read_text())
    if code:
        measured["native"]["exit_code"] = code
        measured["status"] = "command_failed"
        (archive / "done.json").write_text(json.dumps(measured["native"]))
        report.write_text(json.dumps(measured))
    output = archive / "output"
    output.mkdir()
    run = dict(native_argv=["/usr/bin/true"], cwd=str(archive), measurement_directory=str(archive),
        configuration=dict(output=str(output)))
    task = dict(index=0, pair=0, arm="boundary", method="orthohmm_high_sensitivity", proteomes=4, repeat=0, run=run)
    def validate(*args):
        assert not code
        return dict(status="threadripper_native_outputs_checked", native=dict(checked_files=[]))
    def native(run, baseline, job, **kwargs):
        return native_module.audit(run, baseline, job, validate_fn=validate, **kwargs)
    monkeypatch.setattr(module, "audit_native", native)
    monkeypatch.setattr(module, "fingerprint", lambda *a: dict(identity={"synthetic": True}, evidence=[],
        helpers=[], source=module.record(__file__)))
    result = module.audit_task(task, {}, 21816, sys.executable)
    assert result["status"] == ("native_failure" if code else "native_outputs_replayed")
    assert result["cpu"]["cpu_seconds"] > 0 and result["step_memory"]["bytes"] == 200
    assert result["scientific_timings_admitted"] is False
    with pytest.raises(ValueError):
        module.audit_task(task, {}, 21816, "/different/launcher")


@pytest.fixture
def panel(tmp_path, plan):
    actual = Path(module.__file__).parent / "results/threadripper_native_overhead_plan_20260930.json"
    retained = json.loads(actual.read_text())
    checkout = Path(module.__file__).resolve().parents[1]
    historical_root = Path(retained["sources"][0]["path"]).parents[2]

    def copied_pin(pin):
        relative = Path(pin["path"]).relative_to(historical_root)
        source = checkout / relative
        current = module.record(source)
        assert {k: current[k] for k in ("bytes", "sha256")} == {
            k: pin[k] for k in ("bytes", "sha256")}
        target = tmp_path / "evidence" / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(source.read_bytes())
        copied = module.record(target)
        assert {k: copied[k] for k in ("bytes", "sha256")} == {
            k: pin[k] for k in ("bytes", "sha256")}
        return copied

    source = tmp_path / "parent.json"
    source.write_text(json.dumps(parent()))
    plan["sources"] = [module.record(source), *[copied_pin(pin) for pin in retained["sources"][1:]]]
    plan["helpers"] = [copied_pin(pin) for pin in retained["helpers"]]
    plan_path = tmp_path / "plan.json"
    plan_path.write_text(json.dumps(plan))
    history_path = tmp_path / "attempts.json"
    history = dict(schema="threadripper_native_overhead_attempts_v1", plan=module.record(plan_path),
        attempts=[], automatic_retry=False)
    history_path.write_text(json.dumps(history))
    baseline = tmp_path / "baseline.json"
    baseline.write_text("{}")
    return plan_path, history_path, module.record(baseline), history


def audit(panel):
    plan, history, baseline, _ = panel
    return module.audit(plan, module.record(plan)["sha256"], history, baseline,
        "/private/python", "/recipe/run.sh", "/recipe")


def test_real_plan_complete_unrun_inventory(panel):
    result = audit(panel)
    assert len(result["runs"]) == 54 and all(r["status"] == "unrun" for r in result["runs"])
    assert result["comparison"]["complete_panel_numerical_budget_passed"] is None
    assert all(c["median_signed_ratio"] is None for c in result["comparison"]["cells"])
    unbound = module.audit(panel[0], module.record(panel[0])["sha256"], panel[1], None, None, None, None)
    assert unbound["external_bindings"] == dict(launcher=None, scheduler_command=None, scheduler_cwd=None)


def test_audit_reads_copied_pins_not_workstation_paths(panel, monkeypatch):
    plan = json.loads(panel[0].read_text())
    evidence_root = panel[0].parent / "evidence"
    references = [*plan["sources"][1:], *plan["helpers"]]
    assert len(references) == 8
    assert all(Path(pin["path"]).is_relative_to(evidence_root) for pin in references)
    original_check = module.check
    checked = []

    def check_local(pin):
        assert Path(pin["path"]).is_relative_to(panel[0].parent)
        checked.append(pin)
        return original_check(pin)

    monkeypatch.setattr(module, "check", check_local)
    result = audit(panel)
    assert all(pin in checked for pin in references)
    assert all(row["status"] == "unrun" for row in result["runs"])
    assert result["comparison"]["engineering_budget_passed"] is None
    assert result["scientific_timings_admitted"] is False


@pytest.mark.parametrize("kind", ["protocol", "helper"])
def test_changed_copied_evidence_is_not_admitted(panel, kind):
    plan = json.loads(panel[0].read_text())
    pin = plan["sources"][1] if kind == "protocol" else plan["helpers"][0]
    Path(pin["path"]).write_text("Changed copied evidence.\n")
    with pytest.raises(ValueError):
        audit(panel)


@pytest.mark.parametrize("state,expected", [("FAILED", "scheduler_failed"),
    ("TIMEOUT", "scheduler_failed"), ("RUNNING", "nonterminal_requires_live_check")])
def test_scheduler_failures_and_live_records_are_not_inspected_or_retried(panel, tmp_path, monkeypatch, state, expected):
    scheduler = tmp_path / "scheduler.txt"
    scheduler.write_text(RAW.replace("RUNNING", state).replace("1-00:00:00", "1-02:00:00"))
    history = panel[3]
    history["attempts"] = [dict(index=0, job_id=42, scheduler=module.record(scheduler))]
    panel[1].write_text(json.dumps(history))
    monkeypatch.setattr(module, "audit_task", lambda *a: pytest.fail("Must not inspect non-successful native outputs"))
    result = audit(panel)
    assert result["runs"][0]["status"] == expected
    assert result["comparison"]["complete_panel_numerical_budget_passed"] is None


@pytest.mark.parametrize("fault", ["skip", "retry", "duplicate_job", "bool_job", "history_plan", "history_schema", "retry_allowed"])
def test_history_cannot_select_or_retry_attempts(panel, tmp_path, fault):
    history = panel[3]
    attempt = dict(index=0, job_id=42, scheduler={})
    history["attempts"] = [attempt]
    if fault == "skip": attempt["index"] = 1
    elif fault == "retry": history["attempts"].append(deepcopy(attempt))
    elif fault == "duplicate_job": history["attempts"].append(dict(attempt, index=1))
    elif fault == "bool_job": attempt["job_id"] = True
    elif fault == "history_plan": history["plan"]["sha256"] = "0" * 64
    elif fault == "history_schema": history["schema"] = "other"
    else: history["automatic_retry"] = True
    panel[1].write_text(json.dumps(history))
    with pytest.raises(ValueError):
        audit(panel)


def test_unrecorded_artifacts_and_evidence_mutation_fail_closed(panel, monkeypatch):
    plan = json.loads(panel[0].read_text())
    root, _ = module.paths(plan["runs"][0]["run"])
    root.mkdir(parents=True)
    result = audit(panel)
    assert result["runs"][0]["status"] == "invalid_evidence"
    assert result["comparison"]["panel_issues"]
    original = module.summarize
    def mutate(*args):
        result = original(*args)
        panel[1].write_text("changed")
        return result
    monkeypatch.setattr(module, "summarize", mutate)
    with pytest.raises(ValueError):
        audit(panel)


def test_new_unrun_artifacts_during_audit_are_not_a_stable_empty_panel(panel, monkeypatch):
    plan = json.loads(panel[0].read_text())
    root, _ = module.paths(plan["runs"][0]["run"])
    original = module.summarize
    def mutate(*args):
        result = original(*args)
        root.mkdir(parents=True)
        return result
    monkeypatch.setattr(module, "summarize", mutate)
    with pytest.raises(ValueError, match="inventory changed"):
        audit(panel)


@pytest.mark.parametrize("fault", ["none", "mismatch", "overlap", "boot"])
def test_attempt_prefix_pair_binding_and_stopping(panel, tmp_path, monkeypatch, fault):
    plan = json.loads(panel[0].read_text())
    history = panel[3]
    for index in range(3 if fault == "mismatch" else 2):
        path = tmp_path / f"scheduler-{index}.txt"
        path.write_text(RAW.replace("RUNNING", "COMPLETED").replace("JobId=42", f"JobId={42+index}")
            .replace("1-00:00:00", "1-02:00:00"))
        history["attempts"].append(dict(index=index, job_id=42+index, scheduler=module.record(path)))
    panel[1].write_text(json.dumps(history))
    calls = []
    def native(task, *args):
        index = task["index"]
        calls.append(index)
        directory = Path(task["run"]["measurement_directory"])
        directory.mkdir(parents=True)
        row = rows(plan)[index]
        if fault == "mismatch" and index == 1: row["work_identity"] = {"different": True}
        start = 10**12 + index * 2 * 10**11
        if fault == "overlap" and index == 1: start = 10**12
        done = dict(snapshots=[dict(started_monotonic_ns=start, raw=dict(boot_id="changed" if fault == "boot" and index else "same")),
            dict(finished_monotonic_ns=start + row["native_duration_ns"] + 20)])
        (directory / "done.json").write_text(json.dumps(done))
        row["evidence"] = [module.record(directory / "done.json")]
        return row
    monkeypatch.setattr(module, "audit_task", native)
    result = audit(panel)
    assert calls == [0, 1]
    assert result["comparison"]["engineering_budget_passed"] is None
    if fault == "none":
        assert result["comparison"]["pairs"][0]["signed_ratio"] == .03
    elif fault == "mismatch":
        assert result["comparison"]["pairs"][0]["status"] == "output_mismatch"
        assert result["runs"][2]["status"] == "invalid_evidence"
    else:
        assert result["runs"][1]["status"] == "invalid_evidence"
    assert result["comparison"]["complete_panel_numerical_budget_passed"] is None
