from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools import review_threadripper_overhead as module
from benchmark_tools.audit_threadripper_overhead import summarize
from benchmark_tools.prepare_threadripper_overhead import build
from tests.unit.test_audit_threadripper_overhead import panel, plan, rows
from tests.unit.test_prepare_threadripper_overhead import parent, roots
from tests.unit.test_verify_threadripper_controller import RAW


def put(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))
    return module.record(path)


@pytest.fixture
def setup(tmp_path, monkeypatch):
    parent_path = Path(__file__).resolve().parents[2] / "benchmark_tools/results/threadripper_private_commands_20260928.json"
    plan = build(parent(), *roots(tmp_path))
    plan.update(sources=[module.record(parent_path)], helpers=[])
    plan_ref = put(tmp_path / "plan.json", plan)
    native_rows = rows(plan)
    for index, row in enumerate(native_rows):
        row.update(job_id=1000 + index, native_outcome="exited_zero",
                   output_checks={"synthetic_output_validation": True})
    attempts, sessions = [], []
    for task in plan["runs"]:
        index, job = task["index"], 1000 + task["index"]
        folder = tmp_path / f"session_{index:02d}"
        raw = RAW.replace("JobId=42", f"JobId={job}").replace("RUNNING", "COMPLETED")
        raw = raw.replace("1-00:00:00", "1-02:00:00")
        scheduler_path = folder / "scheduler.txt"
        folder.mkdir()
        scheduler_path.write_text(raw)
        attempts.append(dict(index=index, job_id=job, scheduler=module.record(scheduler_path)))
        controller = put(folder / "controller.json", dict(
            command=["scontrol", "show", "job", str(job), "--oneliner"], returncode=0, stdout=raw))
        native = put(folder / "native_audit.json", dict(job_id=job, native_outcome="exited_zero",
            status="native_success_outputs_verified", outputs=native_rows[index]["output_checks"],
            replay=dict(native_wall_s=native_rows[index]["native_wall_s"])))
        supporting = put(folder / "synthetic_support.json", dict(synthetic_test_only=True))
        reviews = {k: put(folder / (k + ".json"), dict(schema="threadripper_overhead_review_v1",
            index=index, pair=task["pair"], arm=task["arm"], job_id=job,
            plan_sha256=plan_ref["sha256"], category=k, decision="passed",
            review_reference="Synthetic unit-test review, not actual admission.",
            evidence=[native] if k == "outputs_or_failure" else [supporting])) for k in module.REVIEWS}
        sessions.append(put(folder / "session.json", dict(schema="threadripper_overhead_session_v1",
            index=index, pair=task["pair"], arm=task["arm"], job_id=job, phase="terminal",
            plan_sha256=plan_ref["sha256"], controller=controller, native_audit=native,
            native_outcome="exited_zero", reviews=reviews)))
    attempts_ref = put(tmp_path / "attempts.json", dict(
        schema="threadripper_native_overhead_attempts_v1", plan=plan_ref, automatic_retry=False, attempts=attempts))
    review_ref = put(tmp_path / "reviews.json", dict(schema="threadripper_overhead_reviews_v1",
        plan=plan_ref, attempt_history=attempts_ref, sessions=sessions))
    raw = dict(plan=plan_ref, attempt_history=attempts_ref, runs=native_rows,
               comparison=summarize(plan, native_rows, []), evidence=[plan_ref, attempts_ref])
    calls = []
    def replay(*args):
        calls.append(args)
        return deepcopy(raw)
    monkeypatch.setattr(module, "audit", replay)
    return dict(root=tmp_path, plan=plan, plan_ref=plan_ref, attempts=attempts, sessions=sessions,
                attempts_ref=attempts_ref, reviews_ref=review_ref, raw=raw, calls=calls)


def run(setup):
    return module.review_panel(setup["plan_ref"], setup["attempts_ref"], setup["reviews_ref"],
        baseline_ref={"synthetic": True}, launcher="/python",
        scheduler_command="/recipe/run.sh", scheduler_cwd="/recipe")


def test_complete_conditional_engineering_decision_is_never_production_admission(setup):
    result = run(setup)
    assert result["status"] == "reviewed_engineering_budget_passed"
    assert result["engineering_budget_passed"] is True
    assert result["reviewed_runtime_environment_complete"]
    assert len(result["runs"]) == 54 and len(setup["calls"]) == 1
    assert not result["review_conclusions_independently_certified"]
    for field in ("scientific_timings_admitted", "publication_ready", "next_submission_authorized", "automatic_retry"):
        assert result[field] is False
    assert result["raw_output_audit"]["comparison"]["engineering_budget_passed"] is None
    for ref in result["evidence"]: module.check(ref)


def test_complete_review_keeps_exact_numerical_failure(setup):
    for row in setup["raw"]["runs"]:
        if row["arm"] == "periodic":
            row.update(native_wall_s=106., native_duration_ns=106 * 10**9)
            folder = setup["root"] / f"session_{row['index']:02d}"
            session = json.loads((folder / "session.json").read_text())
            native = json.loads((folder / "native_audit.json").read_text())
            native["replay"]["native_wall_s"] = 106.
            native_ref = put(folder / "native_audit.json", native)
            session["native_audit"] = native_ref
            output_review = json.loads((folder / "outputs_or_failure.json").read_text())
            output_review["evidence"] = [native_ref]
            session["reviews"]["outputs_or_failure"] = put(folder / "outputs_or_failure.json", output_review)
            setup["sessions"][row["index"]] = put(folder / "session.json", session)
    data = json.loads(Path(setup["reviews_ref"]["path"]).read_text())
    data["sessions"] = setup["sessions"]
    setup["reviews_ref"] = put(Path(setup["reviews_ref"]["path"]), data)
    setup["raw"]["comparison"] = summarize(setup["plan"], setup["raw"]["runs"], [])
    result = run(setup)
    assert result["engineering_budget_passed"] is False
    assert result["status"] == "reviewed_engineering_budget_failed"


@pytest.mark.parametrize("defect", ["review_count", "review_plan", "wrong_history", "scheduler",
    "native_outputs", "native_wall", "review_reference", "unbound_native", "production_schema", "raw_identity", "missing_raw"])
def test_mixed_or_unsupported_review_cannot_produce_a_decision(setup, defect):
    review_path = Path(setup["reviews_ref"]["path"])
    data = json.loads(review_path.read_text())
    folder = setup["root"] / "session_53"
    session = json.loads((folder / "session.json").read_text())
    if defect == "review_count": data["sessions"].pop()
    elif defect == "review_plan": data["plan"] = setup["attempts_ref"]
    elif defect == "wrong_history": data["attempt_history"] = setup["plan_ref"]
    elif defect == "scheduler":
        controller = json.loads((folder / "controller.json").read_text())
        controller["stdout"] += "\n"
        session["controller"] = put(folder / "controller.json", controller)
    elif defect in {"native_outputs", "native_wall"}:
        native = json.loads((folder / "native_audit.json").read_text())
        if defect == "native_outputs": native["outputs"] = {"changed": True}
        else: native["replay"]["native_wall_s"] += 1
        native_ref = put(folder / "native_audit.json", native)
        session["native_audit"] = native_ref
        output_review = json.loads((folder / "outputs_or_failure.json").read_text())
        output_review["evidence"] = [native_ref]
        session["reviews"]["outputs_or_failure"] = put(folder / "outputs_or_failure.json", output_review)
    elif defect in {"review_reference", "unbound_native", "production_schema"}:
        reviewed = json.loads((folder / "outputs_or_failure.json").read_text())
        if defect == "review_reference": reviewed["review_reference"] = " "
        elif defect == "unbound_native": reviewed["evidence"] = [setup["plan_ref"]]
        else: reviewed["schema"] = "threadripper_panel_review_v1"
        session["reviews"]["outputs_or_failure"] = put(folder / "outputs_or_failure.json", reviewed)
    elif defect == "raw_identity": setup["raw"]["runs"][0]["pair"] = 100
    else: setup["raw"]["runs"].pop()
    if defect not in {"review_count", "review_plan", "wrong_history", "raw_identity", "missing_raw"}:
        data["sessions"][-1] = put(folder / "session.json", session)
    setup["reviews_ref"] = put(review_path, data)
    with pytest.raises(ValueError): run(setup)


@pytest.mark.parametrize("decision", [None, "failed", "unresolved"])
def test_last_review_incomplete_or_failed_has_no_budget_decision(setup, decision):
    folder = setup["root"] / "session_53"
    session = json.loads((folder / "session.json").read_text())
    if decision is None:
        session["reviews"] = None
    else:
        reviewed = json.loads((folder / "runtime.json").read_text())
        reviewed["decision"] = decision
        session["reviews"]["runtime"] = put(folder / "runtime.json", reviewed)
    data = json.loads(Path(setup["reviews_ref"]["path"]).read_text())
    data["sessions"][-1] = put(folder / "session.json", session)
    setup["reviews_ref"] = put(Path(setup["reviews_ref"]["path"]), data)
    result = run(setup)
    assert result["engineering_budget_passed"] is None
    assert not result["reviewed_runtime_environment_complete"]
    assert result["status"] == "engineering_panel_incomplete_or_ineligible"


def test_incomplete_raw_panel_cannot_be_rescued_by_passing_reviews(setup):
    setup["raw"]["comparison"]["raw_output_panel_complete"] = False
    result = run(setup)
    assert result["engineering_budget_passed"] is None
    assert not result["reviewed_runtime_environment_complete"]


def unrun_inputs(tmp_path, panel):
    results = Path(__file__).resolve().parents[2] / "benchmark_tools/results"
    plan_path, attempts_path, _, attempts = panel
    copied = json.loads(plan_path.read_text())
    parent_path = Path(copied["sources"][0]["path"])
    parent_path.write_bytes((results / "threadripper_private_commands_20260928.json").read_bytes())
    copied["sources"][0] = module.record(parent_path)
    plan_ref = put(plan_path, copied)
    attempts["plan"] = plan_ref
    attempts_ref = put(attempts_path, attempts)
    reviews_ref = put(tmp_path / "empty_reviews.json", dict(schema="threadripper_overhead_reviews_v1",
        plan=plan_ref, attempt_history=attempts_ref, sessions=[]))
    return plan_ref, attempts_ref, reviews_ref


def test_real_unrun_panel_replays_copied_pins_without_native_or_host_work(tmp_path, panel):
    plan_ref, attempts_ref, reviews_ref = unrun_inputs(tmp_path, panel)
    result = module.review_panel(plan_ref, attempts_ref, reviews_ref)
    assert len(result["runs"]) == 54
    assert all(r["raw_status"] == r["review_status"] == "unrun" for r in result["runs"])
    assert result["engineering_budget_passed"] is None
    assert result["bound_review_history"]["progress"]["index"] == 0
    assert not result["raw_output_audit"]["comparison"]["raw_output_panel_complete"]


@pytest.mark.parametrize("fault", [None, "plan_hash", "attempts_hash", "reviews_hash", "existing", "baseline_missing"])
def test_cli_keeps_external_pins_and_never_overwrites(tmp_path, panel, monkeypatch, fault):
    refs = unrun_inputs(tmp_path, panel)
    output = tmp_path / "result.json"
    argv = ["review"]
    for name, ref in zip(("plan", "attempts", "reviews"), refs):
        argv.extend(["--" + name, ref["path"], "--" + name + "-sha256",
                     "0" * 64 if fault == name + "_hash" else ref["sha256"]])
    argv.extend(["--output", str(output)])
    if fault == "existing": output.write_text("Preserve existing output.")
    if fault == "baseline_missing": argv.extend(["--baseline", refs[0]["path"]])
    monkeypatch.setattr("sys.argv", argv)
    if fault:
        with pytest.raises((ValueError, FileExistsError)): module.main()
        if fault == "existing": assert output.read_text() == "Preserve existing output."
        else: assert not output.exists()
    else:
        module.main()
        result = json.loads(output.read_text())
        assert result["engineering_budget_passed"] is None
        assert not result["scientific_timings_admitted"]
