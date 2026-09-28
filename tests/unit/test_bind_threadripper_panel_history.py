import json

import pytest

from benchmark_tools.bind_threadripper_panel_history import bind, REVIEWS
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.prepare_scaling_inputs import planned_runs
from tests.unit.test_verify_threadripper_controller import RAW


def put(path, data):
    path.write_text(json.dumps(data))
    return record(path)


def setup(tmp_path):
    plan = put(tmp_path / "plan.json", dict(runs=planned_runs()))
    controller = put(tmp_path / "controller.json", dict(
        command=["scontrol", "show", "job", "42", "--oneliner"], returncode=0,
        stdout=RAW.replace("RUNNING", "COMPLETED")))
    support = put(tmp_path / "support.json", dict(diagnostic=True))
    reviews = {k: put(tmp_path / (k+".json"), dict(schema="threadripper_panel_review_v1",
        index=0, job_id=42, plan_sha256=plan["sha256"], category=k,
        decision="passed", evidence=[support])) for k in REVIEWS}
    session = dict(schema="threadripper_panel_session_v1", index=0, job_id=42,
        plan_sha256=plan["sha256"], controller=controller, reviews=reviews,
        phase="terminal", native_outcome="exited_zero")
    session["native_audit"] = put(tmp_path / "native_audit.json", dict(job_id=42,
        native_outcome="exited_zero", status="native_success_outputs_verified"))
    return plan, session


def run(tmp_path, plan, session):
    ref = put(tmp_path / "session.json", session)
    return bind(plan, [ref], command="/recipe/run.sh", cwd="/recipe")


def test_bound_history_advances_but_does_not_authorize(tmp_path):
    plan, session = setup(tmp_path)
    result = run(tmp_path, plan, session)
    assert result["progress"]["index"] == 1
    assert result["scientific_execution_authorized"] is False


def test_missing_review_waits(tmp_path):
    plan, session = setup(tmp_path)
    session["reviews"] = None
    assert run(tmp_path, plan, session)["progress"]["status"] == "terminal_attempt_requires_review"


@pytest.mark.parametrize("field,value", [("job_id", 43), ("index", 1),
    ("category", "wrong"), ("plan_sha256", "wrong"), ("evidence", [])])
def test_rehashed_wrong_review_rejected(tmp_path, field, value):
    plan, session = setup(tmp_path)
    path = tmp_path / "runtime.json"
    data = json.loads(path.read_text())
    data[field] = value
    session["reviews"]["runtime"] = put(path, data)
    with pytest.raises(ValueError):
        run(tmp_path, plan, session)


@pytest.mark.parametrize("name", ["controller.json", "support.json", "runtime.json", "plan.json"])
def test_changed_evidence_rejected(tmp_path, name):
    plan, session = setup(tmp_path)
    (tmp_path / name).write_text("{}")
    with pytest.raises(ValueError):
        run(tmp_path, plan, session)


@pytest.mark.parametrize("field,value", [("returncode", 1), ("returncode", False),
    ("command", ["scontrol", "show", "job", "43", "--oneliner"]),
    ("stdout", RAW.replace("ExitCode=0:0", "ExitCode=1:0").replace("RUNNING", "COMPLETED"))])
def test_controller_contradiction_rejected(tmp_path, field, value):
    plan, session = setup(tmp_path)
    path = tmp_path / "controller.json"
    data = json.loads(path.read_text())
    data[field] = value
    session["controller"] = put(path, data)
    with pytest.raises(ValueError):
        run(tmp_path, plan, session)


@pytest.mark.parametrize("field,value", [("job_id", 43), ("native_outcome", "timed_out"),
                                       ("status", "unverified")])
def test_wrong_native_audit_rejected(tmp_path, field, value):
    plan, session = setup(tmp_path)
    path = tmp_path / "native_audit.json"
    data = json.loads(path.read_text())
    data[field] = value
    session["native_audit"] = put(path, data)
    with pytest.raises(ValueError, match="Native audit"):
        run(tmp_path, plan, session)
