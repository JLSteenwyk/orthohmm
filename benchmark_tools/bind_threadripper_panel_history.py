"""Bind recorded controller/review provenance before evaluating panel position."""

import json
from pathlib import Path

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.threadripper_panel_progress import position
from benchmark_tools.verify_threadripper_controller import validate

REVIEWS = {"runtime", "environment", "resources", "outputs_or_failure"}


def bind(plan_ref, session_refs, *, command, cwd):
    """Read immutable history only; caller must freeze arguments and audit evidence."""
    evidence = []

    def read(ref):
        path = Path(ref["path"])
        if not path.is_absolute() or path.resolve() != path:
            raise ValueError("Require direct absolute evidence paths")
        check(ref)
        data = json.loads(path.read_text())
        check(ref)
        evidence.append(ref)
        return data

    plan = read(plan_ref)
    if not isinstance(session_refs, list):
        raise ValueError("Require ordered session references")
    attempts = []
    for index, ref in enumerate(session_refs):
        session = read(ref)
        job = session.get("job_id")
        if (session.get("schema") != "threadripper_panel_session_v1"
                or type(session.get("index")) is not int or session["index"] != index
                or type(job) is not int or job <= 0
                or session.get("plan_sha256") != plan_ref["sha256"]):
            raise ValueError("Session plan, index or job identity differs")
        controller = read(session["controller"])
        expected_command = ["scontrol", "show", "job", str(job), "--oneliner"]
        if (controller.get("command") != expected_command
                or type(controller.get("returncode")) is not int or controller["returncode"] != 0):
            raise ValueError("Controller observation failed or queried another job")
        allocation = validate(controller["stdout"], job, session["phase"],
                              command=command, cwd=cwd)
        outcome = session.get("native_outcome")
        if outcome is not None:
            native = read(session["native_audit"])
            expected_status = ("native_success_outputs_verified" if outcome == "exited_zero"
                               else "native_failure_requires_review")
            if (type(native.get("job_id")) is not int or native["job_id"] != job
                    or native.get("native_outcome") != outcome
                    or native.get("status") != expected_status):
                raise ValueError("Native audit job or outcome differs from session")
        reviews = session.get("reviews")
        decisions = None
        if reviews is not None:
            if not isinstance(reviews, dict) or set(reviews) != REVIEWS:
                raise ValueError("Require all four review references")
            decisions = {}
            for category in sorted(REVIEWS):
                review = read(reviews[category])
                expected = dict(schema="threadripper_panel_review_v1", index=index,
                    job_id=job, plan_sha256=plan_ref["sha256"], category=category)
                if any(type(review.get(k)) is not type(v) or review[k] != v
                       for k, v in expected.items()):
                    raise ValueError("Review belongs to another run, job, plan or category")
                if not isinstance(review.get("evidence"), list) or not review["evidence"]:
                    raise ValueError("Review lacks pinned supporting evidence")
                for supporting in review["evidence"]:
                    check(supporting)
                    evidence.append(supporting)
                decisions[category] = review["decision"]
        attempts.append(dict(index=index, job_id=job,
            scheduler_state=allocation["scheduler_state"], review=decisions,
            native_outcome=outcome))
        if allocation["scheduler_state"] == "COMPLETED" and allocation["scheduler_exit_code"] != "0:0":
            raise ValueError("Completed scheduler record has contradictory exit status")
    progress = position(plan["runs"], attempts)
    for ref in evidence:
        check(ref)
    return dict(status="threadripper_panel_history_bound", progress=progress,
        attempts=attempts, evidence=evidence, source=record(__file__),
        scientific_execution_authorized=False, scientific_timings_admitted=False,
        limitations=["Hashes and identities bind recorded claims; they do not independently validate review conclusions.",
                     "Caller must freeze plan/session references and command/cwd, authenticate observations and audit native outcomes.",
                     "Fresh live-job observations and quiet-host preflight remain mandatory before any submission."])
