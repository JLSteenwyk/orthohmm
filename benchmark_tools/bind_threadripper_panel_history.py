"""Bind recorded controller/review provenance before evaluating panel position."""

import json
from pathlib import Path

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.threadripper_panel_progress import overhead_position, position, PRE_NATIVE_RESOLUTIONS
from benchmark_tools.verify_threadripper_controller import validate

REVIEWS = {"runtime", "environment", "resources", "outputs_or_failure"}


def bind(plan_ref, session_refs, *, command, cwd, time_limit="1-00:00:00", overhead=False,
         allocation_mode="exclusive"):
    """Read immutable history only; caller must freeze arguments and audit evidence."""
    if type(overhead) is not bool:
        raise ValueError("Require explicit panel kind")
    if allocation_mode not in {"exclusive", "shared"}:
        raise ValueError("Unsupported allocation mode")
    session_schema = "threadripper_overhead_session_v1" if overhead else "threadripper_panel_session_v1"
    review_schema = "threadripper_overhead_review_v1" if overhead else "threadripper_panel_review_v1"
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
    if not isinstance(session_refs, list) or len(session_refs) > len(plan["runs"]):
        raise ValueError("Require ordered session references")
    attempts = []
    for index, ref in enumerate(session_refs):
        session = read(ref)
        job = session.get("job_id")
        if (session.get("schema") != session_schema
                or type(session.get("index")) is not int or session["index"] != index
                or type(job) is not int or job <= 0
                or session.get("plan_sha256") != plan_ref["sha256"]):
            raise ValueError("Session plan, index or job identity differs")
        if overhead:
            task = plan["runs"][index]
            if any(type(session.get(k)) is not type(task[k]) or session[k] != task[k]
                   for k in ("pair", "arm")):
                raise ValueError("Engineering session pair or arm differs")
        controller = read(session["controller"])
        expected_command = ["scontrol", "show", "job", str(job), "--oneliner"]
        if (controller.get("command") != expected_command
                or type(controller.get("returncode")) is not int or controller["returncode"] != 0):
            raise ValueError("Controller observation failed or queried another job")
        allocation = validate(controller["stdout"], job, session["phase"],
                              command=command, cwd=cwd, time_limit=time_limit, allocation_mode=allocation_mode)
        outcome = session.get("native_outcome")
        if outcome is not None:
            if outcome == "not_started":
                if overhead or allocation_mode != "shared" or "native_audit" in session:
                    raise ValueError("Pre-native abort requires explicit shared production evidence")
                native = read(session["pre_native_audit"])
                expected_status = "pre_native_infrastructure_failure_reviewed"
                if (native.get("index") != index or type(native.get("index")) is not int
                        or native.get("execution_scope") != "shared_host_matched_resources"
                        or native.get("resources") is not None
                        or native.get("comparative_timing_eligible") is not False
                        or native.get("automatic_retry") is not False
                        or native.get("scientific_timings_admitted") is not False):
                    raise ValueError("Pre-native review has contradictory identity or resource claims")
                check(native["source"])
                if not isinstance(native.get("evidence"), list) or not native["evidence"]:
                    raise ValueError("Pre-native review lacks pinned evidence")
                for pin in native["evidence"]:
                    check(pin)
                    evidence.append(pin)
            else:
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
                expected = dict(schema=review_schema, index=index,
                    job_id=job, plan_sha256=plan_ref["sha256"], category=category)
                if overhead:
                    expected.update(pair=task["pair"], arm=task["arm"])
                if any(type(review.get(k)) is not type(v) or review[k] != v
                       for k, v in expected.items()):
                    raise ValueError("Review belongs to another run, job, plan or category")
                if not isinstance(review.get("evidence"), list) or not review["evidence"]:
                    raise ValueError("Review lacks pinned supporting evidence")
                for supporting in review["evidence"]:
                    check(supporting)
                    evidence.append(supporting)
                decisions[category] = review["decision"]
        attempt = dict(index=index, job_id=job,
            scheduler_state=allocation["scheduler_state"], review=decisions,
            native_outcome=outcome)
        if "resolution" in session:
            if overhead or allocation_mode != "shared":
                raise ValueError("Monitoring resolution requires the shared production panel")
            resolution = read(session["resolution"])
            if resolution.get("plan_sha256") != plan_ref["sha256"]:
                raise ValueError("Monitoring resolution belongs to another plan")
            original = read(resolution["original_session"])
            if original != {key: value for key, value in session.items() if key != "resolution"}:
                raise ValueError("Resolution changed the original attempt or review decisions")
            if resolution.get("schema") == "threadripper_pre_native_failure_resolution_v1":
                kind = resolution.get("kind")
                if type(kind) is not str or kind not in PRE_NATIVE_RESOLUTIONS:
                    raise ValueError("Pre-native resolution has an unsupported failure kind")
                failure_kind, validation_status = PRE_NATIVE_RESOLUTIONS[kind]
                if (outcome != "not_started" or resolution.get("pre_native_audit") != session.get("pre_native_audit")
                        or resolution.get("environment_review") != reviews["environment"]
                        or native.get("failure_kind") != failure_kind
                        or native.get("scheduler_state") != allocation["scheduler_state"]
                        or native.get("scheduler_exit_code") != allocation["scheduler_exit_code"]):
                    raise ValueError("Pre-native resolution borrows or changes the original failure")
                supporting = resolution.get("evidence")
                if not isinstance(supporting, list) or not supporting or not resolution.get("review_reference"):
                    raise ValueError("Pre-native resolution lacks explicit supporting evidence")
                for pin in supporting:
                    check(pin)
                    evidence.append(pin)
                repair = read(resolution["repair_validation"])
                if (repair.get("status") != validation_status
                        or repair.get("returncode") != 0 or type(repair.get("returncode")) is not int
                        or not repair.get("evidence")):
                    raise ValueError("Pre-native resolution lacks successful repair validation")
                for pin in repair["evidence"]:
                    check(pin)
                    evidence.append(pin)
                pinned = {Path(pin["path"]).name: pin for pin in native["evidence"]}
                gate = read(pinned["go.json"])
                aborted = read(pinned["aborted_before_native.json"])
                release = read(pinned["environment_release.json"])
                preflight = read(pinned["environment_preflight.json"])
                if (gate != {"abort": True} or aborted != {"status": "observer_did_not_release_native"}
                        or release.get("status") != "environment_release_failed"
                        or preflight.get("decision") != "failed"
                        or preflight.get("index") != index or preflight.get("job_id") != job):
                    raise ValueError("Pre-native resolution lacks the actual aborted release evidence")
                if kind == "pre_native_environment_response_deadline" and (
                        preflight.get("publication_deadline_expired") is not True
                        or release.get("error_type") != "TimeoutError"
                        or release.get("error") != "Environmental worker response deadline exceeded"):
                    raise ValueError("Deadline resolution lacks the actual expired response evidence")
                attempt.update(resolution=resolution, scheduler_exit_code=allocation["scheduler_exit_code"])
                attempts.append(attempt)
                continue
            if (resolution["native_audit"] != session["native_audit"]
                    or resolution["environment_review"] != reviews["environment"]):
                raise ValueError("Resolution borrows another native or environment review")
            supporting = resolution.get("evidence")
            if not isinstance(supporting, list) or not supporting or not resolution.get("review_reference"):
                raise ValueError("Resolution lacks explicit supporting evidence")
            for pin in supporting:
                check(pin)
                evidence.append(pin)
            replay = read(resolution["environment_replay"])
            env_review = read(reviews["environment"])
            if resolution["environment_replay"] not in env_review["evidence"]:
                raise ValueError("Resolution replay was not bound by the original independent review")
            process, pressure = replay["processes"], replay["pressure"]
            failures = process["failures"]
            if (set(failures) != {"sample_period_bound_exceeded"}
                    or type(failures["sample_period_bound_exceeded"]) is not int
                    or failures["sample_period_bound_exceeded"] <= 0
                    or process["sampled_process_policy_satisfied"] is not False
                    or process["execution_scope"] != "shared_host_matched_resources"
                    or process["job_scope"].split("/")[-1] != f"job_{job}"
                    or process["command_bracketed_by_samples"] is not True
                    or process["intervals"] != process["records"] - 1
                    or process["policy_matched_intervals"] != process["intervals"]
                    or pressure["sampled_pressure_evidence_satisfied"] is not True
                    or pressure["failures"]):
                raise ValueError("Resolution cannot waive identity, bracketing, pressure or accounting failures")
            result = read(resolution["executor_result"])
            if (result.get("status") != "executor_failed" or result.get("job_id") != job
                    or result.get("index") != index
                    or result.get("error") != "Whole-run sampled environment policy was not satisfied"
                    or result["wrapper"]["status"] != "command_exited_zero"
                    or result["wrapper"]["measurement"]["native"]["exit_code"] != 0
                    or result["wrapper"]["measurement"]["native"]["timed_out"] is not False):
                raise ValueError("Resolution does not explain the actual post-native executor failure")
            attempt.update(resolution=resolution, scheduler_exit_code=allocation["scheduler_exit_code"])
        attempts.append(attempt)
        if allocation["scheduler_state"] == "COMPLETED" and allocation["scheduler_exit_code"] != "0:0":
            raise ValueError("Completed scheduler record has contradictory exit status")
    extra = dict(allow_monitoring_resolution=allocation_mode == "shared") if not overhead else {}
    progress = (overhead_position if overhead else position)(plan["runs"], attempts, **extra)
    for ref in evidence:
        check(ref)
    return dict(status="threadripper_panel_history_bound", progress=progress,
        controller_policy=dict(command=command, cwd=cwd, time_limit=time_limit),
        attempts=attempts, evidence=evidence, source=record(__file__),
        scientific_execution_authorized=False, scientific_timings_admitted=False,
        limitations=["Hashes and identities bind recorded claims; they do not independently validate review conclusions.",
                     "Caller must freeze plan/session references and command/cwd/time limit, authenticate observations and audit native outcomes.",
                     "Fresh live-job observations and the selected environmental preflight remain mandatory before native release."])
