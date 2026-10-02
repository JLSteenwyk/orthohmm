"""Combine native replay and recorded engineering reviews; never launch work."""

import argparse
from pathlib import Path

from benchmark_tools.audit_threadripper_overhead import audit, direct_pin, read_pin, IDENTITY
from benchmark_tools.bind_threadripper_panel_history import bind, REVIEWS
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.verify_lineage_native_provenance import same
from benchmark_tools.run_threadripper_scaling import PRIVATE_PLAN_SHA


def review_panel(plan_ref, attempts_ref, reviews_ref, *, baseline_ref=None,
                 launcher=None, scheduler_command=None, scheduler_cwd=None):
    """Recompute native arithmetic, then bind externally reviewed conclusions.

    Runtime and environment review conclusions are not independently certified
    here. A complete decision is conditional on those recorded reviews.
    """
    plan, attempts, reviews = map(read_pin, (plan_ref, attempts_ref, reviews_ref))
    if plan["sources"][0]["sha256"] != PRIVATE_PLAN_SHA:
        raise ValueError("Require the frozen private scientific parent")
    if (reviews.get("schema") != "threadripper_overhead_reviews_v1"
            or not same(reviews.get("plan"), plan_ref)
            or not same(reviews.get("attempt_history"), attempts_ref)
            or not isinstance(reviews.get("sessions"), list)
            or len(reviews["sessions"]) != len(attempts["attempts"])):
        raise ValueError("Require exactly one engineering review session per recorded attempt")
    raw = audit(Path(plan_ref["path"]), plan_ref["sha256"], Path(attempts_ref["path"]),
                baseline_ref, launcher, scheduler_command, scheduler_cwd)
    if not same(raw["plan"], plan_ref) or not same(raw["attempt_history"], attempts_ref):
        raise ValueError("Raw replay belongs to another plan or attempt history")
    if len(raw["runs"]) != 54:
        raise ValueError("Require all 54 raw outcomes")
    # Empty histories have no allocation envelope to validate.
    history = bind(plan_ref, reviews["sessions"], command=scheduler_command,
                   cwd=scheduler_cwd, time_limit="1-02:00:00", overhead=True)
    evidence = [plan_ref, attempts_ref, reviews_ref, *raw["evidence"], *history["evidence"]]
    outcomes = []
    for task, row in zip(plan["runs"], raw["runs"]):
        if any(not same(row.get(key), task[key]) for key in IDENTITY):
            raise ValueError("Native replay task identity differs")
        outcome = {key: task[key] for key in IDENTITY}
        outcome.update(raw_status=row["status"], review_status="unrun", decisions=None,
                       scientific_timings_admitted=False)
        index = task["index"]
        if index < len(reviews["sessions"]):
            session = read_pin(reviews["sessions"][index])
            attempt = attempts["attempts"][index]
            controller = read_pin(session["controller"])
            scheduler_ref = attempt["scheduler"]
            check(scheduler_ref)
            scheduler = Path(scheduler_ref["path"])
            if scheduler.resolve() != scheduler or scheduler.is_symlink():
                raise ValueError("Require direct scheduler evidence")
            if (controller["stdout"] != scheduler.read_text()
                    or session["job_id"] != attempt["job_id"]):
                raise ValueError("Review controller and raw replay identify different scheduler evidence")
            outcome.update(job_id=attempt["job_id"], session=reviews["sessions"][index],
                           review_status="review_incomplete")
            references = session.get("reviews")
            if references is not None:
                decisions = {}
                for category in sorted(REVIEWS):
                    reviewed = read_pin(references[category])
                    if (not isinstance(reviewed.get("review_reference"), str)
                            or not reviewed["review_reference"].strip()):
                        raise ValueError("Require explicit independent review reference")
                    decisions[category] = reviewed["decision"]
                outcome["decisions"] = decisions
                outcome["review_status"] = "review_failed_or_unresolved"
                if all(value == "passed" for value in decisions.values()):
                    if row["status"] == "native_outputs_replayed":
                        native = read_pin(session["native_audit"])
                        if (native.get("native_outcome") != row["native_outcome"]
                                or not same(native.get("outputs"), row["output_checks"])
                                or not same(native["replay"]["native_wall_s"], row["native_wall_s"])):
                            raise ValueError("Recorded native audit and independent replay disagree")
                        if session["native_audit"] not in read_pin(references["outputs_or_failure"])["evidence"]:
                            raise ValueError("Output review must bind the recorded native audit")
                        outcome["review_status"] = "reviewed_native_success"
                    else:
                        outcome["review_status"] = "raw_outcome_prevents_admission"
        outcomes.append(outcome)
    if len(outcomes) != 54:
        raise ValueError("Require all 54 prescribed outcomes including unrun tasks")
    complete = (history["progress"]["status"] == "all_attempts_reviewed"
                and all(row["review_status"] == "reviewed_native_success" for row in outcomes)
                and raw["comparison"]["raw_output_panel_complete"] is True)
    budget = raw["comparison"]["complete_panel_numerical_budget_passed"] if complete else None
    if complete and type(budget) is not bool:
        raise ValueError("Complete reviewed panel lacks a numerical budget decision")
    for ref in evidence:
        check(ref)
    return dict(schema="threadripper_overhead_panel_review_v1",
        status=("reviewed_engineering_budget_passed" if budget is True
                else "reviewed_engineering_budget_failed" if budget is False
                else "engineering_panel_incomplete_or_ineligible"),
        plan=plan_ref, attempt_history=attempts_ref, review_history=reviews_ref,
        runs=outcomes, raw_output_audit=raw, bound_review_history=history,
        reviewed_runtime_environment_complete=complete,
        engineering_budget_passed=budget, evidence=evidence, source=record(__file__),
        review_conclusions_independently_certified=False,
        scientific_timings_admitted=False, publication_ready=False,
        next_submission_authorized=False, automatic_retry=False,
        limitations=[
            "Native/resource/output replay is recomputed; runtime/environment conclusions are externally recorded reviews, not independently certified here.",
            "A budget decision requires all 54 successful, reviewed and output-consistent tasks; incomplete panels have no passing or failing budget decision.",
            "The complete decision is conditional on the provenance and validity of all four independent review categories.",
            "Engineering-only result: no production identity, causal slowdown certification, timing correction, next-run authorization or publication readiness.",
            "No selective exclusions, partial-cell medians, retries or subtraction of observer overhead."])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("plan", "attempts", "reviews"):
        parser.add_argument("--" + name, type=Path, required=True)
        parser.add_argument("--" + name + "-sha256", required=True)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--baseline-sha256")
    for name in ("launcher", "scheduler-command", "scheduler-cwd"):
        parser.add_argument("--" + name)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    refs = {}
    for name in ("plan", "attempts", "reviews"):
        refs[name] = direct_pin(getattr(args, name))
        if refs[name]["sha256"] != getattr(args, name + "_sha256"):
            raise ValueError("External " + name + " digest differs")
    baseline = None if args.baseline is None else direct_pin(args.baseline)
    if (baseline is None) != (args.baseline_sha256 is None):
        raise ValueError("Baseline path and digest must be supplied together")
    if baseline is not None and baseline["sha256"] != args.baseline_sha256:
        raise ValueError("External baseline digest differs")
    result = review_panel(refs["plan"], refs["attempts"], refs["reviews"], baseline_ref=baseline,
                          launcher=args.launcher, scheduler_command=args.scheduler_command,
                          scheduler_cwd=args.scheduler_cwd)
    save(args.output, result)
    print(result["status"])


if __name__ == "__main__":
    main()
