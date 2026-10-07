"""Review terminal allocated-core inference without changing historical reviewers."""

import argparse
import json
from pathlib import Path

from benchmark_tools.native_factorial_allocated_execution import (
    ROOT, SCRIPT, SCOPE, MEMORY, amendment, validate_request, verify_terminal, native_command)
from benchmark_tools.replay_allocated_threadripper_scaling import replay
from benchmark_tools.review_native_factorial_attempt import (
    DISCLOSURE, save, scheduler_fields, runtime_review, resource_review, classify,
    bind_session as historical_bind_session)
from benchmark_tools.derive_threadripper_resources import SCOPES
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_native_factorial_cost import read, available_memory
from benchmark_tools.validate_native_factorial_outputs import Evidence, require
from benchmark_tools.validate_allocated_native_factorial_outputs import validate as validate_outputs
from benchmark_tools.review_threadripper_process_policy import shared_environment, review as review_process_pair
from benchmark_tools.review_threadripper_process_stream import evaluate as process_evaluate
from benchmark_tools.review_threadripper_pressure_stream import evaluate as pressure_evaluate, native_pressure_role
from benchmark_tools.slurm_resource_snapshot import scoped_path
from benchmark_tools.verify_lineage_native_provenance import same
from benchmark_tools.verify_threadripper_controller import remaining_budget


def bind_session(result, verification, request_ref, request, amendment_ref, run):
    historical_bind_session(result, verification, request_ref, request, request["plan"], run)
    require(result.get("schema") == "allocated_native_factorial_session_v1"
        and result.get("amendment") == amendment_ref
        and result.get("source") == record(ROOT / "benchmark_tools/run_allocated_native_factorial_cost.py"),
        "Allocated session source/amendment differs")


def environment_review(request_ref, request, amendment_ref, run, result, directory, done, evidence):
    # Historical environment_review binds the old batch path and remains frozen.
    plan_ref, policy_ref = request["plan"], request["policy"]
    policy = read(policy_ref)
    require(shared_environment(policy) and native_pressure_role(policy) == "diagnostic_only"
        and policy["plan_sha256"] == plan_ref["sha256"], "Wrong shared-host environment policy")
    preflight_path = directory / "environment_preflight.json"
    preflight = evidence.json(preflight_path)
    require(preflight.get("schema") == "threadripper_environment_preflight_v1"
        and preflight.get("decision") == "passed" and preflight.get("job_id") == request["job_id"]
        and preflight.get("index") == run["index"] and preflight.get("environment_policy") == policy_ref
        and preflight.get("plan_sha256") == plan_ref["sha256"] and preflight.get("amendment") == amendment_ref
        and preflight.get("execution_scope") == SCOPE and preflight.get("background_competition_recorded") is True,
        "Environment preflight does not bind this allocated attempt")
    retained_ref = result["environment_review"]
    require(retained_ref["path"] == str(directory / "process_stream_review.json"), "Wrong retained review path")
    retained = read(retained_ref)
    require(retained.get("job_id") == request["job_id"] and retained.get("index") == run["index"]
        and retained.get("execution_scope") == SCOPE and retained.get("uncontended_timing") is False
        and retained.get("background_cpu_used_for_eligibility") is False
        and retained.get("pressure_thresholds_used_for_eligibility") is False
        and retained.get("source") == record(ROOT / "benchmark_tools/review_threadripper_process_stream.py"),
        "Retained environment scope differs")
    process_policy = read(policy["process_policy"])
    refs = [policy_ref, policy["process_policy"], retained_ref,
        *policy["configuration_files"], *policy["evidence"], *preflight["evidence"], *retained["evidence"]]
    for ref in refs:
        check(ref)
        evidence.bind(ref["path"])
    ready = evidence.json(directory / "ready.json")
    scope = scoped_path(ready["cgroup"], request["job_id"])
    job_scope = next(p for p in scope.parents if p.name == f"job_{request['job_id']}")
    stream_path = evidence.bind(directory / "host_processes.jsonl")
    with stream_path.open() as stream:
        first = json.loads(stream.readline())
        observer = first["observer_pid"]
        stream.seek(0)
        processes = process_evaluate(stream, process_policy, boot_id=preflight["boot_id"],
            job_scope=str(job_scope), observer_pid=observer, launch=done["started_ns"] / 1e9,
            end=done["finished_ns"] / 1e9, maximum_foreign_average_cores=policy["maximum_foreign_average_cores"],
            maximum_sample_period_s=policy["maximum_sample_period_s"])
    launch = evidence.json(directory / "launch_environment_observation.json")
    require(len(launch["process_snapshots"]) == 2 and same(launch["process_snapshots"][0], first["snapshot"]),
            "Launch observer snapshot differs")
    launch_processes = review_process_pair(process_policy, *launch["process_snapshots"],
        boot_id=preflight["boot_id"], job_scope=str(job_scope), observer_pid=observer)
    require(same(launch["process_review"], launch_processes) and launch_processes["process_policy_matched"]
        and launch["background_cpu_used_for_eligibility"] is False, "Launch process attribution differs")
    capacity = available_memory(launch["raw_meminfo"])
    require(capacity == launch["available_memory_bytes"] and capacity >= MEMORY
        and policy["minimum_available_memory_bytes"] == MEMORY, "Unsafe or invalid launch capacity")
    budget = evidence.json(directory / "release_budget.json")
    require(budget["status"] == "release_budget_check_passed" and budget["returncode"] == 0,
            "Release budget observation failed")
    expected_budget = remaining_budget(budget["stdout"], request["job_id"], command=str(SCRIPT), cwd=str(ROOT),
        query_elapsed_s=budget["query_elapsed_s"], allocation_mode="shared")
    require(same(expected_budget, budget["budget"])
        and expected_budget["allocation"]["fields"].get("Comment") == request_ref["sha256"]
        and expected_budget["allocation"]["fields"].get("JobName") == "orthohmm_allocated_factorial",
        "Release scheduler identity/budget differs")
    for key, value in processes.items():
        if key not in {"schema", "limitations"}:
            require(same(retained.get(key), value), "Independent process replay differs: " + key)
    paths = sorted(directory.glob("point_*.json"))
    require(paths == [directory / f"point_{i:06d}.json" for i in range(len(paths))], "Pressure inventory differs")
    pressure = pressure_evaluate((evidence.json(p) for p in paths), boot_id=preflight["boot_id"],
        job_scope=str(job_scope), launch_ns=done["started_ns"], end_ns=done["finished_ns"],
        limits=policy["maximum_pressure_percent"], maximum_period_s=policy["maximum_pressure_sample_period_s"],
        pressure_role="diagnostic_only")
    require(sorted(directory.glob("point_*.json")) == paths, "Pressure inventory changed")
    passed = processes["sampled_process_policy_satisfied"] and pressure["sampled_pressure_evidence_satisfied"]
    require(same(retained["pressure_review"], pressure)
        and retained["sampled_environment_policy_satisfied"] is passed
        and result["sampled_environment_evidence_valid"] is passed, "Independent environment verdict differs")
    return dict(status="shared_environment_replayed", sampled_environment_evidence_valid=passed,
        amendment=amendment_ref, processes=processes, pressure=pressure,
        preflight=record(preflight_path), retained=retained_ref, launch_capacity_bytes=capacity,
        launch_processes=launch_processes, release_budget=expected_budget, uncontended_timing=False,
        background_cpu_used_for_eligibility=False, pressure_thresholds_used_for_eligibility=False)


def review(request_ref, destination):
    request = read(request_ref)
    amendment_ref = request["amendment"]
    execution, plan = amendment(amendment_ref)
    plan_ref = execution["historical_plan"]
    validate_request(request, amendment_ref, execution, request["job_id"])
    terminal = verify_terminal(request["job_id"])
    fields = scheduler_fields(terminal)
    if terminal["source"] == "live_controller":
        require(fields.get("Comment") == request_ref["sha256"], "Terminal request comment differs")
    run = plan["runs"][request["index"]]
    baseline = read(plan["baseline"])
    root = Path(run["output_root"])
    session = Path(plan["panel_root"]) / "sessions" / f"run_{run['index']:02d}"
    directory = root / "measurement"
    destination = Path(destination)
    require(destination.is_absolute() and destination.resolve() == destination
        and destination.is_relative_to(ROOT) and not destination.is_relative_to(root)
        and not destination.is_relative_to(session) and not destination.is_relative_to(Path(baseline["core_root"])),
        "Require separate direct review destination outside inference/runtime")
    destination.mkdir(parents=True, exist_ok=False)
    source = record(__file__)
    scheduler_ref = save(destination / "scheduler.json", terminal)
    evidence = Evidence()
    try:
        result = evidence.json(session / "result.json")
        verification = evidence.json(root / "verification.json")
        bind_session(result, verification, request_ref, request, amendment_ref, run)
        for ref in [request_ref, amendment_ref, plan_ref, plan["baseline"], *execution["new_sources"],
                    *plan["helper_sources"], *plan["evidence"]]:
            check(ref)
            evidence.bind(ref["path"])
        runtime = runtime_review(plan, run, session, verification, evidence)
        runtime_ref = save(destination / "runtime.json", runtime)
        replayed = replay(directory, request["job_id"], native_command(amendment_ref, run, baseline))
        require(same(replayed["measured"], verification["measurement"]), "Collector and runtime wrapper differ")
        for ref in replayed["evidence"]:
            check(ref)
            evidence.bind(ref["path"])
        replay_ref = save(destination / "resource_replay.json", replayed)
        done = evidence.json(directory / "done.json")
        status = classify(terminal, result, verification, replayed)
        resources = resource_review(replayed, done, request["job_id"])
        resources_ref = save(destination / "resources.json", resources)
        environment = environment_review(request_ref, request, amendment_ref, run, result, directory, done, evidence)
        environment_ref = save(destination / "environment.json", environment)
        outputs = None
        if status == "native_success":
            outputs = validate_outputs(request_ref)
            require(outputs["job_id"] == request["job_id"] and outputs["index"] == run["index"]
                and outputs["plan"] == plan_ref and outputs["amendment"] == amendment_ref
                and outputs["schema"] == "allocated_native_factorial_output_review_v1"
                and outputs["native_outputs_validated"] is True, "Output review identity differs")
            for ref in [*outputs["evidence"], *outputs["checked_files"]]:
                check(ref)
                evidence.bind(ref["path"])
        outputs_ref = save(destination / "outputs_or_failure.json", outputs if outputs is not None else
            dict(status=status, native_outcome=replayed["native_outcome"], native_exit_code=replayed["native_exit_code"],
                accuracy_evaluated=False, native_outputs_validated=False, automatic_retry=False))
        sampled_valid = environment["sampled_environment_evidence_valid"]
        checked = evidence.finish()
        check(source)
        return save(destination / "review.json", dict(schema="allocated_native_factorial_terminal_review_v1",
            status=status, index=run["index"], job_id=request["job_id"], dataset=run["dataset"],
            cell=run["cell"], repeat=run["repeat"], plan=plan_ref, amendment=amendment_ref,
            request=request_ref, source=source, scheduler=scheduler_ref,
            common_reviewer_source=record(ROOT / "benchmark_tools/review_native_factorial_attempt.py"),
            scheduler_state=fields.get("JobState", fields.get("State")), scheduler_exit_code=fields["ExitCode"],
            reviews=dict(runtime=runtime_ref, resources=resources_ref, environment=environment_ref, outputs_or_failure=outputs_ref),
            resource_replay=replay_ref, evidence=checked, terminal_reviewed=True, next_identity_authorized=sampled_valid,
            resources=resources["primary"], resource_scopes=SCOPES, native_cpu_ids=replayed["native_cpu_ids"],
            allocated_placement=replayed["allocated_placement"],
            shared_host_resources_reviewed=sampled_valid, primary_resources_replayed=True,
            native_outputs_validated=outputs is not None, accuracy_evaluated=False,
            execution_scope=SCOPE, uncontended_timing=False, scientific_timings_admitted=False,
            contention_distortion="unknown_potentially_tool_dependent", timing_disclosure=DISCLOSURE,
            whole_run_maximum_foreign_average_cores=environment["processes"]["maximum_observed_foreign_average_cores"],
            automatic_retry=False, whole_job_cpu_seconds=None, whole_job_peak_memory_bytes=None,
            publication_ready=False, limitations=["New placement review is not benchmark scoring or publication readiness.",
                "Next identity authorization concerns only a different unrun frozen identity after fresh gates.",
                "No retry, background subtraction, fastest-repeat selection or isolated-performance ranking.",
                "CPU includes native wrapper work; memory is native-step lifetime peak, not algorithm-only RSS.",
                "Periodic evidence/runtime brackets do not establish continuous containment or integrity."]))
    except BaseException as error:
        save(destination / "failure.json", dict(schema="allocated_native_factorial_review_failure_v1",
            status="terminal_factorial_review_failed", index=run["index"], job_id=request["job_id"],
            plan=plan_ref, amendment=amendment_ref, request=request_ref, scheduler=scheduler_ref,
            source=source, error_type=type(error).__name__, error=str(error), terminal_reviewed=False,
            next_identity_authorized=False, automatic_retry=False, accuracy_evaluated=False))
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request", type=Path, required=True)
    parser.add_argument("--request-sha256", required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    ref = record(args.request)
    require(ref["sha256"] == args.request_sha256, "Request checksum differs")
    print(json.dumps(review(ref, args.output_directory.absolute()), sort_keys=True))
