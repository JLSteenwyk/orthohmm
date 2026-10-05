"""Join terminal factorial evidence; never submit, retry or alter an attempt."""

import argparse
import json
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.derive_threadripper_resources import endpoints, SCOPES
from benchmark_tools.inspect_native_python_lookup import compare_lookup
from benchmark_tools.measure_native_scaling_run import check_manifests
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save as persist
from benchmark_tools.replay_threadripper_scaling import replay
from benchmark_tools.review_threadripper_process_policy import shared_environment, review as review_process_pair
from benchmark_tools.review_threadripper_process_stream import evaluate as process_evaluate
from benchmark_tools.review_threadripper_pressure_stream import evaluate as pressure_evaluate, native_pressure_role
from benchmark_tools.run_native_factorial_cost import ROOT, SCOPE, SCRIPT, MEMORY, available_memory, read, validate_plan, validate_request, verify_terminal
from benchmark_tools.slurm_resource_snapshot import scoped_path
from benchmark_tools.validate_native_factorial_outputs import Evidence, require, validate as validate_outputs
from benchmark_tools.verify_lineage_native_provenance import same
from benchmark_tools.verify_threadripper_controller import remaining_budget


DISCLOSURE = ("Shared Threadripper measurements while other analyses were running. CPU, memory-bandwidth "
    "and I/O contention may affect elapsed times by an unknown, potentially tool-dependent amount. "
    "These are observed shared-host timings, not estimates of isolated performance.")


def save(path, value):
    persist(path, value)
    return record(path)


def scheduler_fields(terminal):
    verified = terminal["verified"]
    return verified.get("fields", verified)


def native_command(plan_ref, run, baseline):
    return [baseline["tool_entrypoints"]["orthohmm_python"]["absolute_path"], "-B",
        str(ROOT / "benchmark_tools/run_native_factorial_cost.py"), "--native", "--plan", plan_ref["path"],
        "--plan-sha256", plan_ref["sha256"], "--index", str(run["index"])]


def bind_session(result, verification, request_ref, request, plan_ref, run):
    require(result.get("job_id") == request["job_id"] and result.get("index") == run["index"]
        and result.get("cell") == run["cell"] and result.get("plan") == plan_ref
        and result.get("request") == request_ref and result.get("execution_scope") == SCOPE
        and result.get("automatic_retry") is False and result.get("next_identity_authorized") is False
        and result.get("uncontended_timing") is False and same(result.get("wrapper"), verification),
        "Terminal session and retained wrapper/request differ")
    require(verification.get("scientific_results_admitted") is False
        and verification.get("source_sha256") == record(ROOT / "benchmark_tools/run_verified_slurm_measurement.py")["sha256"],
        "Runtime wrapper source or admission scope differs")


def runtime_review(plan, run, session, verification, evidence, *, tree_checker=check_manifests):
    lookup = read(plan["runtime_lookup"])
    binding = read(lookup["binding"])
    require(lookup["baseline"] == plan["baseline"]
        and lookup["source"] == record(ROOT / "benchmark_tools/inspect_native_python_lookup.py"),
        "Runtime lookup baseline or inspector source differs")
    manifests, expected = [], []
    for path, digest in binding["runtime_specs"]:
        ref = record(path)
        require(ref["sha256"] == digest, "Runtime inventory checksum differs")
        manifest = read(ref)
        evidence.bind(path)
        manifests.append(ref)
        expected.append(dict(path=path, sha256=digest, status="runtime_tree_identity_matches",
            records=len(manifest["records"]), scientific_execution_authorized=False))
    phases = {}
    for number, side in enumerate(("before", "after"), start=1):
        checked = evidence.json(session / f"lookup_checks/checked_{number:02d}.json")
        observed = verification[side]
        require(checked.get("status") == "runtime_and_lookup_checked"
            and checked.get("scientific_execution_authorized") is False
            and same(checked["runtime"], expected) and same(observed["runtime"], checked)
            and observed["original_inputs"] == run["inputs"], "Runtime phase or input binding differs")
        require(type(verification[side + "_check_wall_s"]) in (int, float)
            and verification[side + "_check_wall_s"] > 0, "Missing runtime-check duration")
        if side == "before":
            require(observed["prepared_inputs"] is None, "Preflight reused prepared inputs")
        else:
            prepared = observed["prepared_inputs"]["datasets"]
            require(len(prepared) == 1 and prepared[0]["native_order"] == run["native_order"],
                    "Postflight native enumeration differs")
            copies = prepared[0]["inputs_in_native_order"]
            source_by_name = {Path(r["path"]).name: (r["bytes"], r["sha256"]) for r in run["inputs"]}
            require(len(copies) == len(source_by_name) and
                {Path(r["path"]).name: (r["bytes"], r["sha256"]) for r in copies} == source_by_name,
                "Postflight prepared input bytes differ")
            for ref in copies:
                require(Path(ref["path"]).parent == Path(run["output_root"]) / "input",
                        "Postflight input copy path differs")
                check(ref)
                evidence.bind(ref["path"])
        comparisons = {}
        for name in ("orthohmm", "orthofinder"):
            expected_ref = lookup["interpreters"][name]["reports"][-1]
            observed_ref = checked["lookup"][name]["report"]
            require(observed_ref["path"] == str(session / f"lookup_checks/check_{number:02d}/{name}.json"),
                    "Runtime phase lookup report path differs")
            prior, current = read(expected_ref), read(observed_ref)
            evidence.bind(expected_ref["path"])
            evidence.bind(observed_ref["path"])
            verdict = compare_lookup(prior, current)
            require(same({k: v for k, v in checked["lookup"][name].items() if k != "report"}, verdict),
                    "Retained lookup verdict does not reproduce")
            comparisons[name] = dict(verdict, report=observed_ref)
        phases[side] = dict(lookup=comparisons, inventories=expected)
    started = time.monotonic()
    current = tree_checker(binding["runtime_specs"])
    require(same(current, expected), "Fresh terminal runtime inventory does not match")
    for ref in [plan["runtime_lookup"], lookup["binding"], plan["baseline"], *manifests]:
        check(ref)
        evidence.bind(ref["path"])
    return dict(status="runtime_brackets_and_lookup_replayed", phases=phases,
        fresh_terminal_inventory=current, fresh_terminal_check_wall_s=time.monotonic() - started,
        continuous_runtime_integrity_established=False)


def environment_review(request_ref, request, plan_ref, run, result, directory, done, evidence):
    policy_ref = request["policy"]
    policy = read(policy_ref)
    require(shared_environment(policy) and native_pressure_role(policy) == "diagnostic_only"
        and policy["plan_sha256"] == plan_ref["sha256"], "Wrong shared-host environment policy")
    preflight_path = directory / "environment_preflight.json"
    preflight = evidence.json(preflight_path)
    require(preflight.get("schema") == "threadripper_environment_preflight_v1"
        and preflight.get("decision") == "passed" and preflight.get("job_id") == request["job_id"]
        and preflight.get("index") == run["index"] and preflight.get("environment_policy") == policy_ref
        and preflight.get("plan_sha256") == plan_ref["sha256"] and preflight.get("execution_scope") == SCOPE
        and preflight.get("background_competition_recorded") is True,
        "Environment preflight does not bind this shared attempt")
    retained_ref = result["environment_review"]
    require(retained_ref["path"] == str(directory / "process_stream_review.json"), "Wrong retained environment review path")
    retained = read(retained_ref)
    require(retained.get("job_id") == request["job_id"] and retained.get("index") == run["index"]
        and retained.get("execution_scope") == SCOPE and retained.get("uncontended_timing") is False
        and retained.get("background_cpu_used_for_eligibility") is False
        and retained.get("pressure_thresholds_used_for_eligibility") is False
        and retained.get("source") == record(ROOT / "benchmark_tools/review_threadripper_process_stream.py"),
        "Retained environmental scope differs")
    process_policy = read(policy["process_policy"])
    references = [policy_ref, policy["process_policy"], retained_ref,
        *policy["configuration_files"], *policy["evidence"], *preflight["evidence"], *retained["evidence"]]
    for ref in references:
        check(ref)
        evidence.bind(ref["path"])
    ready = evidence.json(directory / "ready.json")
    scope = scoped_path(ready["cgroup"], request["job_id"])
    job_scope = next(path for path in scope.parents if path.name == f"job_{request['job_id']}")
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
            "Launch observer snapshot differs from raw stream")
    launch_processes = review_process_pair(process_policy, *launch["process_snapshots"],
        boot_id=preflight["boot_id"], job_scope=str(job_scope), observer_pid=observer)
    require(same(launch["process_review"], launch_processes) and launch_processes["process_policy_matched"]
        and launch["background_cpu_used_for_eligibility"] is False, "Launch process attribution differs")
    capacity = available_memory(launch["raw_meminfo"])
    require(capacity == launch["available_memory_bytes"] and capacity >= MEMORY
        and policy["minimum_available_memory_bytes"] == MEMORY, "Unsafe or incorrectly recorded launch capacity")
    budget = evidence.json(directory / "release_budget.json")
    require(budget["status"] == "release_budget_check_passed" and budget["returncode"] == 0,
            "Release budget observation failed")
    replayed_budget = remaining_budget(budget["stdout"], request["job_id"], command=str(SCRIPT), cwd=str(ROOT),
        query_elapsed_s=budget["query_elapsed_s"], allocation_mode="shared")
    require(same(replayed_budget, budget["budget"])
        and replayed_budget["allocation"]["fields"].get("Comment") == request_ref["sha256"],
        "Release scheduler identity or remaining budget differs")
    for key, value in processes.items():
        if key not in {"schema", "limitations"}:
            require(same(retained.get(key), value), "Independent process replay differs: " + key)
    points = sorted(directory.glob("point_*.json"))
    require(points == [directory / f"point_{i:06d}.json" for i in range(len(points))], "Pressure point inventory differs")
    pressure = pressure_evaluate((evidence.json(path) for path in points), boot_id=preflight["boot_id"],
        job_scope=str(job_scope), launch_ns=done["started_ns"], end_ns=done["finished_ns"],
        limits=policy["maximum_pressure_percent"], maximum_period_s=policy["maximum_pressure_sample_period_s"],
        pressure_role="diagnostic_only")
    require(sorted(directory.glob("point_*.json")) == points, "Pressure inventory changed during review")
    passed = processes["sampled_process_policy_satisfied"] and pressure["sampled_pressure_evidence_satisfied"]
    require(same(retained["pressure_review"], pressure)
        and retained["sampled_environment_policy_satisfied"] is passed
        and result["sampled_environment_evidence_valid"] is passed, "Independent environmental verdict differs")
    return dict(status="shared_environment_replayed", sampled_environment_evidence_valid=passed,
        processes=processes, pressure=pressure, preflight=record(preflight_path), retained=retained_ref,
        launch_capacity_bytes=capacity, launch_processes=launch_processes, release_budget=replayed_budget,
        uncontended_timing=False, background_cpu_used_for_eligibility=False,
        pressure_thresholds_used_for_eligibility=False)


def classify(terminal, result, verification, replayed):
    fields = scheduler_fields(terminal)
    state, exit_code = fields.get("JobState", fields.get("State")), fields["ExitCode"]
    outcome = replayed["native_outcome"]
    if outcome == "exited_zero":
        require(state == "COMPLETED" and exit_code == "0:0"
            and verification["status"] == "command_exited_zero"
            and result["status"] == "measurement_returned_pending_independent_review",
            "Successful native outcome and terminal infrastructure disagree")
        return "native_success"
    require(outcome in {"exited_nonzero", "timed_out"} and state == "FAILED" and exit_code == "1:0"
        and verification["status"] == "command_failed"
        and result["status"] == "factorial_attempt_failed_retained",
        "Unclassified infrastructure or native failure requires separate review")
    return "native_timeout_retained" if outcome == "timed_out" else "native_failure_retained"


def resource_review(replayed, done, job):
    derived = endpoints(replayed, done, job)
    require(derived["primary_scopes"] == SCOPES and "violation" not in replayed["affinity_observation_statuses"],
            "Wrong resource scopes or observed affinity violation")
    return dict(derived, shared_host_observation=True, uncontended_timing=False,
        contention_distortion="unknown_potentially_tool_dependent", timing_disclosure=DISCLOSURE,
        affinity_observation_statuses=replayed["affinity_observation_statuses"],
        narrow_flagged_intervals=replayed["narrow_flagged_intervals"],
        whole_job_cpu_seconds=None, whole_job_peak_memory_bytes=None)


def review(request_ref, destination):
    request = read(request_ref)
    plan_ref = request["plan"]
    plan = read(plan_ref)
    runs = validate_plan(plan)
    validate_request(request, plan_ref, request["job_id"])
    terminal = verify_terminal(request["job_id"])
    fields = scheduler_fields(terminal)
    if terminal["source"] == "live_controller":
        require(fields.get("Comment") == request_ref["sha256"], "Terminal request comment differs")
    run = runs[request["index"]]
    baseline = read(plan["baseline"])
    root = Path(run["output_root"])
    session = Path(plan["panel_root"]) / "sessions" / f"run_{run['index']:02d}"
    directory = root / "measurement"
    destination = Path(destination)
    require(destination.is_absolute() and destination.resolve() == destination
        and destination.is_relative_to(ROOT) and not destination.is_relative_to(root)
        and not destination.is_relative_to(session), "Require a separate direct review destination")
    destination.mkdir(parents=True, exist_ok=False)
    source = record(__file__)
    scheduler_ref = save(destination / "scheduler.json", terminal)
    evidence = Evidence()
    try:
        result = evidence.json(session / "result.json")
        verification = evidence.json(root / "verification.json")
        bind_session(result, verification, request_ref, request, plan_ref, run)
        for ref in [request_ref, plan_ref, plan["baseline"], *plan["helper_sources"], *plan["evidence"]]:
            check(ref)
            evidence.bind(ref["path"])
        runtime = runtime_review(plan, run, session, verification, evidence)
        runtime_ref = save(destination / "runtime.json", runtime)
        replayed = replay(directory, request["job_id"], native_command(plan_ref, run, baseline))
        require(same(replayed["measured"], verification["measurement"]), "Collector and runtime wrapper differ")
        for ref in replayed["evidence"]:
            check(ref)
            evidence.bind(ref["path"])
        replay_ref = save(destination / "resource_replay.json", replayed)
        done = evidence.json(directory / "done.json")
        status = classify(terminal, result, verification, replayed)
        resources = resource_review(replayed, done, request["job_id"])
        resources_ref = save(destination / "resources.json", resources)
        environment = environment_review(request_ref, request, plan_ref, run, result, directory, done, evidence)
        environment_ref = save(destination / "environment.json", environment)
        outputs = None
        if status == "native_success":
            outputs = validate_outputs(request_ref)
            require(outputs["job_id"] == request["job_id"] and outputs["index"] == run["index"]
                and outputs["plan"] == plan_ref and outputs["native_outputs_validated"] is True,
                "Output review identity differs")
            for ref in [*outputs["evidence"], *outputs["checked_files"]]:
                check(ref)
                evidence.bind(ref["path"])
        outputs_ref = save(destination / "outputs_or_failure.json", outputs if outputs is not None else
            dict(status=status, native_outcome=replayed["native_outcome"], native_exit_code=replayed["native_exit_code"],
                accuracy_evaluated=False, native_outputs_validated=False, automatic_retry=False))
        sampled_valid = environment["sampled_environment_evidence_valid"]
        checked = evidence.finish()
        check(source)
        return save(destination / "review.json", dict(schema="native_factorial_terminal_review_v1", status=status,
            index=run["index"], job_id=request["job_id"], dataset=run["dataset"], cell=run["cell"], repeat=run["repeat"],
            plan=plan_ref, request=request_ref, source=source, scheduler=scheduler_ref,
            scheduler_state=fields.get("JobState", fields.get("State")), scheduler_exit_code=fields["ExitCode"],
            reviews=dict(runtime=runtime_ref, resources=resources_ref, environment=environment_ref,
                         outputs_or_failure=outputs_ref), resource_replay=replay_ref, evidence=checked,
            terminal_reviewed=True, next_identity_authorized=sampled_valid,
            resources=resources["primary"], resource_scopes=SCOPES,
            shared_host_resources_reviewed=sampled_valid, primary_resources_replayed=True,
            native_outputs_validated=outputs is not None, accuracy_evaluated=False,
            execution_scope=SCOPE, uncontended_timing=False, scientific_timings_admitted=False,
            contention_distortion="unknown_potentially_tool_dependent", timing_disclosure=DISCLOSURE,
            whole_run_maximum_foreign_average_cores=environment["processes"]["maximum_observed_foreign_average_cores"],
            automatic_retry=False, whole_job_cpu_seconds=None, whole_job_peak_memory_bytes=None,
            limitations=["Independent terminal review is not benchmark scoring or publication readiness.",
                "Next identity authorization concerns only a different frozen identity; its fresh launch checks remain mandatory.",
                "No retry, fastest-repeat selection, background subtraction or isolated-performance ranking.",
                "CPU includes wrapper bracket work; memory is native-step lifetime peak, not algorithm-only RSS.",
                "Periodic monitoring and runtime brackets do not establish continuous containment or integrity."]))
    except BaseException as error:
        save(destination / "failure.json", dict(status="terminal_factorial_review_failed", index=run["index"],
            job_id=request["job_id"], plan=plan_ref, request=request_ref, scheduler=scheduler_ref,
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
    print(json.dumps(review(ref, args.output_directory.resolve()), sort_keys=True))
