"""Explicit final-identity review using unchanged scientific/accounting kernels."""

import argparse
import json
from pathlib import Path
import subprocess
import time

from benchmark_tools import run_native12_composed_cost as executor
from benchmark_tools.native12_composed_execution import current_trees
from benchmark_tools.native_factorial_allocated_execution import terminal_accounting, native_command, SCOPE, MEMORY
from benchmark_tools.run_allocated_native_factorial_cost import native_placement
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.replay_allocated_threadripper_scaling import replay
from benchmark_tools.review_native_factorial_attempt import (
    save, bind_session, scheduler_fields, classify, resource_review, DISCLOSURE)
from benchmark_tools.inspect_native_python_lookup import compare_lookup
from benchmark_tools.review_threadripper_process_policy import shared_environment, review as review_process_pair
from benchmark_tools.review_threadripper_process_stream import evaluate as process_evaluate
from benchmark_tools.review_threadripper_pressure_stream import evaluate as pressure_evaluate, native_pressure_role
from benchmark_tools.run_native_factorial_cost import ROOT, read, available_memory
from benchmark_tools.slurm_resource_snapshot import scoped_path
from benchmark_tools.validate_native_factorial_outputs import Evidence, require, validate_semantics
from benchmark_tools.verify_lineage_native_provenance import same
from benchmark_tools.verify_threadripper_controller import validate as controller, remaining_budget
from benchmark_tools.derive_threadripper_resources import SCOPES


SCHEMA = "native12_composed_terminal_review_v1"
DESTINATION = ROOT / "benchmarks/work/native12_composed_terminal_review_20261008_v1"


def verify_terminal(job):
    require(type(job) is int and job > 24034, "Wrong final native job identity")
    query = ["scontrol", "show", "job", str(job), "--oneliner"]
    result = subprocess.run(query, capture_output=True, text=True, timeout=5)
    observation = dict(command=query, returncode=result.returncode, stdout=result.stdout, stderr=result.stderr)
    if result.returncode == 0:
        verified = controller(result.stdout, job, "terminal", command=str(executor.SCRIPT), cwd=str(ROOT),
            time_limit="1-02:00:00", allocation_mode="shared")
        require(verified["fields"].get("JobName") == executor.JOB_NAME, "Final native job name differs")
        return dict(source="live_controller", observation=observation, verified=verified)
    require("Invalid job id specified" in result.stderr, "Controller failed without established expiry")
    query = ["sacct", "-X", "-j", str(job), "-n", "-P",
        "--format=JobIDRaw,State,ExitCode,AllocCPUS,ReqMem,NodeList,Partition,TimelimitRaw,Submit,Start,End,JobName"]
    result = subprocess.run(query, capture_output=True, text=True, check=True, timeout=5)
    return dict(source="fresh_accounting_after_controller_expiry", controller_observation=observation,
        observation=dict(command=query, stdout=result.stdout, stderr=result.stderr),
        verified=terminal_accounting(result.stdout, job),
        limitation="Accounting lacks historical command/comment; exact held request/source evidence is also required.")


def session_gate(result, verification, request_ref, request, run):
    bind_session(result, verification, request_ref, request, request["plan"], run)
    require(result.get("schema") == executor.SESSION_SCHEMA
        and result.get("source") == record(executor.__file__)
        and result.get("amendment") == request["amendment"],
        "Final native session/controller source differs")
    history = result.get("history", {})
    require(history.get("schema") == executor.HISTORY_SCHEMA
        and history.get("original_request") == request["original_request"]
        and history.get("composed_review") == request["composed_review"]
        and history.get("historical_prefix") == request["history"][:-1]
        and history.get("next_unrun_index") == 12
        and history.get("original_review_translated") is False
        and history.get("next_identity_authorized") is False,
        "Retained native12 session history basis differs")


def runtime_review(plan, run, session, verification, request, evidence):
    lookup = read(plan["runtime_lookup"])
    lookup_binding = read(lookup["binding"])
    specs = lookup_binding["runtime_specs"]
    require(lookup["baseline"] == plan["baseline"]
        and lookup["source"] == record(ROOT / "benchmark_tools/inspect_native_python_lookup.py"),
        "Final native lookup baseline/source differs")
    started = time.monotonic()
    fresh = current_trees(specs, request["runtime_basis"])
    terminal_check_wall_s = time.monotonic() - started
    phases = {}
    for number, side in enumerate(("before", "after"), start=1):
        checked = evidence.json(session / f"lookup_checks/checked_{number:02d}.json")
        tree_check = evidence.json(session / f"lookup_checks/tree_check_{number:02d}.json")
        observed = verification[side]
        require(checked.get("status") == "runtime_and_lookup_checked"
            and checked.get("scientific_execution_authorized") is False
            and checked.get("runtime") == fresh and observed.get("runtime") == checked
            and observed.get("original_inputs") == run["inputs"]
            and type(verification.get(side + "_check_wall_s")) in (int, float)
            and verification[side + "_check_wall_s"] > 0
            and tree_check.get("schema") == "native12_current_runtime_check_v1"
            and tree_check.get("source") == record(ROOT / "benchmark_tools/native12_composed_execution.py")
            and tree_check.get("runtime_basis") == request["runtime_basis"]
            and tree_check.get("inventories") == fresh
            and tree_check.get("original_os_inventory_equality") is False
            and tree_check.get("current_identity_revalidated") is True
            and tree_check.get("continuous_runtime_integrity_established") is False
            and type(tree_check.get("check_wall_s")) in (int, float) and tree_check["check_wall_s"] > 0,
            "Final native runtime phase is incomplete or does not reproduce")
        if side == "before":
            require(observed.get("prepared_inputs") is None, "Preflight reused prepared inputs")
        else:
            copied = observed["prepared_inputs"]["datasets"]
            require(len(copied) == 1 and copied[0]["native_order"] == run["native_order"],
                "Final native input order differs")
            rows = copied[0]["inputs_in_native_order"]
            original = {Path(ref["path"]).name: (ref["bytes"], ref["sha256"]) for ref in run["inputs"]}
            require(len(rows) == len(original)
                and {Path(ref["path"]).name: (ref["bytes"], ref["sha256"]) for ref in rows} == original,
                "Final native prepared input bytes differ")
            for ref in rows:
                require(Path(ref["path"]).parent == Path(run["output_root"]) / "input",
                    "Final native copied input path differs")
                check(ref)
                evidence.bind(ref["path"])
        comparisons = {}
        for name in ("orthohmm", "orthofinder"):
            prior_ref = lookup["interpreters"][name]["reports"][-1]
            observed_ref = checked["lookup"][name]["report"]
            require(observed_ref["path"] == str(session / f"lookup_checks/check_{number:02d}/{name}.json"),
                "Final native lookup report path differs")
            verdict = compare_lookup(read(prior_ref), read(observed_ref))
            require({key: value for key, value in checked["lookup"][name].items() if key != "report"} == verdict,
                "Final native import verdict does not reproduce")
            for ref in (prior_ref, observed_ref):
                evidence.bind(ref["path"])
            comparisons[name] = dict(verdict, report=observed_ref)
        phases[side] = dict(inventories=fresh, lookup=comparisons)
    for ref in [plan["runtime_lookup"], lookup["binding"], plan["baseline"], request["runtime_basis"],
                *[row["manifest"] for row in fresh]]:
        check(ref)
        evidence.bind(ref["path"])
    return dict(schema="native12_composed_runtime_review_v1", status="fresh_runtime_brackets_and_lookup_replayed",
        phases=phases, fresh_terminal_inventory=fresh, runtime_basis=request["runtime_basis"],
        fresh_terminal_check_wall_s=terminal_check_wall_s,
        current_original_os_inventory_equality=False, prospective_current_inventory_equality=True,
        continuous_runtime_integrity_established=False)


def environment_review(request_ref, request, run, result, directory, done, evidence):
    policy_ref = request["policy"]
    policy = read(policy_ref)
    require(shared_environment(policy) and native_pressure_role(policy) == "diagnostic_only"
        and policy["plan_sha256"] == request["plan"]["sha256"], "Wrong final native shared-host policy")
    preflight = evidence.json(directory / "environment_preflight.json")
    require(preflight.get("schema") == "threadripper_environment_preflight_v1"
        and preflight.get("decision") == "passed" and preflight.get("job_id") == request["job_id"]
        and preflight.get("index") == 12 and preflight.get("environment_policy") == policy_ref
        and preflight.get("plan_sha256") == request["plan"]["sha256"]
        and preflight.get("amendment") == request["amendment"] and preflight.get("execution_scope") == SCOPE
        and preflight.get("controller_schema") == executor.REQUEST_SCHEMA
        and preflight.get("background_competition_recorded") is True, "Final native launch policy differs")
    retained_ref = result["environment_review"]
    require(retained_ref["path"] == str(directory / "process_stream_review.json"), "Wrong final native stream review")
    retained = read(retained_ref)
    require(retained.get("job_id") == request["job_id"] and retained.get("index") == 12
        and retained.get("execution_scope") == SCOPE and retained.get("uncontended_timing") is False
        and retained.get("background_cpu_used_for_eligibility") is False
        and retained.get("pressure_thresholds_used_for_eligibility") is False
        and retained.get("source") == record(ROOT / "benchmark_tools/review_threadripper_process_stream.py"),
        "Final native retained process/pressure scope differs")
    for ref in [policy_ref, policy["process_policy"], retained_ref, *policy["configuration_files"],
                *policy["evidence"], *preflight["evidence"], *retained["evidence"]]:
        check(ref)
        evidence.bind(ref["path"])
    ready = evidence.json(directory / "ready.json")
    scope = scoped_path(ready["cgroup"], request["job_id"])
    job_scope = next(parent for parent in scope.parents if parent.name == f"job_{request['job_id']}")
    path = evidence.bind(directory / "host_processes.jsonl")
    with path.open() as stream:
        first = json.loads(stream.readline())
        stream.seek(0)
        processes = process_evaluate(stream, read(policy["process_policy"]), boot_id=preflight["boot_id"],
            job_scope=str(job_scope), observer_pid=first["observer_pid"], launch=done["started_ns"] / 1e9,
            end=done["finished_ns"] / 1e9, maximum_foreign_average_cores=policy["maximum_foreign_average_cores"],
            maximum_sample_period_s=policy["maximum_sample_period_s"])
    launch = evidence.json(directory / "launch_environment_observation.json")
    require(len(launch["process_snapshots"]) == 2 and launch["process_snapshots"][0] == first["snapshot"],
        "Final native launch snapshot differs")
    pair = review_process_pair(read(policy["process_policy"]), *launch["process_snapshots"],
        boot_id=preflight["boot_id"], job_scope=str(job_scope), observer_pid=first["observer_pid"])
    capacity = available_memory(launch["raw_meminfo"])
    require(pair == launch["process_review"] and pair["process_policy_matched"]
        and launch["background_cpu_used_for_eligibility"] is False
        and capacity == launch["available_memory_bytes"] and capacity >= MEMORY
        and policy["minimum_available_memory_bytes"] == MEMORY, "Final native capacity/attribution differs")
    budget = evidence.json(directory / "release_budget.json")
    require(budget.get("status") == "release_budget_check_passed" and budget.get("returncode") == 0,
        "Final native release budget failed")
    expected_budget = remaining_budget(budget["stdout"], request["job_id"], command=str(executor.SCRIPT), cwd=str(ROOT),
        query_elapsed_s=budget["query_elapsed_s"], allocation_mode="shared")
    fields = expected_budget["allocation"]["fields"]
    require(same(expected_budget, budget["budget"]) and fields.get("Comment") == request_ref["sha256"]
        and fields.get("JobName") == executor.JOB_NAME, "Final native release identity/budget differs")
    for key, value in processes.items():
        if key not in {"schema", "limitations"}:
            require(same(retained.get(key), value), "Final native process replay differs: " + key)
    paths = sorted(directory.glob("point_*.json"))
    require(paths == [directory / f"point_{index:06d}.json" for index in range(len(paths))],
        "Final native pressure point inventory differs")
    pressure = pressure_evaluate((evidence.json(path) for path in paths), boot_id=preflight["boot_id"],
        job_scope=str(job_scope), launch_ns=done["started_ns"], end_ns=done["finished_ns"],
        limits=policy["maximum_pressure_percent"], maximum_period_s=policy["maximum_pressure_sample_period_s"],
        pressure_role="diagnostic_only")
    require(sorted(directory.glob("point_*.json")) == paths, "Final native pressure inventory changed")
    passed = processes["sampled_process_policy_satisfied"] and pressure["sampled_pressure_evidence_satisfied"]
    require(same(retained["pressure_review"], pressure)
        and retained["sampled_environment_policy_satisfied"] is passed
        and result["sampled_environment_evidence_valid"] is passed, "Final native environment verdict differs")
    return dict(status="shared_environment_replayed", sampled_environment_evidence_valid=passed,
        amendment=request["amendment"], processes=processes, pressure=pressure, retained=retained_ref,
        preflight=record(directory / "environment_preflight.json"), launch_capacity_bytes=capacity,
        launch_processes=pair, release_budget=expected_budget, uncontended_timing=False,
        background_cpu_used_for_eligibility=False, pressure_thresholds_used_for_eligibility=False)


def scientific_output(request_ref, request, execution, plan, run, baseline, terminal):
    fields = scheduler_fields(terminal)
    require(fields.get("JobState", fields.get("State")) == "COMPLETED" and fields.get("ExitCode") == "0:0",
        "Semantic review requires actual successful final native completion")
    root = Path(run["output_root"])
    native_ref, ready_ref = record(root / "native_execution.json"), record(root / "measurement/ready.json")
    native, ready = read(native_ref), read(ready_ref)
    allowed = native_placement(ready, native["placement"], request["job_id"], native["parent_pid"])
    require(native.get("schema") == "allocated_native_factorial_execution_v1"
        and native.get("status") == "native_factorial_completed_pending_output_review"
        and native.get("plan") == request["plan"] and native.get("amendment") == request["amendment"]
        and native.get("index") == 12 and native.get("cell") == run["cell"]
        and native.get("allocated_ready") == ready_ref and native.get("native_cpu_ids") == allowed
        and native.get("native_order") == run["native_order"] and native.get("automatic_retry") is False
        and native.get("source") == record(ROOT / "benchmark_tools/run_allocated_native_factorial_cost.py"),
        "Final native unchanged scientific producer/placement differs")
    refs = [request_ref, request["amendment"], request["plan"], plan["baseline"], native_ref, ready_ref,
        *execution["new_sources"], *request["new_sources"], *plan["helper_sources"], *run["inputs"]]
    for ref in refs:
        check(ref)
    context = dict(run, input_directory=str(root / "input"), cpu=32, threads_per_worker=4,
        command=native_command(request["amendment"], run, baseline, metrics=True), cwd=baseline["core_root"],
        aligner=baseline["tool_entrypoints"]["mafft"]["absolute_path"],
        tree_builder=baseline["tool_entrypoints"]["FastTree"]["absolute_path"])
    result = validate_semantics(context)
    require(result.get("native_outputs_validated") is True
        and result.get("accuracy_evaluated") is False,
        "Original semantic kernel did not validate final native output")
    preparation_ref = record(root / "preparation.json")
    preparation = read(preparation_ref)
    require(preparation.get("status") == "fresh_factorial_inputs_prepared"
        and preparation.get("gene_ownership_sha256") == result["gene_ownership_sha256"]
        and preparation.get("per_species_counts") == result["per_species_counts"]
        and preparation.get("genes") == run["genes"] and native["factors"] == result["factors"],
        "Final native preparation/factors do not reproduce")
    refs.append(preparation_ref)
    for ref in refs:
        check(ref)
    return dict(result, schema="native12_composed_output_review_v1", semantic_validator_source=result["source"],
        source=record(__file__), index=12, job_id=request["job_id"], plan=request["plan"],
        amendment=request["amendment"], request=request_ref, native_cpu_ids=allowed, allocated_ready=ready_ref,
        scheduler=terminal, execution_scope=SCOPE, evidence=refs, terminal_reviewed=False,
        terminal_scheduler_confirmed=True, uncontended_timing=False)


def review(request_ref, destination, expected_source_sha256):
    destination = Path(destination)
    require(destination == DESTINATION and destination.is_absolute() and destination.resolve() == destination
        and not destination.exists() and not destination.is_symlink(), "Require one fresh final native review destination")
    source = record(__file__)
    require(source["sha256"] == expected_source_sha256, "Prospective final native reviewer changed")
    request = read(request_ref)
    request, context, _ = executor.execution_binding(request_ref, request["job_id"])
    execution, plan = context[1:3]
    executor.held_gate(request["held_scheduler"]["stdout"], request["job_id"])
    terminal = verify_terminal(request["job_id"])
    fields = scheduler_fields(terminal)
    if terminal["source"] == "live_controller":
        require(fields.get("Comment") == request_ref["sha256"], "Final native terminal request digest differs")
    run, baseline = plan["runs"][12], read(plan["baseline"])
    root = Path(run["output_root"])
    session = Path(plan["panel_root"]) / "sessions/run_12"
    directory = root / "measurement"
    destination.mkdir(parents=True, exist_ok=False)
    scheduler_ref = save(destination / "scheduler.json", terminal)
    evidence = Evidence()
    try:
        result, verification = evidence.json(session / "result.json"), evidence.json(root / "verification.json")
        session_gate(result, verification, request_ref, request, run)
        for ref in [source, request_ref, request["amendment"], request["plan"], plan["baseline"],
                    *execution["new_sources"], *request["new_sources"], *plan["helper_sources"], *plan["evidence"],
                    *request["history_basis"]["evidence"]]:
            check(ref)
            evidence.bind(ref["path"])
        runtime_ref = save(destination / "runtime.json", runtime_review(plan, run, session, verification, request, evidence))
        replayed = replay(directory, request["job_id"], native_command(request["amendment"], run, baseline))
        require(same(replayed["measured"], verification["measurement"]), "Final native collector/runtime wrapper differs")
        for ref in replayed["evidence"]:
            check(ref)
            evidence.bind(ref["path"])
        replay_ref = save(destination / "resource_replay_summary.json", dict(
            schema="native12_composed_resource_replay_summary_v1", source=record(ROOT / "benchmark_tools/replay_allocated_threadripper_scaling.py"),
            native_outcome=replayed["native_outcome"], native_exit_code=replayed["native_exit_code"],
            measured_matches_retained_wrapper=True, full_replay_executed=True, evidence=replayed["evidence"]))
        done = evidence.json(directory / "done.json")
        status = classify(terminal, result, verification, replayed)
        resources = resource_review(replayed, done, request["job_id"])
        resources_ref = save(destination / "resources.json", resources)
        environment = environment_review(request_ref, request, run, result, directory, done, evidence)
        environment_ref = save(destination / "environment.json", environment)
        outputs = scientific_output(request_ref, request, execution, plan, run, baseline, terminal) if status == "native_success" else None
        if outputs is not None:
            require(outputs.get("schema") == "native12_composed_output_review_v1"
                and outputs.get("index") == 12 and outputs.get("job_id") == request["job_id"]
                and outputs.get("plan") == request["plan"] and outputs.get("amendment") == request["amendment"]
                and outputs.get("request") == request_ref and outputs.get("native_outputs_validated") is True,
                "Final native semantic review identity differs")
            for ref in [*outputs["evidence"], *outputs["checked_files"]]:
                check(ref)
                evidence.bind(ref["path"])
        outputs_ref = save(destination / "outputs_or_failure.json", outputs if outputs is not None else
            dict(status=status, native_outcome=replayed["native_outcome"], native_exit_code=replayed["native_exit_code"],
                native_outputs_validated=False, accuracy_evaluated=False, automatic_retry=False))
        checked = evidence.finish()
        check(source)
        return save(destination / "review.json", dict(schema=SCHEMA, source=source, status=status, index=12,
            job_id=request["job_id"], cell=run["cell"], dataset=run["dataset"], repeat=run["repeat"],
            request=request_ref, plan=request["plan"], amendment=request["amendment"], scheduler=scheduler_ref,
            scheduler_state=fields.get("JobState", fields.get("State")), scheduler_exit_code=fields["ExitCode"],
            reviews=dict(runtime=runtime_ref, resources=resources_ref, environment=environment_ref, outputs_or_failure=outputs_ref),
            resource_replay=replay_ref, resources=resources["primary"], resource_scopes=SCOPES,
            native_cpu_ids=replayed["native_cpu_ids"], allocated_placement=replayed["allocated_placement"],
            terminal_reviewed=True, primary_resources_replayed=True,
            shared_host_resources_reviewed=environment["sampled_environment_evidence_valid"],
            native_outputs_validated=outputs is not None, next_identity_authorized=False,
            original_review_translated=False, current_original_os_inventory_equality=False,
            prospective_current_inventory_equality=True, continuous_runtime_integrity_established=False,
            execution_scope=SCOPE, uncontended_timing=False, scientific_timings_admitted=False,
            accuracy_evaluated=False, automatic_retry=False, publication_ready=False,
            timing_disclosure=DISCLOSURE, contention_distortion="unknown_potentially_tool_dependent",
            historical_review_failures_retained=[23986, 24033], evidence=checked,
            limitations=["Explicit new request/session/runtime lineage; not ordinary-review substitution.",
                "Full original accounting and semantic kernels replayed; expanded host JSON not duplicated.",
                "No scoring, inference retry, continuous integrity or publication-readiness claim."]))
    except BaseException as error:
        save(destination / "failure.json", dict(schema="native12_composed_review_failure_v1", source=source,
            status="terminal_review_failed", request=request_ref, job_id=request["job_id"], index=12,
            scheduler=scheduler_ref, error_type=type(error).__name__, error=str(error), terminal_reviewed=False,
            next_identity_authorized=False, automatic_retry=False, accuracy_evaluated=False))
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request", type=Path, required=True)
    parser.add_argument("--request-sha256", required=True)
    parser.add_argument("--source-sha256", required=True)
    args = parser.parse_args()
    ref = record(args.request)
    require(ref["sha256"] == args.request_sha256, "Final native request digest differs")
    print(review(ref, DESTINATION, args.source_sha256))
