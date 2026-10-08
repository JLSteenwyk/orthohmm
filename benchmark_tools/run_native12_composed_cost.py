"""Prospective controller for the final frozen native identity, not a retry."""

import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from benchmark_tools.native12_composed_execution import (
    history_binding, require_unrun, ProspectiveRuntimeChecker, HISTORY_SCHEMA)
from benchmark_tools.native_factorial_allocated_execution import amendment, native_command, SCOPE, MEMORY
from benchmark_tools.native_factorial_allocated_placement import validate as validate_placement
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.run_native_factorial_cost import (
    CORE_COMMIT, read, available_memory, original_inputs, prepared_inputs, prepare_inputs)
from benchmark_tools.verify_threadripper_controller import validate as controller, ReleaseBudgetGuard
from benchmark_tools.validate_native_factorial_outputs import require


SCRIPT = ROOT / "benchmark_tools/run_native12_composed_cost.sh"
REQUEST = ROOT / "benchmarks/work/native_factorial_launch_20261004/request_12_composed_v1.json"
JOB_NAME = "orthohmm_allocated_factorial"
REQUEST_SCHEMA = "native12_composed_request_v1"
SESSION_SCHEMA = "native12_composed_session_v1"
SOURCE_NAMES = ("run_native12_composed_cost.py", "run_native12_composed_cost.sh",
                "native12_composed_execution.py", "native11_composed_review_binding.py")


def sources():
    return [record(ROOT / "benchmark_tools" / name) for name in SOURCE_NAMES]


def held_gate(raw, job):
    pairs = re.findall(r"(?<!\S)([A-Za-z][^\s=]*)=([^\s]+)", raw)
    fields = dict(pairs)
    expected = dict(JobId=str(job), JobName=JOB_NAME, JobState="PENDING", Reason="JobHeldUser",
        Partition="gpu", ReqNodeList="bizon", NumCPUs="64", NumTasks="1", MinMemoryNode="128G",
        Requeue="0", Restarts="0", TimeLimit="1-02:00:00", Command=str(SCRIPT), WorkDir=str(ROOT))
    expected["CPUs/Task"] = "64"
    require(type(job) is int and job > 24034 and len(fields) == len(pairs)
        and all(fields.get(key) == value for key, value in expected.items())
        and fields.get("NumNodes") in {"1", "1-1"}
        and fields.get("UserId", "").endswith("(" + str(os.getuid()) + ")")
        and not any(key in fields for key in ("ArrayJobId", "ArrayTaskId", "HetJobId", "HetJobOffset")),
        "Held native12 ownership, command or resource envelope differs")
    return fields


def request_gate(request, job, context):
    original, execution, plan, _, review = context[:5]
    require(request.get("schema") == REQUEST_SCHEMA and type(job) is int and job > 24034
        and type(request.get("job_id")) is int and request["job_id"] == job
        and type(request.get("index")) is int and request["index"] == 12
        and request.get("cell") == "p1_c1_r1" and request.get("execution_authorized") is True
        and request.get("plan") == original["plan"]
        and request.get("amendment") == original["amendment"]
        and request.get("policy") == original["policy"]
        and request.get("scheduler_command") == str(SCRIPT) and request.get("allocation_cwd") == str(ROOT)
        and request.get("automatic_retry") is False
        and request.get("scientific_method_modified") is False
        and request.get("original_review_translated") is False
        and request.get("scientific_timings_admitted") is False
        and request.get("publication_ready") is False
        and request.get("source") == record(__file__)
        and request.get("history") == original["history"] + [request.get("composed_review")]
        and request.get("runtime_basis") == review["reviews"]["runtime"]
        and request.get("new_sources") == sources(),
        "Prospective native12 request differs from exact frozen scope/history/sources")
    history = request.get("history_basis", {})
    require(history.get("schema") == HISTORY_SCHEMA
        and history.get("status") == "prefix_revalidated_for_native12_preparation"
        and history.get("original_request") == request.get("original_request")
        and history.get("composed_review") == request.get("composed_review")
        and history.get("historical_prefix") == original["history"]
        and history.get("next_unrun_index") == 12
        and history.get("original_review_translated") is False
        and history.get("historical_failures_retained") == [23986, 24033]
        and history.get("next_identity_authorized") is False,
        "Request lacks explicit new-type reviewed-history basis")
    require(plan["runs"][12]["index"] == 12 and plan["runs"][12]["cell"] == request["cell"]
        and execution["historical_plan"] == request["plan"], "Frozen native12 recipe differs")


def prepare(job, original_ref, composed_ref, expected_source_sha256):
    source = record(__file__)
    require(source["sha256"] == expected_source_sha256, "Prospective controller source changed")
    context, history = history_binding(original_ref, composed_ref)
    original, _, plan, _, review = context[:5]
    require_unrun(plan["runs"][12], plan, REQUEST)
    raw_memory = Path("/proc/meminfo").read_text()
    capacity = available_memory(raw_memory)
    require(capacity >= MEMORY, "Unsafe available memory for native12 preparation")
    query = ["scontrol", "show", "job", str(job), "--oneliner"]
    result = subprocess.run(query, capture_output=True, text=True, check=True, timeout=5)
    held_gate(result.stdout, job)
    request = dict(schema=REQUEST_SCHEMA, execution_authorized=True, job_id=job, index=12,
        cell="p1_c1_r1", plan=original["plan"], amendment=original["amendment"], policy=original["policy"],
        original_request=original_ref, composed_review=composed_ref,
        history=original["history"] + [composed_ref], history_basis=history,
        runtime_basis=review["reviews"]["runtime"], new_sources=sources(),
        scheduler_command=str(SCRIPT), allocation_cwd=str(ROOT), automatic_retry=False,
        scientific_method_modified=False, original_review_translated=False,
        scientific_timings_admitted=False, publication_ready=False,
        source=source, source_commit=subprocess.check_output(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
        prepared_unix_ns=time.time_ns(),
        held_scheduler=dict(command=query, stdout=result.stdout, stderr=result.stderr),
        capacity_precheck=dict(raw_meminfo=raw_memory, available_memory_bytes=capacity,
            capacity_guaranteed_through_run=False, background_cpu_used_for_eligibility=False))
    request_gate(request, job, context)
    for ref in [original_ref, composed_ref, request["runtime_basis"], *request["new_sources"], *history["evidence"]]:
        check(ref)
    require_unrun(plan["runs"][12], plan, REQUEST)
    save(REQUEST, request)
    return record(REQUEST)


def execution_binding(request_ref, job):
    require(request_ref["path"] == str(REQUEST), "Wrong native12 request destination")
    request = read(request_ref)
    context, history = history_binding(request["original_request"], request["composed_review"])
    request_gate(request, job, context)
    for ref in [request_ref, request["runtime_basis"], *request["new_sources"], *history["evidence"],
                *request["history_basis"]["evidence"]]:
        check(ref)
    return request, context, history


class ReleaseGuard:
    def __init__(self, request_ref, request, run, baseline):
        self.request_ref, self.request = request_ref, request
        self.run, self.baseline = run, baseline

    def __call__(self, directory):
        from benchmark_tools.probe_host_counters import snapshot
        from benchmark_tools.observe_threadripper_process_identity import enriched_snapshot
        from benchmark_tools.review_threadripper_process_policy import review
        from benchmark_tools.slurm_resource_snapshot import scoped_path

        request = read(self.request_ref)
        require(request == self.request, "Prepared native12 request changed before release")
        amendment(request["amendment"])
        for ref in request["new_sources"]:
            check(ref)
        prepared_inputs(self.run, self.baseline)
        ready = read(record(directory / "ready.json"))
        allowed = validate_placement(ready["allocated_placement"], request["job_id"])
        require(ready["placement"]["affinity"] == allowed, "Ready native12 worker affinity differs")
        scope = scoped_path(ready["cgroup"], request["job_id"])
        job_scope = next(parent for parent in scope.parents if parent.name == f"job_{request['job_id']}")
        with (directory / "host_processes.jsonl").open() as stream:
            first = json.loads(stream.readline())
        second = enriched_snapshot()
        boot = Path("/proc/sys/kernel/random/boot_id").read_text().strip()
        policy = read(request["policy"])
        verdict = review(read(policy["process_policy"]), first["snapshot"], second,
            boot_id=boot, job_scope=str(job_scope), observer_pid=first["observer_pid"])
        raw_memory = Path("/proc/meminfo").read_text()
        capacity = available_memory(raw_memory)
        save(directory / "launch_environment_observation.json", dict(host_counters=snapshot(),
            raw_meminfo=raw_memory, available_memory_bytes=capacity,
            process_snapshots=[first["snapshot"], second], process_review=verdict,
            background_cpu_used_for_eligibility=False, capacity_guaranteed_through_run=False))
        require(capacity >= MEMORY and verdict["process_policy_matched"],
            "Unsafe launch memory or invalid process attribution")
        budget = ReleaseBudgetGuard(request["job_id"], command=str(SCRIPT), cwd=str(ROOT),
            allocation_mode="shared")(directory)
        fields = budget["allocation"]["fields"]
        require(fields.get("Comment") == self.request_ref["sha256"] and fields.get("JobName") == JOB_NAME,
            "Release scheduler request/name differs")
        save(directory / "environment_preflight.json", dict(schema="threadripper_environment_preflight_v1",
            decision="passed", job_id=request["job_id"], index=12, plan_sha256=request["plan"]["sha256"],
            amendment=request["amendment"], environment_policy=request["policy"], execution_scope=SCOPE,
            boot_id=boot, background_competition_recorded=True, uncontended_timing=False,
            controller_schema=REQUEST_SCHEMA,
            evidence=[record(directory / "launch_environment_observation.json"), record(directory / "release_budget.json")]))
        return budget


def execute(request_ref):
    from benchmark_tools.isolated_numba_cache import fresh_cache
    from benchmark_tools.isolated_native_tmp import fresh_tmp
    from benchmark_tools.measure_native_scaling_run import run_checked
    from benchmark_tools.measure_threadripper_run import assert_environment
    from benchmark_tools.measure_allocated_threadripper_scaling import measure
    from benchmark_tools.review_threadripper_process_stream import audit, shared_environment
    from benchmark_tools.run_simulation_methods import execution_environment, verify_environment

    job = int(os.environ["SLURM_JOB_ID"])
    request, context, history = execution_binding(request_ref, job)
    execution, plan = context[1:3]
    require(Path.cwd() == ROOT and os.uname().nodename == "bizon" and sys.dont_write_bytecode
        and bool(sys.pycache_prefix) and not Path(sys.pycache_prefix).exists()
        and not Path(sys.pycache_prefix).is_symlink()
        and os.environ.get("PYTHONHASHSEED") == "0" and os.environ.get("PYTHONNOUSERSITE") == "1",
        "Wrong native12 controller bootstrap")
    lookup = read(plan["runtime_lookup"])
    lookup_binding = read(lookup["binding"])
    baseline = read(plan["baseline"])
    require(lookup["baseline"] == plan["baseline"] and baseline["core_commit"] == CORE_COMMIT
        and os.path.abspath(sys.executable) == lookup_binding["controller_python"]["path"],
        "Native12 runtime/private controller differs")
    check(lookup_binding["controller_python"])
    policy = read(request["policy"])
    require(shared_environment(policy) and policy["plan_sha256"] == request["plan"]["sha256"]
        and policy["minimum_available_memory_bytes"] == MEMORY, "Wrong shared-host native12 policy")
    query = ["scontrol", "show", "job", str(job), "--oneliner"]
    observation = subprocess.run(query, check=True, capture_output=True, text=True, timeout=5)
    fields = controller(observation.stdout, job, "running", command=str(SCRIPT), cwd=str(ROOT),
        time_limit="1-02:00:00", allocation_mode="shared")["fields"]
    require(fields.get("Comment") == request_ref["sha256"] and fields.get("JobName") == JOB_NAME
        and fields.get("UserId", "").endswith("(" + str(os.getuid()) + ")"),
        "Running allocation does not bind native12 request")
    run = plan["runs"][12]
    root = Path(run["output_root"])
    session = Path(plan["panel_root"]) / "sessions/run_12"
    require(not root.exists() and not root.is_symlink() and not session.exists() and not session.is_symlink(),
        "Native12 attempt already exists")
    session.mkdir(parents=True, exist_ok=False)
    outcome = dict(schema=SESSION_SCHEMA, status="factorial_executor_started", job_id=job,
        index=12, cell=run["cell"], request=request_ref, plan=request["plan"], amendment=request["amendment"],
        source=record(__file__), execution_scope=SCOPE, history=history,
        automatic_retry=False, scientific_timings_admitted=False, next_identity_authorized=False,
        uncontended_timing=False, contention_distortion="unknown_potentially_method_dependent")
    save(session / "started.json", outcome)
    previous_cwd = Path.cwd()
    try:
        env, _ = execution_environment(baseline)
        env.update(PYTHONDONTWRITEBYTECODE="1", PYTHONPYCACHEPREFIX=str(session / "absent_python_cache"))
        os.environ.update(env)
        os.chdir(baseline["core_root"])
        checker = ProspectiveRuntimeChecker(plan, request["runtime_basis"], session / "lookup_checks")

        def check_all(specifications):
            assert_environment(dict(cwd=baseline["core_root"]), baseline)
            amendment(request["amendment"])
            for ref in [request_ref, request["policy"], *request["new_sources"]]:
                check(ref)
            runtime = checker(specifications)
            with fresh_cache(session / f"verification_cache_{checker.count}"):
                verify_environment(baseline)
            original_inputs(run)
            copied = prepared_inputs(run, baseline) if (root / "input").exists() else None
            return dict(runtime=runtime, original_inputs=run["inputs"], prepared_inputs=copied)

        def measurement(directory):
            prepare_inputs(run, baseline)
            guard = ReleaseGuard(request_ref, request, run, baseline)
            with fresh_tmp(root / "native_tmp"), fresh_cache(root / "native_numba_cache"):
                return measure(native_command(request["amendment"], run, baseline), directory, job,
                    32, MEMORY, 85800, 1., monitor_host=True, host_interval_s=30., release_guard=guard)

        outcome["wrapper"] = run_checked(lookup_binding["runtime_specs"], root, measurement, checker=check_all)
        outcome["status"] = "measurement_returned_pending_independent_review"
        if (root / "measurement/done.json").exists() and (root / "measurement/environment_preflight.json").exists():
            ref, reviewed = audit(root / "measurement", request["policy"],
                record(root / "measurement/environment_preflight.json"), job_id=job, index=12)
            outcome["environment_review"] = ref
            outcome["sampled_environment_evidence_valid"] = reviewed["sampled_environment_policy_satisfied"]
        if outcome["wrapper"]["status"] != "command_exited_zero":
            outcome["status"] = "factorial_attempt_failed_retained"
    except BaseException as error:
        outcome.update(status="factorial_executor_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        save(session / "result.json", outcome)
        os.chdir(previous_cwd)
    return outcome


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--job", type=int)
    parser.add_argument("--original-request", type=Path)
    parser.add_argument("--original-request-sha256")
    parser.add_argument("--composed-review", type=Path)
    parser.add_argument("--composed-review-sha256")
    parser.add_argument("--source-sha256")
    parser.add_argument("--request", type=Path)
    parser.add_argument("--request-sha256")
    args = parser.parse_args()
    if args.prepare:
        if (not all((args.job, args.original_request, args.original_request_sha256,
                     args.composed_review, args.composed_review_sha256, args.source_sha256))
                or args.request is not None or args.request_sha256 is not None):
            parser.error("Preparation requires held job, original request, composed review and source bindings only")
        original_ref, composed_ref = record(args.original_request), record(args.composed_review)
        require(original_ref["sha256"] == args.original_request_sha256
            and composed_ref["sha256"] == args.composed_review_sha256, "Preparation input digest differs")
        print(prepare(args.job, original_ref, composed_ref, args.source_sha256))
    else:
        if (not args.request or not args.request_sha256 or any(value is not None for value in (
                args.job, args.original_request, args.original_request_sha256,
                args.composed_review, args.composed_review_sha256, args.source_sha256))):
            parser.error("Execution requires request and digest only")
        ref = record(args.request)
        require(ref["sha256"] == args.request_sha256, "Native12 request checksum differs")
        result = execute(ref)
        raise SystemExit(0 if result["wrapper"]["status"] == "command_exited_zero" else 1)


if __name__ == "__main__":
    main()
