"""Execute one new frozen identity through the allocation-aware placement route.

Separate from the hash-bound historical controller. Scientific entrypoint,
arguments, inputs, runtime and accounting helpers are reused unchanged.
"""

import argparse
import importlib
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from benchmark_tools.native_factorial_allocated_execution import (
    SCRIPT, SCOPE, MEMORY, amendment, validate_request, native_command, reviewed_history)
from benchmark_tools.native_factorial_allocated_placement import validate as validate_placement
from benchmark_tools.native_factorial_adapter import entrypoint
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.run_native_factorial_cost import (
    CORE_COMMIT, read, prepare_inputs, prepared_inputs, original_inputs, native_kwargs,
    verify_metrics, update_native_receipt, available_memory)
from benchmark_tools.verify_lineage_native_provenance import same


def native_placement(ready, current, job, parent_pid):
    from benchmark_tools.probe_threadripper_allocation import validate
    allowed = validate_placement(ready["allocated_placement"], job)
    validate(ready["placement"], current, 64)
    if (not same(ready["placement"], ready["allocated_placement"]["bound"])
            or type(ready["pid"]) is not int or ready["pid"] != ready["placement"]["pid"]
            or type(parent_pid) is not int or ready["pid"] != parent_pid
            or type(current["pid"]) is not int or current["pid"] <= 0 or current["pid"] == ready["pid"]
            or ready["cgroup"] != current["cgroup"]
            or not same(ready["placement"]["ancestors"], current["ancestors"])
            or current["affinity"] != allowed
            or current["slurm"]["SLURM_JOB_ID"] != str(job)):
        raise ValueError("Native child did not inherit its actual bound worker placement")
    return allowed


def native(amendment_ref, index):
    from benchmark_tools.probe_threadripper_allocation import inspect
    execution, plan = amendment(amendment_ref)
    if type(index) is not int or index not in execution["allowed_indices"]:
        raise ValueError("Only genuinely unrun native identities are supported")
    plan_ref, run = execution["historical_plan"], plan["runs"][index]
    baseline = read(plan["baseline"])
    root = Path(run["output_root"])
    ready_ref = record(root / "measurement/ready.json")
    ready = read(ready_ref)
    job = int(os.environ["SLURM_JOB_ID"])
    placement = inspect()
    allowed = native_placement(ready, placement, job, os.getppid())
    if (not sys.dont_write_bytecode or not sys.pycache_prefix
            or Path(sys.pycache_prefix).exists() or Path(sys.pycache_prefix).is_symlink()
            or any(os.environ.get(k) != v for k, v in baseline["environment_overrides"].items())
            or os.path.abspath(sys.executable) != baseline["tool_entrypoints"]["orthohmm_python"]["absolute_path"]):
        raise ValueError("Native bytecode, environment or private interpreter differs")
    core = Path(baseline["core_root"])
    sys.path.insert(0, str(core))
    module = importlib.import_module("orthohmm.orthohmm")
    if (Path(module.__file__).resolve() != core / "orthohmm/orthohmm.py"
            or module.fetch_fasta_files(str(root / "input")) != run["native_order"]):
        raise ValueError("Frozen pipeline origin or native enumeration changed")
    adapted, flags = entrypoint(module, run["cell"])
    if (root / "native").exists() or (root / "native").is_symlink():
        raise FileExistsError("Native output must be fresh")
    receipt = dict(schema="allocated_native_factorial_execution_v1", status="native_factorial_running",
        index=index, cell=run["cell"], plan=plan_ref, amendment=amendment_ref,
        factors=flags, placement=placement, parent_pid=os.getppid(), allocated_ready=ready_ref, native_cpu_ids=allowed,
        native_order=run["native_order"], source=record(__file__), automatic_retry=False,
        cold_cache_claim=False, accuracy_evaluated=False)
    path = root / "native_execution.json"
    save(path, receipt)
    initial = dict(receipt)
    try:
        try:
            adapted(**native_kwargs(module, run, baseline))
        except SystemExit as error:
            if error.code not in (None, 0):
                raise
        metrics = json.loads((root / "metrics.json").read_text())
        verify_metrics(metrics, run, flags)
        receipt.update(status="native_factorial_completed_pending_output_review", stages=sorted(metrics["stages"]),
            counts=metrics["counts"], exact_historical_output_equivalence_established=False)
    except BaseException as error:
        receipt.update(status="native_factorial_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        check(ready_ref)
        check(amendment_ref)
        update_native_receipt(path, initial, receipt)


class ReleaseGuard:
    def __init__(self, request_ref, amendment_ref, execution, plan, run, baseline):
        self.request_ref, self.amendment_ref, self.execution = request_ref, amendment_ref, execution
        self.plan, self.run, self.baseline = plan, run, baseline

    def __call__(self, directory):
        from benchmark_tools.probe_host_counters import snapshot
        from benchmark_tools.observe_threadripper_process_identity import enriched_snapshot
        from benchmark_tools.review_threadripper_process_policy import review
        from benchmark_tools.slurm_resource_snapshot import scoped_path
        from benchmark_tools.verify_threadripper_controller import ReleaseBudgetGuard
        request = read(self.request_ref)
        execution, _ = amendment(self.amendment_ref)
        if not same(execution, self.execution):
            raise ValueError("Execution amendment changed before release")
        job = request["job_id"]
        validate_request(request, self.amendment_ref, execution, job)
        prepared_inputs(self.run, self.baseline)
        ready = read(record(directory / "ready.json"))
        allowed = validate_placement(ready["allocated_placement"], job)
        if ready["placement"]["affinity"] != allowed:
            raise ValueError("Ready worker affinity differs")
        scope = scoped_path(ready["cgroup"], job)
        job_scope = next(p for p in scope.parents if p.name == f"job_{job}")
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
        if capacity < MEMORY or not verdict["process_policy_matched"]:
            raise ValueError("Unsafe available RAM or invalid process attribution")
        budget = ReleaseBudgetGuard(job, command=str(SCRIPT), cwd=str(ROOT), allocation_mode="shared")(directory)
        if (budget["allocation"]["fields"].get("Comment") != self.request_ref["sha256"]
                or budget["allocation"]["fields"].get("JobName") != "orthohmm_allocated_factorial"):
            raise ValueError("Scheduler differs from prepared allocated request")
        save(directory / "environment_preflight.json", dict(schema="threadripper_environment_preflight_v1",
            decision="passed", job_id=job, index=self.run["index"], plan_sha256=request["plan"]["sha256"],
            amendment=self.amendment_ref, environment_policy=request["policy"], execution_scope=SCOPE,
            boot_id=boot, background_competition_recorded=True, uncontended_timing=False,
            evidence=[record(directory / "launch_environment_observation.json"), record(directory / "release_budget.json")]))
        return budget


def execute(request_ref):
    from benchmark_tools.check_threadripper_runtime import RuntimeChecker
    from benchmark_tools.isolated_numba_cache import fresh_cache
    from benchmark_tools.isolated_native_tmp import fresh_tmp
    from benchmark_tools.measure_native_scaling_run import run_checked
    from benchmark_tools.measure_threadripper_run import assert_environment
    from benchmark_tools.measure_allocated_threadripper_scaling import measure
    from benchmark_tools.review_threadripper_process_stream import audit, shared_environment
    from benchmark_tools.run_simulation_methods import execution_environment, verify_environment
    request = read(request_ref)
    amendment_ref = request["amendment"]
    execution, plan = amendment(amendment_ref)
    plan_ref = execution["historical_plan"]
    job = int(os.environ["SLURM_JOB_ID"])
    validate_request(request, amendment_ref, execution, job)
    if (Path.cwd() != ROOT or os.uname().nodename != "bizon" or not sys.dont_write_bytecode
            or not sys.pycache_prefix or Path(sys.pycache_prefix).exists() or Path(sys.pycache_prefix).is_symlink()
            or os.environ.get("PYTHONHASHSEED") != "0" or os.environ.get("PYTHONNOUSERSITE") != "1"):
        raise ValueError("Wrong controller host, cwd or bytecode/environment bootstrap")
    lookup = read(plan["runtime_lookup"])
    binding = read(lookup["binding"])
    baseline = read(plan["baseline"])
    if lookup["baseline"] != plan["baseline"] or baseline["core_commit"] != CORE_COMMIT:
        raise ValueError("Runtime lookup/core baseline differs")
    if os.path.abspath(sys.executable) != binding["controller_python"]["path"]:
        raise ValueError("Wrong controller interpreter")
    check(binding["controller_python"])
    policy = read(request["policy"])
    if (not shared_environment(policy) or policy["plan_sha256"] != plan_ref["sha256"]
            or policy["minimum_available_memory_bytes"] != MEMORY):
        raise ValueError("Policy does not bind this shared-host plan")
    history_evidence = reviewed_history(request, amendment_ref, execution, plan)
    for ref in [plan_ref, request_ref, request["policy"], *plan["evidence"]]:
        check(ref)
    run = plan["runs"][request["index"]]
    root = Path(run["output_root"])
    session = Path(plan["panel_root"]) / "sessions" / f"run_{run['index']:02d}"
    session.mkdir(parents=True, exist_ok=False)
    outcome = dict(schema="allocated_native_factorial_session_v1", status="factorial_executor_started",
        job_id=job, index=run["index"], cell=run["cell"], request=request_ref, plan=plan_ref,
        amendment=amendment_ref, source=record(__file__), execution_scope=SCOPE, automatic_retry=False,
        scientific_timings_admitted=False, next_identity_authorized=False, uncontended_timing=False,
        contention_distortion="unknown_potentially_method_dependent", history_scheduler_evidence=history_evidence)
    save(session / "started.json", outcome)
    env, _ = execution_environment(baseline)
    env.update(PYTHONDONTWRITEBYTECODE="1", PYTHONPYCACHEPREFIX=str(session / "absent_python_cache"))
    os.environ.update(env)
    os.chdir(baseline["core_root"])
    checker = RuntimeChecker(Path(plan["runtime_lookup"]["path"]), plan["runtime_lookup"]["sha256"], session / "lookup_checks")
    context = dict(cwd=baseline["core_root"])

    def check_all(specifications):
        assert_environment(context, baseline)
        amendment(amendment_ref)
        for ref in [request_ref, request["policy"]]:
            check(ref)
        runtime = checker(specifications)
        with fresh_cache(session / f"verification_cache_{checker.count}"):
            verify_environment(baseline)
        original_inputs(run)
        copied = prepared_inputs(run, baseline) if (root / "input").exists() else None
        return dict(runtime=runtime, original_inputs=run["inputs"], prepared_inputs=copied)

    def measurement(directory):
        prepare_inputs(run, baseline)
        guard = ReleaseGuard(request_ref, amendment_ref, execution, plan, run, baseline)
        with fresh_tmp(root / "native_tmp"), fresh_cache(root / "native_numba_cache"):
            return measure(native_command(amendment_ref, run, baseline), directory, job, 32, MEMORY,
                85800, 1., monitor_host=True, host_interval_s=30., release_guard=guard)

    try:
        outcome["wrapper"] = run_checked(binding["runtime_specs"], root, measurement, checker=check_all)
        outcome["status"] = "measurement_returned_pending_independent_review"
        if (root / "measurement/done.json").exists() and (root / "measurement/environment_preflight.json").exists():
            ref, review = audit(root / "measurement", request["policy"], record(root / "measurement/environment_preflight.json"),
                                job_id=job, index=run["index"])
            outcome["environment_review"] = ref
            outcome["sampled_environment_evidence_valid"] = review["sampled_environment_policy_satisfied"]
        if outcome["wrapper"]["status"] != "command_exited_zero":
            outcome["status"] = "factorial_attempt_failed_retained"
    except BaseException as error:
        outcome.update(status="factorial_executor_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        save(session / "result.json", outcome)
        os.chdir(ROOT)
    return outcome


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", action="store_true")
    parser.add_argument("--amendment", type=Path)
    parser.add_argument("--amendment-sha256")
    parser.add_argument("--index", type=int)
    parser.add_argument("--request", type=Path)
    parser.add_argument("--request-sha256")
    args = parser.parse_args()
    if args.native:
        if args.request is not None or args.request_sha256 or args.amendment is None or not args.amendment_sha256 or args.index is None:
            parser.error("Native mode requires only amendment, digest and index")
        ref = record(args.amendment)
        if ref["sha256"] != args.amendment_sha256:
            raise ValueError("Native amendment checksum differs")
        native(ref, args.index)
    else:
        if args.amendment is not None or args.amendment_sha256 or args.index is not None or args.request is None or not args.request_sha256:
            parser.error("Controller mode requires only request and digest")
        ref = record(args.request)
        if ref["sha256"] != args.request_sha256:
            raise ValueError("Request checksum differs")
        result = execute(ref)
        raise SystemExit(0 if result["wrapper"]["status"] == "command_exited_zero" else 1)


if __name__ == "__main__":
    main()
