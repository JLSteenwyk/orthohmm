"""Execute one externally reviewed local scaling identity; never submit or retry."""

import argparse
import json
import os
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.bind_threadripper_panel_history import bind
from benchmark_tools.check_threadripper_runtime import RuntimeChecker
from benchmark_tools.measure_threadripper_run import measure_run
from benchmark_tools.measure_threadripper_scaling import measure
from benchmark_tools.measure_threadripper_boundary import measure as measure_boundary
from benchmark_tools.audit_threadripper_overhead import design as overhead_design
from benchmark_tools.manage_threadripper_environment_worker import EnvironmentWorker
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save, wait_file
from benchmark_tools.run_simulation_methods import read_frozen, execution_environment
from benchmark_tools.threadripper_panel_progress import overhead_position, position
from benchmark_tools.verify_threadripper_controller import ReleaseBudgetGuard
from benchmark_tools.review_threadripper_process_stream import audit as audit_process_stream
from benchmark_tools.review_threadripper_process_policy import number, shared_environment, SHARED_SCOPE
from benchmark_tools.review_threadripper_pressure_stream import native_pressure_role

PLAN_SHA = "c384e27730e3802b39ba14a42f7f50e84da5ce6deb9de9b2c32a74a745aed296"
LOOKUP_SHA = "d5f26d4c31346f6d710ca911a5220a0448a0a58e1997d708ceb821a8e2a18e34"
PRIVATE_PLAN_SHA = "a9358ac4f3c2f3eb9c1d7ce6a32528f8dd4bffef25d43cf95ac9765c5f2f3d05"
PRIVATE_LOOKUP_SHA = "5996f36dad39c7a7e38c38f134cbe5e13796c4a23dae0444a77a473f12cfdb54"
TIME_LIMIT = "1-02:00:00"


def deployment(request):
    """Select an explicit retained deployment, without changing historical pins."""
    name = request.get("deployment", "shared_v3_20260928")
    scope = request.get("execution_scope", "isolated_controlled")
    if not isinstance(scope, str) or scope not in {"isolated_controlled", SHARED_SCOPE}:
        raise ValueError("Unsupported execution scope")
    shared = scope == SHARED_SCOPE
    if shared and (name != "private_v2_20260928" or "runtime_lookup" not in request):
        raise ValueError("Shared-host timing requires explicit private runtime lookup")
    if "overhead_plan" in request and (name != "private_v2_20260928" or "runtime_lookup" not in request):
        raise ValueError("Engineering overhead requires explicit private runtime lookup")
    if "runtime_lookup" in request and name != "private_v2_20260928":
        raise ValueError("Explicit runtime lookup requires the retained private deployment")
    if name == "shared_v3_20260928":
        return dict(name=name, plan_sha=PLAN_SHA, lookup_sha=LOOKUP_SHA,
                    plan="threadripper_scaling_commands_20260928.json",
                    lookup="threadripper_python_lookup_v3_20260928.json",
                    script="run_threadripper_scaling.sh")
    if name == "private_v2_20260928":
        selected = dict(name=name, plan_sha=PRIVATE_PLAN_SHA, lookup_sha=PRIVATE_LOOKUP_SHA,
                    plan="threadripper_private_commands_20260928.json",
                    lookup="threadripper_private_lookup_v2_20260928.json",
                    script="run_threadripper_private_scaling.sh")
        if "runtime_lookup" in request:
            ref = request["runtime_lookup"]
            if (not isinstance(ref, dict) or set(ref) != {"path", "bytes", "sha256"}
                    or not isinstance(ref["path"], str)
                    or type(ref["bytes"]) is not int or ref["bytes"] <= 0
                    or not isinstance(ref["sha256"], str) or len(ref["sha256"]) != 64
                    or any(c not in "0123456789abcdef" for c in ref["sha256"])):
                raise ValueError("Require an externally pinned runtime lookup file")
            path = Path(ref["path"])
            if not path.is_absolute() or path.resolve() != path or path.is_symlink():
                raise ValueError("Require a direct absolute runtime lookup path")
            selected.update(lookup=ref["path"], lookup_sha=ref["sha256"])
        if "overhead_plan" in request:
            ref = request["overhead_plan"]
            if (not isinstance(ref, dict) or set(ref) != {"path", "bytes", "sha256"}
                    or not isinstance(ref["path"], str)
                    or type(ref["bytes"]) is not int or ref["bytes"] <= 0
                    or not isinstance(ref["sha256"], str) or len(ref["sha256"]) != 64
                    or any(c not in "0123456789abcdef" for c in ref["sha256"])):
                raise ValueError("Require an externally pinned engineering plan")
            path = Path(ref["path"])
            if not path.is_absolute() or path.resolve() != path or path.is_symlink():
                raise ValueError("Require a direct absolute engineering plan path")
            selected.update(parent_plan=selected["plan"], parent_plan_sha=selected["plan_sha"],
                            plan=ref["path"], plan_sha=ref["sha256"], purpose="native_overhead")
        if shared:
            selected.update(script="run_threadripper_shared_scaling.sh", execution_scope=SHARED_SCOPE,
                            allocation_mode="shared")
        return selected
    raise ValueError("Unknown retained Threadripper deployment")


def read(ref):
    path = Path(ref["path"])
    if not path.is_absolute() or path.resolve() != path or path.is_symlink():
        raise ValueError("Require direct absolute evidence paths")
    check(ref)
    data = json.loads(path.read_text())
    check(ref)
    return data


def expect(data, fields):
    if any(type(data.get(k)) is not type(v) or data[k] != v for k, v in fields.items()):
        raise ValueError("Execution evidence identity or decision differs")


def review(ref, expected):
    data = read(ref)
    expect(data, expected)
    if (not isinstance(data.get("review_reference"), str) or not data["review_reference"].strip()
            or not isinstance(data.get("evidence"), list) or not data["evidence"]):
        raise ValueError("Require explicit review reference and supporting evidence")
    for item in data["evidence"]:
        check(item)
    return data


def select(request, root, job):
    selected = deployment(request)
    overhead = "overhead_plan" in request
    expect(request, dict(schema="threadripper_execution_request_v1", job_id=job,
                         execution_authorized=True, plan_sha256=selected["plan_sha"],
                         lookup_sha256=selected["lookup_sha"]))
    index = request.get("index")
    if type(index) is not int or not 0 <= index < (54 if overhead else 27):
        raise ValueError("Invalid panel index")
    if (request["allocation_cwd"] != str(root)
            or request["scheduler_command"] != str(root / "benchmark_tools" / selected["script"])):
        raise ValueError("Wrong local allocation entry point")
    recipe = read(request["recipe"])
    expect(recipe, dict(schema="threadripper_executor_recipe_v1", root=str(root)))
    sources = recipe["sources"]
    required = set((root / "benchmark_tools").glob("*.py")) | {Path(request["scheduler_command"])}
    if (len(sources) != len(required) or {Path(r["path"]) for r in sources} != required):
        raise ValueError("Require complete executor Python and submission source inventory")
    for item in sources:
        check(item)
    results = root / "benchmark_tools/results"
    plan_path = results / selected["plan"]
    plan = read_frozen(plan_path, selected["plan_sha"])
    plan_evidence = []
    if overhead:
        if read(request["overhead_plan"]) != plan:
            raise ValueError("Engineering plan reference differs")
        parent_path = results / selected["parent_plan"]
        parent = read_frozen(parent_path, selected["parent_plan_sha"])
        overhead_design(plan, parent)
        parent_ref = record(parent_path)
        if not plan.get("sources") or plan["sources"][0] != parent_ref:
            raise ValueError("Engineering plan must bind the frozen private parent")
        plan_evidence = [request["overhead_plan"], parent_ref, *plan["sources"], *plan["helpers"]]
        for ref in plan_evidence:
            check(ref)
        overhead_position(plan["runs"], [])
    else:
        position(plan["runs"], [])
    lookup_path = results / selected["lookup"]
    lookup = read_frozen(lookup_path, selected["lookup_sha"])
    lookup_evidence = []
    if "runtime_lookup" in request:
        explicit = read(request["runtime_lookup"])
        expect(explicit, dict(status="native_lookup_repeated_identity_match",
                              scientific_execution_authorized=False))
        retained_path = results / "threadripper_private_lookup_v2_20260928.json"
        retained = read_frozen(retained_path, PRIVATE_LOOKUP_SHA)
        if explicit != lookup or explicit.get("baseline") != retained["baseline"]:
            raise ValueError("Explicit runtime lookup must preserve the frozen private baseline")
        binding = read(explicit["binding"])
        retained_binding = read(retained["binding"])
        if binding.get("controller_python") != retained_binding["controller_python"]:
            raise ValueError("Explicit runtime lookup must preserve the private controller")
        if any(binding.get(key) != retained_binding[key]
               for key in ("baseline", "command_plan", "baseline_paths", "retired_roots")):
            raise ValueError("Explicit runtime lookup must preserve scientific deployment bindings")
        specs = binding.get("runtime_specs")
        retained_specs = retained_binding["runtime_specs"]
        if (not isinstance(specs, list) or len(specs) != len(retained_specs)
                or specs[1:] != retained_specs[1:]):
            raise ValueError("Explicit runtime lookup must preserve private runtime manifests")
        lookup_evidence = [request["runtime_lookup"], record(retained_path),
                           retained["binding"], explicit["binding"], explicit["baseline"]]
    extra = dict(overhead=True) if overhead else {}
    shared = selected.get("execution_scope") == SHARED_SCOPE
    if shared:
        extra["allocation_mode"] = "shared"
    history = bind(record(plan_path), request["history"], command=request["scheduler_command"],
                   cwd=str(root), time_limit=TIME_LIMIT, **extra)
    if (history["progress"]["status"] != "next_identity_requires_preflight"
            or history["progress"]["index"] != index):
        raise ValueError("Panel is live, unresolved, complete or at another identity")
    ready_fields = dict(schema="threadripper_overhead_readiness_review_v1" if overhead else "threadripper_readiness_review_v1",
        decision="passed", plan_sha256=selected["plan_sha"], lookup_sha256=selected["lookup_sha"],
        recipe_sha256=request["recipe"]["sha256"], full_scale_observer_validated=True,
        environment_policy_frozen=True)
    if overhead:
        ready_fields.update(boundary_control_validated=True, incremental_overhead_measurement_pending=True)
    if shared:
        ready_fields.pop("full_scale_observer_validated")
        ready_fields.update(schema="threadripper_shared_readiness_review_v1",
            execution_scope=SHARED_SCOPE, observer_accounting_validated=True,
            contention_annotation_required=True, isolation_required=False)
    ready = review(request["readiness_review"], ready_fields)
    if shared and not shared_environment(read(ready["environment_policy"])):
        raise ValueError("Shared readiness must bind a shared-host environmental policy")
    if overhead:
        task = plan["runs"][index]
        history["engineering_task"] = {k: task[k] for k in ("index", "pair", "arm", "method", "proteomes", "repeat")}
        run = task["run"]
    else:
        run = plan["runs"][index]
    session = Path(run["measurement_directory"]).parent.parent / "sessions" / f"run_{index:02d}"
    if request["environment_preflight_path"] != str(session / "environment_preflight.json"):
        raise ValueError("Require run-specific environmental review location")
    return run, lookup, history, [*sources, *ready["evidence"], *lookup_evidence, *plan_evidence]


class EnvironmentalReleaseGuard:
    """Bind a separately reviewed preflight, then perform the live budget check.

    This validates a review's identity and freshness, not its conclusions.
    The review must be supplied by a workflow that inspects actual host evidence.
    """

    def __init__(self, request, request_ref, budget, *, evidence=(), clock=time.time_ns, waiter=wait_file):
        self.request, self.request_ref, self.budget, self.clock = request, request_ref, budget, clock
        self.evidence = evidence
        self.waiter = waiter

    def __call__(self, directory):
        report = dict(status="environment_release_failed", scientific_timings_admitted=False)
        try:
            check(self.request_ref)
            for item in self.evidence:
                check(item)
            path = Path(self.request["environment_preflight_path"])
            if path.exists() or path.is_symlink():
                raise FileExistsError("Environmental review must follow this release request")
            requested = self.clock()
            save(directory / "environment_review_requested.json", dict(
                job_id=self.request["job_id"], index=self.request["index"],
                request=self.request_ref, requested_unix_ns=requested,
                review_path=str(path), wait_seconds=20, native_released=False))
            # The parked native worker's unchanged gate timeout is 45 seconds.
            self.waiter(path, seconds=20)
            ref = record(path)
            expected = dict(schema="threadripper_environment_preflight_v1", decision="passed",
                job_id=self.request["job_id"], index=self.request["index"], plan_sha256=deployment(self.request)["plan_sha"],
                recipe_sha256=self.request["recipe"]["sha256"],
                readiness_review_sha256=self.request["readiness_review"]["sha256"],
                whole_run_observer_ready=True)
            if self.request.get("execution_scope") == SHARED_SCOPE:
                expected.update(execution_scope=SHARED_SCOPE, background_competition_recorded=True,
                                uncontended_timing=False, foreign_cpu_used_for_eligibility=False)
            else:
                expected["unrelated_scientific_work_present"] = False
            ready = review(ref, expected)
            start, end = ready["observation_started_unix_ns"], ready["observation_finished_unix_ns"]
            now = self.clock()
            if (type(start) is not int or type(end) is not int or not 0 < requested <= start <= end <= now
                    or now - start > 120_000_000_000):
                raise ValueError("Environmental observation is invalid or older than 120 seconds")
            if ready["boot_id"] != Path("/proc/sys/kernel/random/boot_id").read_text().strip():
                raise ValueError("Environmental observation belongs to another boot")
            report.update(status="environment_review_bound", review=ref, checked_unix_ns=now,
                          observational_validity_independently_established=False)
            result = self.budget(directory)
            check(ref)
            if not end <= self.clock() <= start + 120_000_000_000:
                raise ValueError("Environmental review expired during release checks")
            return result
        except BaseException as error:
            report.update(status="environment_release_failed", error_type=type(error).__name__, error=str(error))
            raise
        finally:
            save(directory / "environment_release.json", report)


def execute(request_path, request_sha):
    root = Path(__file__).resolve().parent.parent
    request_ref = record(request_path)
    if request_ref["sha256"] != request_sha:
        raise ValueError("Execution request digest differs")
    request = read(request_ref)
    selected = deployment(request)
    if os.environ.get("ORTHOHMM_THREADRIPPER_DEPLOYMENT", "shared_v3_20260928") != selected["name"]:
        raise ValueError("Submission bootstrap and requested deployment differ")
    job = int(os.environ["SLURM_JOB_ID"])
    cache = Path(f"/dev/shm/orthohmm_scaling_driver_{job}")
    if (os.uname().nodename != "bizon" or Path.cwd() != root
            or os.environ.get("SLURM_CPUS_PER_TASK") != "64"
            or os.environ.get("SLURM_MEM_PER_NODE") != "131072"
            or not sys.dont_write_bytecode or sys.pycache_prefix != str(cache)
            or cache.exists() or cache.is_symlink() or os.environ.get("PYTHONHASHSEED") != "0"):
        raise ValueError("Require local matched allocation, repository cwd and disabled bytecode writes")
    run, lookup, history, sources = select(request, root, job)
    readiness = read(request["readiness_review"])
    policy_ref = readiness["environment_policy"]
    policy = review(policy_ref, dict(decision="reviewed",
                                    host="bizon", plan_sha256=selected["plan_sha"]))
    pressure_role = native_pressure_role(policy)
    shared = shared_environment(policy)
    if shared != (request.get("execution_scope") == SHARED_SCOPE):
        raise ValueError("Request and policy execution scopes differ")
    if selected["name"] == "private_v2_20260928" and pressure_role != "diagnostic_only":
        raise ValueError("Private timing requires explicit v2 diagnostic-only native pressure policy")
    number(policy["maximum_foreign_average_cores"])
    number(policy["maximum_sample_period_s"], positive=True)
    number(policy["maximum_pressure_sample_period_s"], positive=True)
    limits = policy["maximum_pressure_percent"]
    if set(limits) != {"cpu", "memory", "io"} or any(number(v) > 100 for v in limits.values()):
        raise ValueError("Require prospective CPU, memory and I/O pressure bounds")
    baseline = read(lookup["baseline"])
    binding = read(lookup["binding"])
    if selected["name"] == "private_v2_20260928":
        controller = binding.get("controller_python")
        if not isinstance(controller, dict) or os.path.abspath(sys.executable) != controller["path"]:
            raise ValueError("Private deployment requires its bound controller interpreter")
        check(controller)
    session = Path(run["measurement_directory"]).parent.parent / "sessions" / f"run_{run['index']:02d}"
    session.mkdir(parents=True, exist_ok=False)
    result = dict(status="executor_started", index=run["index"], job_id=job,
                  request=request_ref, history=history, deployment=selected, automatic_retry=False,
                  scientific_timings_admitted=False, next_submission_authorized=False,
                  limitations=["External reviews are bound, not independently certified by this executor.",
                               "Terminal scheduler, runtime, environment, resource and native-output audits remain required.",
                               "A returned measurement never authorizes the next identity."])
    if shared:
        result.update(execution_scope=SHARED_SCOPE, uncontended_timing=False,
                      contention_distortion="unknown_potentially_method_dependent")
    task = history.get("engineering_task")
    if task:
        result.update(purpose="native_overhead", engineering_task=task, production_identity=False)
    save(session / "started.json", result)
    try:
        env, _ = execution_environment(baseline)
        env.update(PYTHONDONTWRITEBYTECODE="1", PYTHONPYCACHEPREFIX=str(session / "python_cache"))
        os.environ.update(env)
        os.chdir(run["cwd"])
        checker = RuntimeChecker(root / "benchmark_tools/results" / selected["lookup"],
            selected["lookup_sha"], session / "lookup_checks")
        budget_extra = dict(allocation_mode="shared") if shared else {}
        budget = ReleaseBudgetGuard(job, command=request["scheduler_command"], cwd=str(root), **budget_extra)
        evidence = [request["recipe"], request["readiness_review"], policy_ref, *sources, *history["evidence"]]
        with EnvironmentWorker(session, request_ref, policy_ref, root) as worker:
            # Finish observer setup before parking a native worker on its short gate.
            prepared_ref = worker.wait_prepared()
            def finished_worker_budget(directory):
                worker.finish()
                return budget(directory)
            guard = EnvironmentalReleaseGuard(request, request_ref, finished_worker_budget,
                evidence=[*evidence, prepared_ref], waiter=worker.wait_response)
            collector = measure_boundary if task and task["arm"] == "boundary" else measure
            result["wrapper"] = measure_run(run, baseline, binding["runtime_specs"], collector, job,
                                            runtime_checker=checker, release_guard=guard)
        stream_extra = dict(collector_arm=task["arm"]) if task else {}
        stream_ref, stream_review = audit_process_stream(run["measurement_directory"], policy_ref,
            record(request["environment_preflight_path"]), job_id=job, index=run["index"], **stream_extra)
        result["process_stream_review"] = stream_ref
        if not stream_review["sampled_environment_policy_satisfied"]:
            raise ValueError("Whole-run sampled environment policy was not satisfied")
        for item in [request_ref, *evidence]:
            check(item)
        result["status"] = "measurement_returned_pending_independent_review"
    except BaseException as error:
        result.update(status="executor_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        lifecycle = session / "environment_worker_lifecycle.json"
        if lifecycle.is_file():
            result["environment_worker"] = record(lifecycle)
        save(session / "result.json", result)
        os.chdir(root)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request", type=Path, required=True)
    parser.add_argument("--request-sha256", required=True)
    args = parser.parse_args()
    result = execute(args.request.absolute(), args.request_sha256)
    raise SystemExit(0 if result["wrapper"]["status"] == "command_exited_zero" else 1)
