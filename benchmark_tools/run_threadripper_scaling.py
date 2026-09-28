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
from benchmark_tools.manage_threadripper_environment_worker import EnvironmentWorker
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save, wait_file
from benchmark_tools.run_simulation_methods import read_frozen, execution_environment
from benchmark_tools.threadripper_panel_progress import position
from benchmark_tools.verify_threadripper_controller import ReleaseBudgetGuard
from benchmark_tools.review_threadripper_process_stream import audit as audit_process_stream
from benchmark_tools.review_threadripper_process_policy import number

PLAN_SHA = "c384e27730e3802b39ba14a42f7f50e84da5ce6deb9de9b2c32a74a745aed296"
LOOKUP_SHA = "d5f26d4c31346f6d710ca911a5220a0448a0a58e1997d708ceb821a8e2a18e34"
TIME_LIMIT = "1-02:00:00"


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
    expect(request, dict(schema="threadripper_execution_request_v1", job_id=job,
                         execution_authorized=True, plan_sha256=PLAN_SHA,
                         lookup_sha256=LOOKUP_SHA))
    index = request.get("index")
    if type(index) is not int or not 0 <= index < 27:
        raise ValueError("Invalid panel index")
    if (request["allocation_cwd"] != str(root)
            or request["scheduler_command"] != str(root / "benchmark_tools/run_threadripper_scaling.sh")):
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
    plan_path = results / "threadripper_scaling_commands_20260928.json"
    plan = read_frozen(plan_path, PLAN_SHA)
    position(plan["runs"], [])
    lookup_path = results / "threadripper_python_lookup_v3_20260928.json"
    lookup = read_frozen(lookup_path, LOOKUP_SHA)
    history = bind(record(plan_path), request["history"], command=request["scheduler_command"],
                   cwd=str(root), time_limit=TIME_LIMIT)
    if (history["progress"]["status"] != "next_identity_requires_preflight"
            or history["progress"]["index"] != index):
        raise ValueError("Panel is live, unresolved, complete or at another identity")
    ready = review(request["readiness_review"], dict(schema="threadripper_readiness_review_v1",
        decision="passed", plan_sha256=PLAN_SHA, lookup_sha256=LOOKUP_SHA,
        recipe_sha256=request["recipe"]["sha256"], full_scale_observer_validated=True,
        environment_policy_frozen=True))
    run = plan["runs"][index]
    session = Path(run["measurement_directory"]).parent.parent / "sessions" / f"run_{index:02d}"
    if request["environment_preflight_path"] != str(session / "environment_preflight.json"):
        raise ValueError("Require run-specific environmental review location")
    return run, lookup, history, [*sources, *ready["evidence"]]


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
            ready = review(ref, dict(schema="threadripper_environment_preflight_v1", decision="passed",
                job_id=self.request["job_id"], index=self.request["index"], plan_sha256=PLAN_SHA,
                recipe_sha256=self.request["recipe"]["sha256"],
                readiness_review_sha256=self.request["readiness_review"]["sha256"],
                whole_run_observer_ready=True, unrelated_scientific_work_present=False))
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
    policy = review(policy_ref, dict(schema="threadripper_environment_policy_v1", decision="reviewed",
                                    host="bizon", plan_sha256=PLAN_SHA))
    number(policy["maximum_foreign_average_cores"])
    number(policy["maximum_sample_period_s"], positive=True)
    baseline = read(lookup["baseline"])
    binding = read(lookup["binding"])
    session = Path(run["measurement_directory"]).parent.parent / "sessions" / f"run_{run['index']:02d}"
    session.mkdir(parents=True, exist_ok=False)
    result = dict(status="executor_started", index=run["index"], job_id=job,
                  request=request_ref, history=history, automatic_retry=False,
                  scientific_timings_admitted=False, next_submission_authorized=False,
                  limitations=["External reviews are bound, not independently certified by this executor.",
                               "Terminal scheduler, runtime, environment, resource and native-output audits remain required.",
                               "A returned measurement never authorizes the next identity."])
    save(session / "started.json", result)
    try:
        env, _ = execution_environment(baseline)
        env.update(PYTHONDONTWRITEBYTECODE="1", PYTHONPYCACHEPREFIX=str(session / "python_cache"))
        os.environ.update(env)
        os.chdir(run["cwd"])
        checker = RuntimeChecker(root / "benchmark_tools/results/threadripper_python_lookup_v3_20260928.json",
            LOOKUP_SHA, session / "lookup_checks")
        budget = ReleaseBudgetGuard(job, command=request["scheduler_command"], cwd=str(root))
        evidence = [request["recipe"], request["readiness_review"], policy_ref, *sources, *history["evidence"]]
        with EnvironmentWorker(session, request_ref, policy_ref, root) as worker:
            def finished_worker_budget(directory):
                worker.finish()
                return budget(directory)
            guard = EnvironmentalReleaseGuard(request, request_ref, finished_worker_budget,
                evidence=evidence, waiter=worker.wait_response)
            result["wrapper"] = measure_run(run, baseline, binding["runtime_specs"], measure, job,
                                            runtime_checker=checker, release_guard=guard)
        stream_ref, stream_review = audit_process_stream(run["measurement_directory"], policy_ref,
            record(request["environment_preflight_path"]), job_id=job, index=run["index"])
        result["process_stream_review"] = stream_ref
        if not stream_review["sampled_process_policy_satisfied"]:
            raise ValueError("Whole-run sampled process policy was not satisfied")
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
