"""Execute one prospectively bound factorial cost identity; no submission/retry.

The native wrapper and controller are separate processes. Preparation and
runtime/input hashing remain outside the collector's native interval.
"""

import argparse
import csv
import hashlib
import importlib
import json
import os
from pathlib import Path
import shutil
import re
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from benchmark_tools.native_factorial_adapter import entrypoint, factors
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save

SCOPE = "shared_host_matched_resources"
MEMORY = 128 * 1024**3
CELLS = ("p0_c0_r0", "p0_c0_r1", "p0_c1_r0", "p0_c1_r1", "p1_c0_r1", "p1_c1_r0")
IDENTITIES = tuple(("orthobench", c) for c in CELLS) + tuple(
    ("qfo_corrected", c) for c in (*CELLS, "p1_c1_r1"))
CORE_COMMIT = "7f3a9e40dd7e79f842cc2c11fb8b548f9a802806"
SCRIPT = ROOT / "benchmark_tools/run_native_factorial_cost.sh"


def read(ref):
    path = Path(ref["path"])
    if not path.is_absolute() or path.resolve() != path or path.is_symlink():
        raise ValueError("Require a direct absolute evidence path")
    check(ref)
    result = json.loads(path.read_text())
    check(ref)
    return result


def validate_plan(plan):
    if (plan.get("schema") != "native_factorial_cost_plan_v1"
            or plan.get("root") != str(ROOT) or plan.get("execution_scope") != SCOPE
            or plan.get("core_commit") != CORE_COMMIT
            or plan.get("resources") != dict(native_cpu_ids=list(range(32)), slurm_slots=64,
                memory_bytes=MEMORY, timeout_s=85800, sample_period_s=1., host_period_s=30.,
                minimum_available_memory_bytes=MEMORY)
            or plan.get("automatic_retry") is not False):
        raise ValueError("Wrong factorial execution scope or resources")
    runs = plan.get("runs")
    if not isinstance(runs, list) or len(runs) != len(IDENTITIES):
        raise ValueError("Require all thirteen prospective identities")
    panel = Path(plan["panel_root"])
    if (not panel.is_absolute() or panel.resolve() != panel
            or not panel.is_relative_to(ROOT / "benchmarks/results") or panel.exists() and panel.is_symlink()):
        raise ValueError("Require a direct persistent factorial panel root")
    for i, (run, identity) in enumerate(zip(runs, IDENTITIES)):
        if (type(run.get("index")) is not int or run["index"] != i
                or (run.get("dataset"), run.get("cell")) != identity
                or type(run.get("repeat")) is not int or run["repeat"] != 0
                or run.get("output_root") != str(panel / f"run_{i:02d}")
                or run.get("genes") != (251378 if i < 6 else 984137)
                or run.get("proteomes") != (12 if i < 6 else 78)):
            raise ValueError("Factorial identity, order or universe differs")
        names = [Path(r["path"]).name for r in run["inputs"]]
        if (len(names) != run["proteomes"] or len(set(names)) != len(names)
                or run.get("input_creation_order") != sorted(names)
                or len(run.get("native_order", [])) != len(names) or set(run["native_order"]) != set(names)
                or any(Path(r["path"]).parent != Path(run["input_directory"]) for r in run["inputs"])):
            raise ValueError("Input ownership, inventory or prospective order differs")
    required = {str(ROOT / "benchmark_tools" / n) for n in (
        "run_native_factorial_cost.py", "prepare_native_factorial_cost.py", "native_factorial_adapter.py")}
    source_paths = [r["path"] for r in plan["helper_sources"]]
    if not required <= set(source_paths) or str(SCRIPT) not in source_paths or len(set(source_paths)) != len(source_paths):
        raise ValueError("Missing or duplicated execution source bindings")
    return runs


def validate_request(request, plan_ref, job):
    if (request.get("schema") != "native_factorial_cost_request_v1"
            or request.get("execution_authorized") is not True
            or type(request.get("job_id")) is not int or request["job_id"] != job
            or request.get("plan") != plan_ref or request.get("scheduler_command") != str(SCRIPT)
            or request.get("allocation_cwd") != str(ROOT)
            or type(request.get("index")) is not int or not 0 <= request["index"] < 13
            or not isinstance(request.get("history"), list)
            or len(request["history"]) != request["index"]):
        raise ValueError("Wrong explicit factorial request or sequential history")


def terminal_accounting(raw, job):
    from benchmark_tools.capture_array_scheduler import TERMINAL
    rows = list(csv.reader(raw.strip().splitlines(), delimiter="|"))
    if len(rows) != 1 or len(rows[0]) != 12:
        raise ValueError("Require one full non-array accounting row")
    fields = dict(zip(("JobIDRaw", "State", "ExitCode", "AllocCPUS", "ReqMem", "NodeList", "Partition",
                       "TimelimitRaw", "Submit", "Start", "End", "JobName"), rows[0]))
    if (fields["JobIDRaw"] != str(job) or fields["State"] not in TERMINAL
            or not re.fullmatch(r"[0-9]+:[0-9]+", fields["ExitCode"])
            or fields["AllocCPUS"] != "64" or fields["ReqMem"] not in {"128G", "128Gn", "131072M", "131072Mn"}
            or fields["NodeList"] != "bizon" or fields["Partition"] != "gpu"
            or fields["TimelimitRaw"] != "1560" or fields["JobName"] != "orthohmm_factorial_cost"
            or any(not re.fullmatch(r"[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}", fields[k])
                   for k in ("Submit", "Start", "End"))):
        raise ValueError("Accounting identity, resource envelope or terminal state differs")
    return fields


def verify_terminal(job):
    from benchmark_tools.verify_threadripper_controller import validate
    command = ["scontrol", "show", "job", str(job), "--oneliner"]
    result = subprocess.run(command, capture_output=True, text=True, timeout=5)
    observation = dict(command=command, returncode=result.returncode, stdout=result.stdout, stderr=result.stderr)
    if result.returncode == 0:
        verified = validate(result.stdout, job, "terminal", command=str(SCRIPT), cwd=str(ROOT),
                            time_limit="1-02:00:00", allocation_mode="shared")
        return dict(source="live_controller", observation=observation, verified=verified)
    if "Invalid job id specified" not in result.stderr:
        raise ValueError("Controller observation failed, not an expired terminal record")
    command = ["sacct", "-X", "-j", str(job), "-n", "-P", "--format=JobIDRaw,State,ExitCode,AllocCPUS,ReqMem,NodeList,Partition,TimelimitRaw,Submit,Start,End,JobName"]
    account = subprocess.run(command, capture_output=True, text=True, timeout=5, check=True)
    return dict(source="fresh_accounting_after_controller_expiry", controller_observation=observation,
        observation=dict(command=command, stdout=account.stdout, stderr=account.stderr),
        verified=terminal_accounting(account.stdout, job),
        limitation="Accounting confirms terminal identity/resources, not every historical controller field; prior bound request/review remain required.")


def native_kwargs(module, run, baseline):
    flag = factors(run["cell"])
    root = Path(run["output_root"])
    return dict(fasta_directory=str(root / "input"), output_directory=str(root / "native"),
        phmmer="unused", cpu=32, threads_per_worker=4, single_copy_threshold=.5,
        mcl="unused", inflation_value=1.5, start=None, stop=module.StopStep.infer,
        substitution_matrix=module.SubstitutionMatrix.blosum62, evalue_threshold=.0001,
        search_mode="builtin", clustering="leiden", cpm_resolution=.1,
        refinement_profile="default", accuracy_profile="high_sensitivity", metrics_json=str(root / "metrics.json"),
        phylogeny="reconcile" if flag["reconciliation"] else "off",
        phylogeny_candidates="satellite_v2" if flag["candidate_expansion"] else "seed",
        species_tree_mode="infer", species_tree_rooting="min_variance",
        phylogeny_root_rule="species_overlap", phylogeny_pair_rule="positive_paralogy",
        aligner=baseline["tool_entrypoints"]["mafft"]["absolute_path"],
        tree_builder=baseline["tool_entrypoints"]["FastTree"]["absolute_path"])


def verify_metrics(metrics, run, flags):
    expected = {"search", "edge_thresholds", "network_edges", "clustering", "refinement", "orthogroup_materialization"}
    for factor, stage in (("profile_expansion", "profile_expansion"),
                          ("candidate_expansion", "phylogeny_candidates"), ("reconciliation", "phylogeny")):
        if flags[factor]:
            expected.add(stage)
    if (metrics.get("status") != "complete" or metrics["metadata"].get("native_factorial") != flags
            or metrics["counts"].get("genes") != run["genes"] or set(metrics["stages"]) != expected):
        raise ValueError("Native stage/factor/universe evidence differs")
    if flags["reconciliation"] and (metrics["counts"].get("phylogeny_checkpoint_hits") != 0
            or metrics["counts"].get("phylogeny_species_tree_checkpoint_hit") is not False):
        raise ValueError("Native run reused reconciliation checkpoints")


def original_inputs(run):
    for ref in run["inputs"]:
        check(ref)


def prepared_inputs(run, baseline):
    from benchmark_tools.snapshot_orthohmm_input_order import snapshot
    target = Path(run["output_root"]) / "input"
    expected = {Path(r["path"]).name: r for r in run["inputs"]}
    entries = list(target.iterdir())
    if set(p.name for p in entries) != set(expected) or any(p.is_symlink() or not p.is_file() for p in entries):
        raise ValueError("Prepared input inventory changed")
    rows = [record(p) for p in entries]
    if any(any(row[k] != expected[Path(row["path"]).name][k] for k in ("bytes", "sha256")) for row in rows):
        raise ValueError("Prepared input bytes changed")
    runtime = dict(records=[dict(path=r["absolute_path"], bytes=r["bytes"], sha256=r["sha256"])
        for r in baseline["core_sources"]])
    observed = snapshot(baseline["core_root"], dict(datasets=[dict(proteomes=run["proteomes"],
        input_directory=str(target), inputs=rows)]), runtime)
    if observed["datasets"][0]["native_order"] != run["native_order"]:
        raise ValueError("Prepared frozen enumeration differs from the prospective order")
    return observed


def prepare_inputs(run, baseline):
    start = time.monotonic_ns()
    root = Path(run["output_root"])
    original_inputs(run)
    target = root / "input"
    target.mkdir(exist_ok=False)
    by_name = {Path(r["path"]).name: r for r in run["inputs"]}
    for name in run["input_creation_order"]:
        with Path(by_name[name]["path"]).open("rb") as src, (target / name).open("xb") as dst:
            shutil.copyfileobj(src, dst)
    owners = hashlib.sha256()
    seen, counts = set(), {}
    for name in run["native_order"]:
        counts[name] = 0
        with (target / name).open() as handle:
            for line in handle:
                if line.startswith(">"):
                    gene = line[1:].split()[0]
                    if gene in seen:
                        raise ValueError("Duplicate gene ownership")
                    seen.add(gene)
                    counts[name] += 1
                    owners.update(json.dumps([name, gene], separators=(",", ":")).encode() + b"\n")
    if len(seen) != run["genes"]:
        raise ValueError("Wrong prepared gene universe")
    observed = prepared_inputs(run, baseline)
    original_inputs(run)
    result = dict(status="fresh_factorial_inputs_prepared", input_snapshot=observed,
        gene_ownership_sha256=owners.hexdigest(), genes=len(seen), per_species_counts=counts,
        started_ns=start, finished_ns=time.monotonic_ns(), inference_started=False,
        storage="persistent_filesystem_private_copy", cold_cache_claim=False)
    save(root / "preparation.json", result)
    return result


def native(plan_ref, index):
    from benchmark_tools.probe_threadripper_allocation import inspect, validate
    plan = read(plan_ref)
    runs = validate_plan(plan)
    if type(index) is not int or not 0 <= index < len(runs):
        raise ValueError("Invalid native identity")
    run = runs[index]
    baseline = read(plan["baseline"])
    placement = inspect()
    validate(placement, placement, 64)
    if (placement["affinity"] != list(range(32)) or not sys.dont_write_bytecode
            or not sys.pycache_prefix or Path(sys.pycache_prefix).exists()
            or Path(sys.pycache_prefix).is_symlink()
            or any(os.environ.get(k) != v for k, v in baseline["environment_overrides"].items())):
        raise ValueError("Native placement, bytecode or environment differs")
    if os.path.abspath(sys.executable) != baseline["tool_entrypoints"]["orthohmm_python"]["absolute_path"]:
        raise ValueError("Wrong private native interpreter")
    core = Path(baseline["core_root"])
    sys.path.insert(0, str(core))
    module = importlib.import_module("orthohmm.orthohmm")
    if Path(module.__file__).resolve() != core / "orthohmm/orthohmm.py":
        raise ValueError("Wrong frozen pipeline import")
    if module.fetch_fasta_files(str(Path(run["output_root"]) / "input")) != run["native_order"]:
        raise ValueError("Actual native enumeration changed at entry")
    adapted, flags = entrypoint(module, run["cell"])
    output = Path(run["output_root"]) / "native"
    if output.exists() or output.is_symlink():
        raise FileExistsError("Native output must be fresh")
    receipt = dict(status="native_factorial_running", index=index, cell=run["cell"], plan=plan_ref,
        factors=flags, placement=placement, native_order=run["native_order"],
        automatic_retry=False, cold_cache_claim=False, accuracy_evaluated=False)
    report = Path(run["output_root"]) / "native_execution.json"
    save(report, receipt)
    try:
        try:
            adapted(**native_kwargs(module, run, baseline))
        except SystemExit as error:
            if error.code not in (None, 0):
                raise
        metrics = json.loads((Path(run["output_root"]) / "metrics.json").read_text())
        verify_metrics(metrics, run, flags)
        receipt.update(status="native_factorial_completed_pending_output_review", stages=sorted(metrics["stages"]),
                       counts=metrics["counts"], exact_historical_output_equivalence_established=False)
    except BaseException as error:
        receipt.update(status="native_factorial_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        save(report, receipt)


def available_memory(raw):
    fields = {row.split(":", 1)[0]: row.split(":", 1)[1].split() for row in raw.splitlines() if ":" in row}
    value = fields.get("MemAvailable")
    if not value or len(value) != 2 or value[1] != "kB" or not value[0].isdigit():
        raise ValueError("Missing valid MemAvailable")
    return int(value[0]) * 1024


class ReleaseGuard:
    def __init__(self, request_ref, plan_ref, plan, run, baseline, policy_ref):
        self.request_ref, self.plan_ref, self.plan = request_ref, plan_ref, plan
        self.run, self.baseline, self.policy_ref = run, baseline, policy_ref

    def __call__(self, directory):
        from benchmark_tools.probe_host_counters import snapshot
        from benchmark_tools.observe_threadripper_process_identity import enriched_snapshot
        from benchmark_tools.review_threadripper_process_policy import review
        from benchmark_tools.slurm_resource_snapshot import scoped_path
        from benchmark_tools.verify_threadripper_controller import ReleaseBudgetGuard
        request = read(self.request_ref)
        check(self.plan_ref)
        for pin in self.plan["helper_sources"]:
            check(pin)
        prepared_inputs(self.run, self.baseline)
        ready = json.loads((directory / "ready.json").read_text())
        job = request["job_id"]
        scope = scoped_path(ready["cgroup"], job)
        job_scope = next(p for p in scope.parents if p.name == f"job_{job}")
        with (directory / "host_processes.jsonl").open() as stream:
            first = json.loads(stream.readline())
        second = enriched_snapshot()
        boot = Path("/proc/sys/kernel/random/boot_id").read_text().strip()
        policy = read(self.policy_ref)
        process_policy = read(policy["process_policy"])
        verdict = review(process_policy, first["snapshot"], second, boot_id=boot,
                         job_scope=str(job_scope), observer_pid=first["observer_pid"])
        raw_memory = Path("/proc/meminfo").read_text()
        capacity = available_memory(raw_memory)
        observation = dict(host_counters=snapshot(), raw_meminfo=raw_memory, available_memory_bytes=capacity,
            process_snapshots=[first["snapshot"], second], process_review=verdict,
            background_cpu_used_for_eligibility=False, capacity_guaranteed_through_run=False)
        save(directory / "launch_environment_observation.json", observation)
        if capacity < MEMORY or not verdict["process_policy_matched"]:
            raise ValueError("Unsafe available RAM or invalid process attribution at handoff")
        budget = ReleaseBudgetGuard(job, command=str(SCRIPT), cwd=str(ROOT), allocation_mode="shared")
        result = budget(directory)
        if result["allocation"]["fields"].get("Comment") != self.request_ref["sha256"]:
            raise ValueError("Scheduler comment differs from explicit request")
        save(directory / "environment_preflight.json", dict(schema="threadripper_environment_preflight_v1",
            decision="passed", job_id=job, index=self.run["index"], plan_sha256=self.plan_ref["sha256"],
            environment_policy=self.policy_ref, execution_scope=SCOPE, boot_id=boot,
            background_competition_recorded=True, uncontended_timing=False,
            evidence=[record(directory / "launch_environment_observation.json"), record(directory / "release_budget.json")]))
        return result


def execute(request_ref):
    from benchmark_tools.check_threadripper_runtime import RuntimeChecker
    from benchmark_tools.isolated_numba_cache import fresh_cache
    from benchmark_tools.isolated_native_tmp import fresh_tmp
    from benchmark_tools.measure_native_scaling_run import run_checked
    from benchmark_tools.measure_threadripper_run import assert_environment
    from benchmark_tools.measure_threadripper_scaling import measure
    from benchmark_tools.review_threadripper_process_stream import audit
    from benchmark_tools.run_simulation_methods import execution_environment, verify_environment
    request = read(request_ref)
    plan_ref = request["plan"]
    plan = read(plan_ref)
    runs = validate_plan(plan)
    job = int(os.environ["SLURM_JOB_ID"])
    validate_request(request, plan_ref, job)
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
    from benchmark_tools.review_threadripper_process_policy import shared_environment
    if (not shared_environment(policy) or policy["plan_sha256"] != plan_ref["sha256"]
            or policy["minimum_available_memory_bytes"] != MEMORY):
        raise ValueError("Policy does not bind this shared-host plan")
    history_evidence = []
    for i, prior_ref in enumerate(request["history"]):
        prior = read(prior_ref)
        if (prior.get("index") != i or prior.get("plan") != plan_ref
                or prior.get("terminal_reviewed") is not True or prior.get("next_identity_authorized") is not True):
            raise ValueError("Previous factorial identity unresolved")
        history_evidence.append(dict(review=prior_ref, scheduler=verify_terminal(prior["job_id"])))
    for ref in [plan_ref, request_ref, *plan["evidence"], *plan["helper_sources"]]:
        check(ref)
    run = runs[request["index"]]
    root = Path(run["output_root"])
    session = Path(plan["panel_root"]) / "sessions" / f"run_{run['index']:02d}"
    session.mkdir(parents=True, exist_ok=False)
    outcome = dict(status="factorial_executor_started", job_id=job, index=run["index"], cell=run["cell"],
        request=request_ref, plan=plan_ref, execution_scope=SCOPE, automatic_retry=False,
        scientific_timings_admitted=False, next_identity_authorized=False, uncontended_timing=False,
        contention_distortion="unknown_potentially_method_dependent")
    outcome["history_scheduler_evidence"] = history_evidence
    save(session / "started.json", outcome)
    env, _ = execution_environment(baseline)
    env.update(PYTHONDONTWRITEBYTECODE="1", PYTHONPYCACHEPREFIX=str(session / "absent_python_cache"))
    os.environ.update(env)
    os.chdir(baseline["core_root"])
    checker = RuntimeChecker(Path(plan["runtime_lookup"]["path"]), plan["runtime_lookup"]["sha256"], session / "lookup_checks")
    context = dict(cwd=baseline["core_root"])

    def check_all(specifications):
        assert_environment(context, baseline)
        for ref in [plan_ref, request_ref, request["policy"], *plan["helper_sources"], *plan["evidence"]]:
            check(ref)
        runtime = checker(specifications)
        with fresh_cache(session / f"verification_cache_{checker.count}"):
            verify_environment(baseline)
        original_inputs(run)
        copied = prepared_inputs(run, baseline) if (root / "input").exists() else None
        return dict(runtime=runtime, original_inputs=run["inputs"], prepared_inputs=copied)

    def measurement(directory):
        prepare_inputs(run, baseline)
        argv = [baseline["tool_entrypoints"]["orthohmm_python"]["absolute_path"], "-B", str(Path(__file__).resolve()),
            "--native", "--plan", plan_ref["path"], "--plan-sha256", plan_ref["sha256"], "--index", str(run["index"])]
        guard = ReleaseGuard(request_ref, plan_ref, plan, run, baseline, request["policy"])
        with fresh_tmp(root / "native_tmp"), fresh_cache(root / "native_numba_cache"):
            return measure(argv, directory, job, 32, MEMORY, 85800, 1., monitor_host=True,
                           host_interval_s=30., release_guard=guard)

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
    parser.add_argument("--request", type=Path)
    parser.add_argument("--request-sha256")
    parser.add_argument("--plan", type=Path)
    parser.add_argument("--plan-sha256")
    parser.add_argument("--index", type=int)
    args = parser.parse_args()
    if args.native:
        if args.request is not None or args.request_sha256 is not None or args.plan is None or not args.plan_sha256:
            parser.error("Native mode requires only plan, digest and index")
        ref = record(args.plan)
        if ref["sha256"] != args.plan_sha256:
            raise ValueError("Native plan digest differs")
        native(ref, args.index)
    else:
        if args.plan is not None or args.index is not None or args.plan_sha256 is not None or args.request is None or not args.request_sha256:
            parser.error("Controller mode requires only request and digest")
        ref = record(args.request)
        if ref["sha256"] != args.request_sha256:
            raise ValueError("Request digest differs")
        result = execute(ref)
        raise SystemExit(0 if result["wrapper"]["status"] == "command_exited_zero" else 1)


if __name__ == "__main__":
    main()
