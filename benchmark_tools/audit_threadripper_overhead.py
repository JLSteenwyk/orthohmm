"""Replay paired native receipts/outputs; no environment or timing admission."""

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
from statistics import median

from benchmark_tools.audit_threadripper_native_outcome import audit as audit_native
from benchmark_tools.capture_job_scheduler import terminal_record
from benchmark_tools.fingerprint_native_overhead_outputs import fingerprint
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.prepare_scaling_inputs import planned_runs
from benchmark_tools.prepare_threadripper_overhead import arm_order, relocate
from benchmark_tools.prepare_threadripper_run import paths
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.replay_threadripper_boundary import replay as replay_boundary
from benchmark_tools.replay_threadripper_scaling import replay as replay_periodic
from benchmark_tools.summarize_fixture_memory_scopes import native_cpu, peak
from benchmark_tools.verify_lineage_native_provenance import same
from benchmark_tools.verify_threadripper_controller import validate as validate_controller

IDENTITY = ("index", "pair", "arm", "method", "proteomes", "repeat")
RESOURCES = dict(host="bizon", native_workers=32, native_affinity=list(range(32)),
    scheduler_cpu_slots=64, memory_bytes=128 * 1024**3, native_timeout_s=85800)
BUDGET = dict(every_pair_max=.10, per_method_size_median_max=.05,
              required_pairs_per_method_size=3)
STATUSES = {"unrun", "scheduler_failed", "native_failure", "invalid_evidence",
            "native_outputs_replayed", "nonterminal_requires_live_check"}


def direct_pin(path):
    path = Path(path).absolute()
    if path.resolve() != path or path.is_symlink() or not path.is_file():
        raise ValueError("Require direct absolute evidence file")
    return record(path)


def read_pin(pin):
    path = Path(pin["path"])
    if not path.is_absolute() or path.resolve() != path or path.is_symlink():
        raise ValueError("Require direct absolute evidence")
    check(pin)
    value = json.loads(path.read_text())
    check(pin)
    return value


def design(plan, parent):
    expected = planned_runs()
    if (plan.get("schema") != "threadripper_native_overhead_plan_v1"
            or plan.get("status") != "prospective_native_overhead_unrun"
            or not same(plan.get("resources"), RESOURCES)
            or not same(plan.get("engineering_budget"), BUDGET)
            or not same(plan.get("parent_identities"), expected)
            or plan.get("execution_authorized") is not False
            or plan.get("scientific_timings_admitted") is not False
            or plan.get("publication_ready") is not False
            or plan.get("statistic") != "periodic_native_wall_s / boundary_native_wall_s - 1"
            or plan.get("ordering_seed") != "threadripper-native-observer-20260930-v1"
            or plan.get("failure_policy") != "retain_all_attempts_stop_after_failure_no_selective_retry"
            or not same(plan.get("common_host_observation"), dict(interval_s=30, required_in_both_arms=True))):
        raise ValueError("Frozen paired design or budget differs")
    if (len(parent["runs"]) != 27 or len(plan["runs"]) != 54
            or parent.get("status") != "threadripper_commands_prepared_unrun"
            or parent.get("scientific_execution_authorized") is not False
            or not same(parent["resources"], RESOURCES)
            or not same([{k: r[k] for k in ("index", "method", "proteomes", "repeat")}
                         for r in parent["runs"]], expected)):
        raise ValueError("Require complete frozen parent identities")
    output, inputs = Path(plan["output_root"]), Path(plan["input_root"])
    if (any(not p.is_absolute() or p.resolve() != p or p.is_symlink() for p in (output, inputs))
            or output.is_relative_to("/dev/shm") or not inputs.is_relative_to("/dev/shm")
            or inputs == Path("/dev/shm")):
        raise ValueError("Require direct persistent output and private tmpfs input roots")
    for key in ("environment_overrides", "prepend_path", "resolved_executables"):
        if not same(plan[key], parent[key]):
            raise ValueError("Paired execution environment differs from parent")
    for original in parent["runs"]:
        order = arm_order(original["method"], original["proteomes"], original["repeat"])
        for offset, arm in enumerate(order):
            index = 2 * original["index"] + offset
            expected_task = dict(index=index, pair=original["index"], arm=arm,
                method=original["method"], proteomes=original["proteomes"], repeat=original["repeat"],
                run=relocate(original, Path(plan["output_root"]), Path(plan["input_root"]), index))
            if not same(plan["runs"][index], expected_task):
                raise ValueError("Task order or full native work differs from frozen parent")


def summarize(plan, rows, issues):
    if len(rows) != 54 or len(plan["runs"]) != 54 or not same(plan["engineering_budget"], BUDGET):
        raise ValueError("Require all 54 outcomes and unchanged budgets")
    for task, row in zip(plan["runs"], rows):
        if (any(not same(row.get(k), task[k]) for k in IDENTITY)
                or row.get("status") not in STATUSES or row.get("scientific_timings_admitted") is not False):
            raise ValueError("Outcome identity, order, state or admission differs")
    pairs, ratios = [], {}
    for pair in range(27):
        selected = rows[2 * pair:2 * pair + 2]
        arms = {r["arm"]: r for r in selected}
        if (set(arms) != {"boundary", "periodic"} or any(r["pair"] != pair for r in selected)
                or any(not same(selected[0][k], selected[1][k]) for k in ("method", "proteomes", "repeat"))):
            raise ValueError("Invalid pair membership")
        result = dict(pair=pair, method=selected[0]["method"], proteomes=selected[0]["proteomes"],
            repeat=selected[0]["repeat"], order=[r["arm"] for r in selected],
            task_indices=[r["index"] for r in selected], statuses=[r["status"] for r in selected],
            status="incomplete_or_invalid_pair", signed_ratio=None, numerical_pair_budget_passed=None)
        if all(r["status"] == "native_outputs_replayed" for r in selected):
            if not same(arms["boundary"]["work_identity"], arms["periodic"]["work_identity"]):
                result["status"] = "output_mismatch"
            else:
                walls = [arms[arm]["native_wall_s"] for arm in ("boundary", "periodic")]
                times = [arms[arm]["native_duration_ns"] for arm in ("boundary", "periodic")]
                if (all(type(t) is int and 0 < t <= 85830 * 10**9 for t in times)
                        and all(type(w) in (int, float) and math.isfinite(w) and w == t / 1e9
                                for w, t in zip(walls, times))):
                    ratio = Fraction(times[1] - times[0], times[0])
                    ratios[pair] = ratio
                    result.update(status="raw_output_pair_consistent", boundary_native_wall_s=walls[0],
                        periodic_native_wall_s=walls[1], signed_ratio=float(ratio),
                        ratio_numerator_ns=times[1] - times[0], ratio_denominator_ns=times[0],
                        numerical_pair_budget_passed=ratio <= Fraction(1, 10))
        pairs.append(result)
    cells = []
    keys = sorted({(p["method"], p["proteomes"]) for p in pairs})
    if len(keys) != 9:
        raise ValueError("Require nine method-size cells")
    for method, size in keys:
        selected = [p for p in pairs if (p["method"], p["proteomes"]) == (method, size)]
        if len(selected) != 3 or sorted(p["repeat"] for p in selected) != [0, 1, 2]:
            raise ValueError("Require all three prescribed repeats")
        valid = [p for p in selected if p["status"] == "raw_output_pair_consistent"]
        value = median(ratios[p["pair"]] for p in valid) if len(valid) == 3 else None
        cells.append(dict(method=method, proteomes=size, consistent_pairs=len(valid),
            pair_indices=[p["pair"] for p in selected],
            median_signed_ratio=None if value is None else float(value),
            numerical_median_budget_passed=None if value is None else value <= Fraction(1, 20)))
    complete = not issues and all(p["status"] == "raw_output_pair_consistent" for p in pairs)
    numerical = (all(p["numerical_pair_budget_passed"] for p in pairs)
                 and all(c["numerical_median_budget_passed"] for c in cells)) if complete else None
    return dict(status="native_point_incremental_elapsed_arithmetic", pairs=pairs, cells=cells,
        raw_output_panel_complete=complete, complete_panel_numerical_budget_passed=numerical,
        engineering_budget_passed=None, panel_issues=issues,
        runtime_environment_admission_complete=False, scientific_timings_admitted=False, publication_ready=False,
        limitations=["Raw/output consistency and elapsed arithmetic do not replace full runtime/source/environment review.",
            "Only complete three-pair medians; no selective exclusion, retries or overhead subtraction.",
            "Signed ratios retain negative values; engineering thresholds are not confidence intervals.",
            "Common host-monitor cost is not isolated; no causal slowdown, timing admission or execution approval."])


def audit_task(task, baseline, job_id, launcher):
    run = task["run"]
    directory = Path(run["measurement_directory"])
    worker = str(Path(__file__).with_name("measure_threadripper_scaling.py").resolve())
    expected_launch = ["srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1",
        "--cpus-per-task=64", "--cpu-bind=mask_cpu:0xffffffff", launcher,
        "-B", worker, "--worker", str(directory)]
    if task["arm"] == "boundary":
        def replay_fn(path, job, command):
            return replay_boundary(path, job, command, expected_launcher=launcher, expected_worker=worker)
    elif task["arm"] == "periodic":
        replay_fn = replay_periodic
    else:
        raise ValueError("Unknown collector arm")
    native = audit_native(run, baseline, job_id, replay_fn=replay_fn)
    raw = native["replay"]
    measured = raw["measured"]
    if (not same(measured["launched"], expected_launch)
            or task["arm"] == "periodic" and measured["schema"] != "threadripper_scaling_v5"):
        raise ValueError("Native worker launch/schema differs from external bindings")
    evidence = [*native["evidence"], *raw["evidence"]]
    result = {k: task[k] for k in IDENTITY}
    result.update(job_id=job_id, status="native_failure", native_outcome=native["native_outcome"],
        native_wall_s=raw["native_wall_s"],
        native_duration_ns=measured["native"]["finished_ns"] - measured["native"]["started_ns"],
        cpu=native_cpu(measured["native"], job_id, raw["native_completion"]),
        step_memory=peak(measured["step_memory"]), scientific_timings_admitted=False,
        runtime_environment_admission_complete=False)
    if native["outputs"] is not None:
        output = Path(run["configuration"]["output"])
        if output.resolve() != output or output.is_symlink() or any(p.is_symlink() for p in output.rglob("*")):
            raise ValueError("Indirect native outputs require separate review")
        identity = fingerprint(run, {"/": "/"})
        result.update(status="native_outputs_replayed", work_identity=identity["identity"],
            output_checks=native["outputs"], output_fingerprint=identity)
        evidence.extend(identity["evidence"])
        evidence.extend([identity["source"], *identity["helpers"]])
        evidence.extend(native["outputs"]["native"]["checked_files"])
    for item in evidence:
        check(item)
    result["evidence"] = evidence
    return result


def audit(plan_path, digest, attempts_path, baseline_ref, launcher, scheduler_command, scheduler_cwd):
    plan_ref, attempts_ref = direct_pin(plan_path), direct_pin(attempts_path)
    if plan_ref["sha256"] != digest:
        raise ValueError("Overhead plan digest differs")
    plan, attempts = read_pin(plan_ref), read_pin(attempts_ref)
    parent = read_pin(plan["sources"][0])
    design(plan, parent)
    for pin in [*plan["sources"], *plan["helpers"]]:
        check(pin)
    if (attempts.get("schema") != "threadripper_native_overhead_attempts_v1"
            or not same(attempts.get("plan"), plan_ref)
            or attempts.get("automatic_retry") is not False or not isinstance(attempts.get("attempts"), list)
            or len(attempts["attempts"]) > 54):
        raise ValueError("Require bound one-attempt prefix")
    for path in (launcher, scheduler_command, scheduler_cwd):
        if path is None and not attempts["attempts"]:
            continue
        if not isinstance(path, str) or not Path(path).is_absolute() or ".." in Path(path).parts:
            raise ValueError("Require external absolute launch and scheduler bindings")
    jobs = set()
    for index, attempt in enumerate(attempts["attempts"]):
        if (type(attempt.get("index")) is not int or attempt["index"] != index
                or type(attempt.get("job_id")) is not int or attempt["job_id"] <= 0 or attempt["job_id"] in jobs):
            raise ValueError("Attempts skip, reorder, duplicate or retry an identity")
        jobs.add(attempt["job_id"])
    if attempts["attempts"] and baseline_ref is None:
        raise ValueError("Attempted work requires an externally bound baseline")
    baseline = read_pin(baseline_ref) if attempts["attempts"] else None
    rows, issues, unrun_paths, previous_finish, previous_boot = [], [], [], None, None
    failed = None
    evidence = [plan_ref, attempts_ref, *plan["sources"], *plan["helpers"]]
    if baseline is not None:
        evidence.append(baseline_ref)
    for task in plan["runs"]:
        row = {k: task[k] for k in IDENTITY}
        row.update(status="unrun", scientific_timings_admitted=False)
        try:
            index = task["index"]
            if index >= len(attempts["attempts"]):
                root, inputs = paths(task["run"])
                if any(p.exists() or p.is_symlink() for p in (root, inputs.parent)):
                    raise ValueError("Artifacts exist for a task absent from attempt history")
                unrun_paths.extend((root, inputs.parent))
            else:
                attempt = attempts["attempts"][index]
                scheduler_pin = attempt["scheduler"]
                check(scheduler_pin)
                scheduler = Path(scheduler_pin["path"])
                if scheduler.resolve() != scheduler or scheduler.is_symlink():
                    raise ValueError("Indirect scheduler evidence")
                evidence.append(scheduler_pin)
                if failed is not None:
                    raise ValueError("Later attempt follows an invalid predecessor")
                raw_scheduler = scheduler.read_text()
                if terminal_record(raw_scheduler, attempt["job_id"]) is None:
                    row.update(status="nonterminal_requires_live_check", job_id=attempt["job_id"], scheduler=scheduler_pin)
                else:
                    allocation = validate_controller(raw_scheduler, attempt["job_id"], "terminal",
                        command=scheduler_command, cwd=scheduler_cwd, time_limit="1-02:00:00")
                    if allocation["scheduler_state"] != "COMPLETED" or allocation["scheduler_exit_code"] != "0:0":
                        row.update(status="scheduler_failed", job_id=attempt["job_id"], scheduler=allocation)
                    else:
                        row = audit_task(task, baseline, attempt["job_id"], launcher)
                        row["scheduler"] = allocation
                        evidence.extend(row["evidence"])
                        done = read_pin(record(Path(task["run"]["measurement_directory"]) / "done.json"))
                        start, finish = done["snapshots"][0]["started_monotonic_ns"], done["snapshots"][1]["finished_monotonic_ns"]
                        boot = done["snapshots"][0]["raw"]["boot_id"]
                        if previous_finish is not None and (start <= previous_finish or boot != previous_boot):
                            raise ValueError("Native observation brackets overlap or boot identity changed")
                        previous_finish, previous_boot = finish, boot
                        if index % 2 and rows[-1]["status"] == row["status"] == "native_outputs_replayed":
                            if not same(rows[-1]["work_identity"], row["work_identity"]):
                                failed = index
                                issues.append(dict(index=index, stage="paired_output", reason="Canonical outputs differ"))
                if row["status"] != "native_outputs_replayed":
                    failed = index
        except (OSError, ValueError, KeyError, TypeError, IndexError) as error:
            row.update(status="invalid_evidence", reason=str(error), error_type=type(error).__name__)
            if failed is None:
                failed = task["index"]
            issues.append(dict(index=task["index"], stage="raw_output_task", error_type=type(error).__name__, reason=str(error)))
        rows.append(row)
    summary = summarize(plan, rows, issues)
    for pin in evidence:
        check(pin)
    if any(p.exists() or p.is_symlink() for p in unrun_paths):
        raise ValueError("Unrun artifact inventory changed during audit")
    return dict(status="threadripper_overhead_raw_output_audit_pending_full_admission", plan=plan_ref,
        attempt_history=attempts_ref, runs=rows, comparison=summary, evidence=evidence,
        external_bindings=dict(launcher=launcher, scheduler_command=scheduler_command, scheduler_cwd=scheduler_cwd),
        sources=[record(__file__)], scientific_timings_admitted=False, publication_ready=False,
        limitations=summary["limitations"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("plan", "attempts", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--plan-sha256", required=True)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--baseline-sha256")
    for name in ("launcher", "scheduler-command", "scheduler-cwd"):
        parser.add_argument("--" + name)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    baseline = None if args.baseline is None else direct_pin(args.baseline)
    if baseline is not None and baseline["sha256"] != args.baseline_sha256:
        parser.error("Baseline digest differs")
    result = audit(args.plan, args.plan_sha256, args.attempts, baseline,
                   args.launcher, args.scheduler_command, args.scheduler_cwd)
    save(args.output, result)
    print(result["status"])
