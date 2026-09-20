"""Audit the complete paired collector panel without scientific timing admission."""

import argparse
import gzip
import json
from pathlib import Path

from benchmark_tools.audit_frontier_overhead import recipe_evidence
from benchmark_tools.audit_lineage_native_diagnostics import inventory
from benchmark_tools.audit_root_context_overhead_session import audit as audit_session
from benchmark_tools.fingerprint_native_overhead_outputs import fingerprint
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_root_context_overhead import PROTOCOL_SHA, PRIOR_AUDIT_SHA
from benchmark_tools.replay_lineage_native_measurement import replay as replay_lineage
from benchmark_tools.replay_native_root_context import replay as replay_root
from benchmark_tools.summarize_root_context_overhead import summarize
from benchmark_tools.validate_scaling_outputs import validate
from benchmark_tools.validate_simulation_outputs import NativeOutputFailure
from benchmark_tools.verify_root_context_overhead_provenance import (
    ROOT, RECIPE_ROOT, PLAN_SHA, load_context, same, scheduler_identity, verify)

OUTPUT = "root_context_overhead_v1"
IDENTITY = ("index", "block", "pair", "arm", "method")
ERRORS = (OSError, ValueError, KeyError, TypeError, IndexError, NativeOutputFailure)


def read(path):
    return json.loads(path.read_text())


def bind_panel(directory, context, job, allocation):
    panel, launch = read(directory / "result.json"), read(directory / "launch.json")
    expected = dict(job_id=job, plan_sha256=PLAN_SHA, recipe_sha256=context["recipe_sha256"],
        scientific_timings_admitted=False, native_outputs_validated=False, publication_ready=False)
    if any(not same(panel[k], v) for k, v in expected.items()) or len(panel["tasks"]) != 18:
        raise ValueError("Panel identity or task inventory differs")
    expected = dict(job_id=job, plan_sha256=PLAN_SHA, recipe_sha256=context["recipe_sha256"],
        protocol_sha256=PROTOCOL_SHA, executable=context["plan"]["launcher_python"],
        source_directory=str(RECIPE_ROOT), scientific_timings_admitted=False)
    if any(not same(launch[k], v) for k, v in expected.items()):
        raise ValueError("Launch identity differs")
    affinity = launch["affinity"]
    if (len(affinity) != 20 or any(type(cpu) is not int or cpu < 0 for cpu in affinity)
            or affinity != sorted(set(affinity)) or launch["uname"][1] != "spark-7ff0"):
        raise ValueError("Launch host or affinity differs")
    allowed = {"launch.json", "result.json"}
    failed = None
    for index, (task, row) in enumerate(zip(context["plan"]["runs"], panel["tasks"])):
        expected = {key: task[key] for key in IDENTITY}
        expected.update(task=task, job_id=job, scientific_timings_admitted=False)
        if any(not same(row[k], v) for k, v in expected.items()):
            raise ValueError("Task, arm or pair identity differs")
        run_path = directory.parent / Path(task["run"]["measurement_directory"]).parent.relative_to(ROOT)
        if failed is not None:
            expected.update(status="not_run_after_failure", failed_index=failed)
            if not same(row, expected) or run_path.exists() or run_path.is_symlink():
                raise ValueError("Work or artifacts after failure")
        else:
            allowed.add(run_path.name)
            if row["status"] == "failed":
                if not isinstance(row["error"], str) or not isinstance(row["error_type"], str):
                    raise ValueError("Missing failure evidence")
                failed = index
            elif row["status"] != "measurement_completed":
                raise ValueError("Unknown task status")
        task_name, checkpoint_name = f"task_{index:02d}.json", f"progress_{index:02d}.json"
        allowed.update((task_name, checkpoint_name))
        if not same(read(directory / task_name), row):
            raise ValueError("Task receipt differs from panel")
        checkpoint = dict(tasks=panel["tasks"][:index+1], failed_index=failed, job_id=job)
        if not same(read(directory / checkpoint_name), checkpoint):
            raise ValueError("Cumulative checkpoint differs")
    expected_status = "root_context_overhead_completed" if failed is None else "root_context_overhead_stopped"
    state, exit_code = ("COMPLETED", "0:0") if failed is None else ("FAILED", "1:0")
    if panel["status"] != expected_status or allocation["JobState"] != state or allocation["ExitCode"] != exit_code:
        raise ValueError("Panel and terminal scheduler outcomes differ")
    if any(path.name not in allowed for path in directory.iterdir()):
        raise ValueError("Unexpected panel artifact")
    return panel


def successful_task(archive, context, task, receipt, scheduler, job, prior):
    directory = archive / Path(task["run"]["measurement_directory"]).parent.relative_to(ROOT)
    binding = verify(context, task["index"], directory, receipt, scheduler, job, require_completed=False)
    if task["arm"] == "root_context":
        replayed = replay_root(directory / "measurement", job, binding["measured_argv"], expected_timeout_s=900)
        lineage, root_context = replayed["lineage"], replayed["context"]
    else:
        replayed = replay_lineage(directory / "measurement", job, binding["measured_argv"], expected_timeout_s=900)
        lineage, root_context = replayed, None
    measured = binding["measurement"]
    native = measured["native"]
    checked = validate(binding["run"], dict(command=binding["measured_argv"], cwd=task["run"]["cwd"],
        exit_code=native["exit_code"], timed_out=native["timed_out"]), {ROOT: archive})
    canonical = fingerprint(binding["run"], {ROOT: archive})
    previous = [r for r in prior["runs"] if r["index"] == task["native_parent_index"]]
    if len(previous) != 1 or previous[0]["status"] != "validated" or previous[0]["method"] != task["method"]:
        raise ValueError("Missing validated native baseline")
    equivalent = same(canonical["identity"], previous[0]["work_identity"])
    evidence = [*binding["evidence"], *replayed["evidence"], *checked["checked_files"],
                *canonical["evidence"], checked["gnu_time_companion"]["source"]]
    for item in evidence:
        check(item)
    points = lineage["measured"]["points"]
    first, last = points[0], points[-1]
    common_identity = dict(boot=first["lineage"]["boot_before"], ticks=first["ticks"],
        scopes={scope: first["lineage"]["identities_before"][scope] for scope in ("/", "/system.slice")})
    finish = last["host"][1]["finished_monotonic_ns"]
    root_identity = None
    if task["arm"] == "root_context":
        root_identity = first["root_context"]["identities_before"]
        finish = last["root_context"]["host_after"]["finished_ns"]
    result = {key: task[key] for key in IDENTITY}
    result.update(job_id=job, status="validated" if equivalent else "output_mismatch", output_equivalent=equivalent,
        work_identity=canonical["identity"], prior_work_identity=previous[0]["work_identity"],
        native_wall_s=lineage["native_wall_s"], memory=lineage["memory"],
        original_flagged_intervals=lineage["original_flagged_intervals"],
        narrow_flagged_intervals=lineage["narrow_flagged_intervals"],
        observation_intervals=len(points)-1, root_context=root_context,
        pressure_observation_window=lineage["screening"]["observation_window"]["native_pressure"],
        gnu_time=checked["gnu_time_companion"]["accounting"],
        native_counts={k: checked[k] for k in ("input_genes", "orthogroups", "root_hogs", "checkpoint_groups", "native_pair_rows") if k in checked},
        observation_start_ns=first["host"][0]["started_monotonic_ns"],
        observation_finish_ns=max(finish, lineage["memory"]["finished_ns"]),
        common_scope_identity=common_identity, root_scope_identity=root_identity,
        evidence=evidence, scientific_timings_admitted=False)
    return result


def temporal_issues(rows):
    observed = [r for r in rows if r["status"] in {"validated", "output_mismatch"}]
    issues = []
    for a, b in zip(observed, observed[1:]):
        if (not same(a["common_scope_identity"], b["common_scope_identity"])
                or a["observation_finish_ns"] >= b["observation_start_ns"]):
            issues.append(dict(stage="temporal_identity", left=a["index"], right=b["index"]))
    roots = [r for r in observed if r["arm"] == "root_context"]
    for a, b in zip(roots, roots[1:]):
        if not same(a["root_scope_identity"], b["root_scope_identity"]):
            issues.append(dict(stage="root_scope_identity", left=a["index"], right=b["index"]))
    return issues


def audit(archive, results, recipe_path, recipe_sha, scheduler_path, job, prior_path,
          session_directory, scheduler_timezone):
    sources = [record(scheduler_path)]
    scheduler = scheduler_path.read_text()
    allocation = scheduler_identity(scheduler, job, require_completed=False)
    issues = []
    try:
        session = audit_session(session_directory, scheduler_path, recipe_path, recipe_sha, job, scheduler_timezone)
        sources.extend(session["evidence"])
    except ERRORS as error:
        session = dict(status="failed_or_invalid", error_type=type(error).__name__, reason=str(error))
        issues.append(dict(stage="waiting_session", **session))
    sources.extend((record(recipe_path), record(prior_path)))
    if sources[-1]["sha256"] != PRIOR_AUDIT_SHA:
        raise ValueError("Prior native audit hash differs")
    prior = json.loads(gzip.decompress(prior_path.read_bytes()))
    if prior["all_tasks_validated"] is not True or prior["issues"] or len(prior["runs"]) != 3:
        raise ValueError("Prior native baseline is incomplete or invalid")
    baseline_evidence = prior["inventory"]
    for item in baseline_evidence:
        check(item)
    context = load_context(results, recipe_path, recipe_sha)
    frozen = [record(results / name) for name in ("dgx_root_context_overhead_plan_20260919.json",
        "dgx_root_context_native_plan_20260919.json", "ROOT_CONTEXT_OVERHEAD_PROTOCOL_20260919.md")]
    for item in frozen:
        target = str(RECIPE_ROOT / "benchmark_tools/results" / Path(item["path"]).name)
        entries = [r for r in context["recipe"]["records"] if r["kind"] == "file" and r["path"] == target]
        if len(entries) != 1 or entries[0]["sha256"] != item["sha256"]:
            raise ValueError("Recipe does not pin frozen plan inputs")
    sources.extend(frozen)
    sources.extend(recipe_evidence(archive, context["recipe"], RECIPE_ROOT.name))
    sources.extend(record(Path(__file__).with_name(name)) for name in (
        "verify_root_context_overhead_provenance.py", "replay_native_root_context.py", "replay_lineage_native_measurement.py",
        "validate_scaling_outputs.py", "fingerprint_native_overhead_outputs.py", "summarize_root_context_overhead.py",
        "audit_root_context_overhead_session.py", "submit_root_context_overhead_session.py"))
    directory = archive / OUTPUT
    evidence = inventory(directory)
    rows = []
    try:
        panel = bind_panel(directory, context, job, allocation)
    except ERRORS as error:
        panel = None
        issues.append(dict(stage="panel_binding", error_type=type(error).__name__, reason=str(error)))
    for task in context["plan"]["runs"]:
        row = {key: task[key] for key in IDENTITY}
        row.update(job_id=job, scientific_timings_admitted=False)
        try:
            if panel is None:
                raise ValueError("Panel incomplete or invalid")
            receipt = panel["tasks"][task["index"]]
            if receipt["status"] == "measurement_completed":
                row = successful_task(archive, context, task, receipt, scheduler, job, prior)
            else:
                row.update(status="retained_unrun" if receipt["status"] == "not_run_after_failure" else "retained_failure",
                           retained=receipt)
        except ERRORS as error:
            row.update(status="failed_or_invalid", error_type=type(error).__name__, reason=str(error))
        rows.append(row)
    issues.extend(temporal_issues(rows))
    comparison = summarize(context["plan"], rows, issues)
    for item in [*sources, *evidence, *baseline_evidence]:
        check(item)
    if evidence != inventory(directory):
        raise ValueError("Panel evidence changed during audit")
    return dict(status="root_context_overhead_audited_not_scientific_admission", runs=rows, issues=issues,
        validated_tasks=sum(r["status"] == "validated" for r in rows), comparison=comparison,
        waiting_session=session,
        sources=sources, inventory=evidence, baseline_evidence=baseline_evidence, source=record(__file__),
        scientific_timings_admitted=False, environmental_validity_established=False, publication_ready=False,
        limitations=["All prescribed outcomes and flags retained; raw interval pressure/screens remain in hashed reports.",
            "Missing/invalid waiting receipts prevent a panel-wide engineering-budget conclusion.",
            "Engineering elapsed-time budgets do not establish causal overhead, isolation or timing eligibility."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("archive", "results", "recipe", "scheduler", "prior-audit", "session-directory", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--recipe-sha", required=True)
    parser.add_argument("--job", type=int, required=True)
    parser.add_argument("--scheduler-timezone", required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.archive.resolve(), args.results.resolve(), args.recipe.resolve(), args.recipe_sha,
                   args.scheduler.resolve(), args.job, args.prior_audit.resolve(),
                   args.session_directory.resolve(), args.scheduler_timezone)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
