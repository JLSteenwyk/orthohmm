"""Audit all native root-context outcomes without admitting scientific timing."""

import argparse
import gzip
import json
from pathlib import Path
import re

from benchmark_tools.audit_frontier_overhead import recipe_evidence
from benchmark_tools.audit_lineage_native_diagnostics import inventory
from benchmark_tools.capture_job_scheduler import terminal_record
from benchmark_tools.fingerprint_native_overhead_outputs import fingerprint
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_root_context_native import PROTOCOL_SHA
from benchmark_tools.replay_native_root_context import replay
from benchmark_tools.validate_scaling_outputs import validate
from benchmark_tools.validate_simulation_outputs import NativeOutputFailure
from benchmark_tools.verify_root_context_native_provenance import (
    ROOT, RECIPE_ROOT, PLAN_SHA, load_context, same, scheduler_identity, verify)

PRIOR_SHA = "0e3a77f785b645f44e0ca5cb8bbaeffbdb7106e2104350045e299c658de9050a"
OUTPUT = "root_context_native_v1"
ERRORS = (OSError, ValueError, KeyError, TypeError, IndexError, NativeOutputFailure)


def read(path):
    return json.loads(path.read_text())


def bind_panel(directory, context, job, scheduler):
    panel, launch = read(directory / "result.json"), read(directory / "launch.json")
    expected = dict(job_id=job, plan_sha256=PLAN_SHA, recipe_sha256=context["recipe_sha256"],
        scientific_timings_admitted=False, native_outputs_validated=False, publication_ready=False)
    if any(not same(panel[k], v) for k, v in expected.items()) or len(panel["tasks"]) != 3:
        raise ValueError("Panel identity or task inventory differs")
    expected_launch = dict(job_id=job, plan_sha256=PLAN_SHA, recipe_sha256=context["recipe_sha256"],
        protocol_sha256=PROTOCOL_SHA, executable=context["plan"]["launcher_python"],
        source_directory=str(RECIPE_ROOT), scientific_timings_admitted=False)
    if any(not same(launch[k], v) for k, v in expected_launch.items()):
        raise ValueError("Launch identity differs")
    affinity = launch["affinity"]
    if (len(affinity) != 20 or any(type(cpu) is not int or cpu < 0 for cpu in affinity)
            or affinity != sorted(set(affinity)) or launch["uname"][1] != "spark-7ff0"):
        raise ValueError("Launch host or affinity differs")
    failed = None
    for index, (row, task) in enumerate(zip(panel["tasks"], context["plan"]["runs"])):
        expected = dict(index=index, task=task, job_id=job, scientific_timings_admitted=False)
        if any(not same(row[k], v) for k, v in expected.items()):
            raise ValueError("Task order or identity differs")
        if failed is not None:
            expected.update(status="not_run_after_failure", failed_index=failed)
            if not same(row, expected):
                raise ValueError("Work launched after failure")
            path = directory.parent / Path(task["run"]["measurement_directory"]).parent.relative_to(ROOT)
            if path.exists():
                raise ValueError("Unrun task has native artifacts")
        elif row["status"] == "failed":
            if not isinstance(row["error"], str) or not isinstance(row["error_type"], str):
                raise ValueError("Failure evidence missing")
            failed = index
        elif row["status"] != "measurement_completed":
            raise ValueError("Unknown task status")
        if not same(read(directory / f"task_{index:02d}.json"), row):
            raise ValueError("Task receipt differs from panel")
        expected = dict(tasks=panel["tasks"][:index+1], failed_index=failed, job_id=job)
        if not same(read(directory / f"progress_{index:02d}.json"), expected):
            raise ValueError("Task checkpoint differs")
    expected_status = "root_context_native_completed" if failed is None else "root_context_native_stopped"
    fields = dict(re.findall(r"(?<!\S)([A-Za-z][^\s=]*)=([^\s]+)", scheduler))
    state, exit_code = ("COMPLETED", "0:0") if failed is None else ("FAILED", "1:0")
    if panel["status"] != expected_status or fields["JobState"] != state or fields["ExitCode"] != exit_code:
        raise ValueError("Panel and scheduler outcomes differ")
    allowed = {"result.json", "launch.json"} | {f"{prefix}_{i:02d}.json" for i in range(3)
        for prefix in ("task", "progress")}
    allowed.update(Path(task["run"]["measurement_directory"]).parent.name
                   for task, row in zip(context["plan"]["runs"], panel["tasks"])
                   if row["status"] != "not_run_after_failure")
    if any(path.name not in allowed for path in directory.iterdir()):
        raise ValueError("Unexpected panel artifact")
    return panel


def successful_task(archive, context, task, receipt, scheduler, job, prior):
    directory = archive / Path(task["run"]["measurement_directory"]).parent.relative_to(ROOT)
    measured = read(directory / "measurement/lineage_report.json")
    binding = verify(context, task["index"], read(directory / "preparation.json"),
        directory / "verification.json", receipt, measured, scheduler, job, require_completed=False)
    replayed = replay(directory / "measurement", job, binding["measured_argv"], expected_timeout_s=900)
    lineage = replayed["lineage"]
    native = measured["native"]
    adapted = dict(command=binding["measured_argv"], cwd=task["run"]["cwd"],
                   exit_code=native["exit_code"], timed_out=native["timed_out"])
    checked = validate(binding["run"], adapted, {ROOT: archive})
    canonical = fingerprint(binding["run"], {ROOT: archive})
    previous = [row for row in prior["runs"] if row["index"] == task["index"]]
    if len(previous) != 1 or previous[0]["status"] != "validated" or previous[0]["method"] != task["method"]:
        raise ValueError("Missing validated same-method lineage reference")
    previous = previous[0]
    equivalent = same(canonical["identity"], previous["work_identity"])
    evidence = [*checked["checked_files"], *canonical["evidence"], checked["gnu_time_companion"]["source"],
                *replayed["evidence"], *previous["inventory"], *previous["checked_evidence"]]
    for item in evidence:
        check(item)
    points = lineage["measured"]["points"]
    first, last = points[0]["root_context"], points[-1]["root_context"]
    return dict(index=task["index"], method=task["method"], job_id=job,
        status="validated" if equivalent else "output_mismatch", output_equivalent=equivalent,
        work_identity=canonical["identity"], prior_work_identity=previous["work_identity"],
        native_wall_s=lineage["native_wall_s"], memory=lineage["memory"], screening=lineage["screening"],
        original_flagged_intervals=lineage["original_flagged_intervals"],
        narrow_flagged_intervals=lineage["narrow_flagged_intervals"], root_context=replayed["context"],
        gnu_time=checked["gnu_time_companion"]["accounting"],
        native_counts={k: checked[k] for k in ("input_genes", "orthogroups", "root_hogs", "checkpoint_groups", "native_pair_rows") if k in checked},
        observation_start_ns=points[0]["host"][0]["started_monotonic_ns"],
        observation_finish_ns=last["host_after"]["finished_ns"],
        scope_identity={k: first[k] for k in ("boot_before", "ticks", "identities_before")},
        evidence=evidence, scientific_timings_admitted=False)


def audit(archive, results, recipe_path, recipe_sha, scheduler_path, job, prior_path):
    scheduler_source = record(scheduler_path)
    scheduler = scheduler_path.read_text()
    if terminal_record(scheduler, job) is None:
        raise ValueError("Require terminal allocation before archive inspection")
    scheduler_identity(scheduler, job, require_completed=False)
    sources = [scheduler_source, record(recipe_path), record(prior_path)]
    if sources[-1]["sha256"] != PRIOR_SHA:
        raise ValueError("Prior lineage audit hash differs")
    prior = json.loads(gzip.decompress(prior_path.read_bytes()))
    context = load_context(results, recipe_path, recipe_sha)
    frozen = [record(results / name) for name in ("dgx_root_context_native_plan_20260919.json",
        "dgx_lineage_native_plan_20260919.json", "ROOT_CONTEXT_NATIVE_PROTOCOL_20260919.md")]
    for source in frozen:
        target = str(RECIPE_ROOT / "benchmark_tools/results" / Path(source["path"]).name)
        entries = [r for r in context["recipe"]["records"] if r["path"] == target and r["kind"] == "file"]
        if len(entries) != 1 or entries[0]["sha256"] != source["sha256"]:
            raise ValueError("Archived recipe does not pin the frozen plan inputs")
    sources.extend(frozen)
    sources.extend(recipe_evidence(archive, context["recipe"], RECIPE_ROOT.name))
    sources.extend(record(Path(__file__).with_name(name)) for name in (
        "verify_root_context_native_provenance.py", "replay_native_root_context.py",
        "replay_lineage_native_measurement.py", "validate_scaling_outputs.py",
        "fingerprint_native_overhead_outputs.py", "audit_frontier_overhead.py"))
    directory = archive / OUTPUT
    evidence = inventory(directory)
    issues, rows = [], []
    try:
        panel = bind_panel(directory, context, job, scheduler)
    except ERRORS as error:
        panel = None
        issues.append(dict(stage="panel_binding", error_type=type(error).__name__, reason=str(error)))
    for task in context["plan"]["runs"]:
        row = dict(index=task["index"], method=task["method"], job_id=job, scientific_timings_admitted=False)
        try:
            if panel is None:
                raise ValueError("Panel incomplete or invalid; no task admitted")
            receipt = panel["tasks"][task["index"]]
            if receipt["status"] == "measurement_completed":
                row = successful_task(archive, context, task, receipt, scheduler, job, prior)
            else:
                row.update(status="retained_unrun" if receipt["status"] == "not_run_after_failure" else "retained_failure",
                           retained=receipt)
        except ERRORS as error:
            row.update(status="failed_or_invalid", error_type=type(error).__name__, reason=str(error))
        rows.append(row)
    observed = [r for r in rows if r["status"] in {"validated", "output_mismatch"}]
    for a, b in zip(observed, observed[1:]):
        if not same(a["scope_identity"], b["scope_identity"]) or a["observation_finish_ns"] >= b["observation_start_ns"]:
            issues.append(dict(stage="temporal_identity", left=a["index"], right=b["index"]))
    for item in [*sources, *evidence]:
        check(item)
    if evidence != inventory(directory):
        raise ValueError("Panel inventory changed during audit")
    return dict(status="root_context_native_audited_not_scientific_admission", runs=rows, issues=issues,
        validated_tasks=sum(r["status"] == "validated" for r in rows),
        all_outputs_equivalent=(all(r["output_equivalent"] for r in observed) if len(observed) == 3 else None),
        all_tasks_validated=not issues and all(r["status"] == "validated" for r in rows),
        sources=sources, inventory=evidence, source=record(__file__),
        scientific_timings_admitted=False, environmental_validity_established=False, publication_ready=False,
        limitations=["Failed, unrun and invalid outcomes are retained; missing outcomes are not successes.",
            "Held-session receipt binding requires a separate audit.",
            "Output equivalence does not establish accuracy, identical internal work or collector overhead.",
            "Root-context observations do not establish causal attribution or absence of interference."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("archive", "results", "recipe", "scheduler", "prior-audit", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--recipe-sha", required=True)
    parser.add_argument("--job", type=int, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.archive.resolve(), args.results.resolve(), args.recipe.resolve(), args.recipe_sha,
                   args.scheduler.resolve(), args.job, args.prior_audit.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
