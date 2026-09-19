"""Audit a complete terminal native overhead panel without scientific admission."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.verify_frontier_overhead_provenance import load_context, verify, ROOT, panel_spec, PANELS
from benchmark_tools.replay_frontier_overhead_measurement import replay
from benchmark_tools.fingerprint_native_overhead_outputs import fingerprint
from benchmark_tools.validate_scaling_outputs import validate
from benchmark_tools.validate_simulation_outputs import NativeOutputFailure
from benchmark_tools.summarize_frontier_overhead import terminal_scheduler_rows, summarize
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


def recipe_evidence(archive, recipe, recipe_root="frontier_overhead_recipe_v1"):
    root = archive / recipe_root
    if root.is_symlink():
        raise ValueError("Symlinked archive recipe root")
    expected = {}
    for row in recipe["records"]:
        if row["kind"] == "file":
            path = archive / Path(row["path"]).relative_to(ROOT)
            expected[path] = row
    observed = []
    for path in root.rglob("*"):
        if path.is_symlink():
            raise ValueError("Symlink in archived recipe")
        if path.is_file():
            observed.append(path)
    if set(observed) != set(expected):
        raise ValueError("Archived recipe inventory differs")
    evidence = []
    for path, original in expected.items():
        item = record(path)
        if (item["bytes"], item["sha256"]) != (original["bytes"], original["sha256"]):
            raise ValueError("Archived recipe source differs: " + str(path))
        evidence.append(item)
    return evidence


def successful_task(archive, context, index, scheduler):
    task = context["plan"]["runs"][index]
    directory = archive / Path(task["run"]["measurement_directory"]).parent.relative_to(ROOT)
    report = "boundary_report.json" if task["mode"] == "boundary" else "frontier_report.json"
    paths = [directory / name for name in ("preparation.json", "verification.json", "overhead_task.json")]
    paths += [directory / "measurement" / report, archive / f"scheduler_{index}.txt"]
    evidence = [record(path) for path in paths]
    preparation, verified, receipt, measured = [json.loads(path.read_text()) for path in paths[:4]]
    binding = verify(context, index, preparation, verified, receipt, measured, paths[4].read_text())
    if (scheduler["allocated_cpus"] != "20" or scheduler["requested_memory"] not in {"96G", "96Gn"}
            or scheduler["node"] != "spark-7ff0"):
        raise ValueError("Accounting resource allocation differs")
    pressure = panel_spec(context.get("panel", "frontier_21838"))["pressure"]
    options = {"expected_native_pressure": True} if pressure else {}
    replayed = replay(directory / "measurement", task["mode"], binding["job_id"], binding["measured_argv"], **options)
    native = measured["native"]
    adapted = dict(command=binding["measured_argv"], cwd=task["run"]["cwd"],
                   exit_code=native["exit_code"], timed_out=native["timed_out"])
    checked = validate(binding["run"], adapted, {ROOT: archive})
    canonical = fingerprint(binding["run"], {ROOT: archive})
    evidence.extend(replayed["evidence"])
    evidence.extend(checked["checked_files"])
    evidence.extend(canonical["evidence"])
    evidence.append(checked["gnu_time_companion"]["source"])
    for item in evidence:
        check(item)
    result = dict(index=index, method=task["method"], mode=task["mode"], pair=task["pair"], status="validated",
        scheduler=scheduler, native_wall_s=replayed["native_wall_s"], work_identity=canonical["identity"],
        whole_command_screen_passed=replayed["whole_command_screen_passed"],
        flagged_intervals=replayed["flagged_intervals"], evidence=evidence,
        native_counts={key: checked[key] for key in ("input_genes", "orthogroups", "root_hogs", "checkpoint_groups", "native_pair_rows") if key in checked},
        gnu_time=checked["gnu_time_companion"]["accounting"], memory=replayed["memory"],
        observation_points=len(measured["points"]), started_ns=native["started_ns"], finished_ns=native["finished_ns"],
        boot_id=measured["points"][0]["frontier"]["boot_id"], scientific_timings_admitted=False)
    if pressure:
        result["native_pressure"] = dict(
            whole_command=replayed["screening"]["native_pressure_whole_command"],
            intervals=replayed["screening"].get("native_pressure_intervals"),
            diagnostic_only=True, environmental_validity_established=False)
    return result


def failure_evidence(archive, task):
    directory = archive / Path(task["run"]["measurement_directory"]).parent.relative_to(ROOT)
    paths = [directory / name for name in ("preparation.json", "verification.json", "overhead_task.json", "native.time.tsv")]
    paths += [directory / "measurement" / name for name in (
        "command.json", "done.json", "frontier_report.json", "boundary_report.json", "step_memory.json", "native.log", "step.log")]
    paths.append(archive / f"scheduler_{task['index']}.txt")
    return [record(path) for path in paths if path.is_file()]


def audit(archive, results, accounting_path, *, panel="frontier_21838"):
    spec = panel_spec(panel)
    accounting_record = record(accounting_path)
    # This gate precedes even loading native archive metadata.
    scheduler = terminal_scheduler_rows(accounting_path.read_text(), spec["array_id"])
    sources = [record(results / name) for name in (
        spec["plan_file"], spec["recipe_file"], spec["auth_file"], *spec["protocols"])]
    context = load_context(results, panel)
    for source in sources[3:]:
        protocols = [row for row in context["recipe"]["records"]
                     if row["path"].endswith("/" + Path(source["path"]).name)]
        if len(protocols) != 1 or source["sha256"] != protocols[0]["sha256"]:
            raise ValueError("Prospective protocol differs from frozen recipe")
    recipes = recipe_evidence(archive, context["recipe"], spec["recipe_root"])
    rows = []
    for task in context["plan"]["runs"]:
        index = task["index"]
        job = scheduler[index]
        base = dict(index=index, method=task["method"], mode=task["mode"], pair=task["pair"], scheduler=job)
        if job["state"] != "COMPLETED" or job["exit_code"] != "0:0":
            rows.append(dict(base, status="failed", evidence=failure_evidence(archive, task),
                             reason="Terminal scheduler task did not complete successfully"))
            continue
        try:
            rows.append(successful_task(archive, context, index, job))
        except FileNotFoundError as error:
            rows.append(dict(base, status="missing_evidence", evidence=failure_evidence(archive, task),
                             error_type=type(error).__name__, reason=str(error)))
        except (ValueError, KeyError, TypeError, OSError, NativeOutputFailure) as error:
            rows.append(dict(base, status="invalid_evidence", evidence=failure_evidence(archive, task),
                             error_type=type(error).__name__, reason=str(error)))
    valid = [row for row in rows if row["status"] == "validated"]
    temporal_issues = []
    for left, right in zip(valid, valid[1:]):
        if left["boot_id"] != right["boot_id"] or left["finished_ns"] >= right["started_ns"]:
            temporal_issues.append(dict(left=left["index"], right=right["index"], reason="changed clock domain or nonsequential spans"))
    arithmetic = summarize(context["plan"], rows)
    for item in [accounting_record, *sources, *recipes, *[item for row in rows for item in row["evidence"]]]:
        check(item)
    return dict(status="native_overhead_panel_audited_not_scientific_admission", runs=rows, paired=arithmetic,
        validated_tasks=len(valid), failed_or_unvalidated_tasks=len(rows)-len(valid), temporal_issues=temporal_issues,
        observed_screens_all_pass=(all(row["whole_command_screen_passed"] and not row["flagged_intervals"] for row in valid)
                                  if len(valid) == 18 else None),
        accounting=accounting_record, frozen_sources=sources, archived_recipe=recipes, source=record(__file__),
        helpers=[record(Path(__file__).with_name(name)) for name in (
            "verify_frontier_overhead_provenance.py", "replay_frontier_overhead_measurement.py",
            "fingerprint_native_overhead_outputs.py", "summarize_frontier_overhead.py", "validate_scaling_outputs.py")],
        scientific_timings_admitted=False, environmental_validity_established=False, publication_ready=False,
        limitations=["Every terminal task is retained; no selective retries, exclusions or overhead subtraction.",
            "Numerical budgets are separate from output equivalence, duration, temporal and screening checks.",
            "Boundary-only tasks lack interval screening even when observed whole-span screens pass.",
            "Runtime identity is supported by retained before/after checks, not a rerun of the remote environment.",
            "Smallest scaling input only; no causal overhead bound or scientific runtime comparison is established."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("archive", "results", "accounting", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--panel", choices=sorted(PANELS), default="frontier_21838")
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.archive.resolve(), args.results.resolve(), args.accounting.resolve(), panel=args.panel)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
