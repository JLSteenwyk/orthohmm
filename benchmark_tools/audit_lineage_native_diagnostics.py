"""Audit all three terminal lineage-native diagnostics; never admit scientific timing."""

import argparse
import gzip
import json
from pathlib import Path

from benchmark_tools.verify_lineage_native_provenance import load_context, verify, ROOT, JOBS
from benchmark_tools.replay_lineage_native_measurement import replay
from benchmark_tools.audit_frontier_overhead import recipe_evidence
from benchmark_tools.capture_job_scheduler import terminal_record
from benchmark_tools.fingerprint_native_overhead_outputs import fingerprint
from benchmark_tools.validate_scaling_outputs import validate
from benchmark_tools.validate_simulation_outputs import NativeOutputFailure
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check

PRIOR_SHA = "d656ae37dcb64d617a132c82391745e33bfd5f078f6c1df7186156f0022ea006"


def inventory(directory):
    if directory.is_symlink():
        raise ValueError("Symlinked native archive root")
    paths = sorted(directory.rglob("*"))
    if any(path.is_symlink() for path in paths):
        raise ValueError("Symlink in native archive")
    return [record(path) for path in paths if path.is_file()]


def prior_results(path):
    source = record(path)
    if source["sha256"] != PRIOR_SHA:
        raise ValueError("Prior comparison audit differs from pinned panel")
    result = json.loads(gzip.decompress(path.read_bytes()))
    check(source)
    return result, source


def successful_task(archive, context, task, scheduler, prior):
    directory = archive / Path(task["run"]["measurement_directory"]).parent.relative_to(ROOT)
    snapshots = inventory(directory)
    preparation, verification, receipt = [json.loads((directory / name).read_text()) for name in (
        "preparation.json", "verification.json", "lineage_diagnostic_task.json")]
    measured = json.loads((directory / "measurement/lineage_report.json").read_text())
    binding = verify(context, task["index"], preparation, verification, receipt, measured, scheduler)
    replayed = replay(directory / "measurement", binding["job_id"], binding["measured_argv"])
    native = measured["native"]
    adapted = dict(command=binding["measured_argv"], cwd=task["run"]["cwd"],
                   exit_code=native["exit_code"], timed_out=native["timed_out"])
    checked = validate(binding["run"], adapted, {ROOT: archive})
    canonical = fingerprint(binding["run"], {ROOT: archive})
    checked_evidence = [*checked["checked_files"], *canonical["evidence"], checked["gnu_time_companion"]["source"]]
    previous = [r for r in prior["runs"] if r["index"] == task["original_index"]]
    if (len(previous) != 1 or previous[0]["status"] != "validated"
            or previous[0]["method"] != task["method"] or previous[0]["mode"] != "periodic"):
        raise ValueError("Prior same-method periodic comparison is missing or invalid")
    previous = previous[0]
    for item in previous["evidence"]:
        check(item)
    equivalent = canonical["identity"] == previous["work_identity"]
    for item in [*snapshots, *checked_evidence]:
        check(item)
    if snapshots != inventory(directory):
        raise ValueError("Native archive inventory changed during audit")
    original = replayed["screening"]["original_screening"]["original_threshold_screen"]
    return dict(index=task["index"], method=task["method"], job_id=binding["job_id"],
        status="validated" if equivalent else "output_mismatch", output_equivalent=equivalent,
        prior_index=task["original_index"], prior_work_identity=previous["work_identity"],
        work_identity=canonical["identity"], native_wall_s=replayed["native_wall_s"],
        original_whole_screen_passed=original["whole_command_screen"]["screen_passed"],
        original_flagged_intervals=replayed["original_flagged_intervals"],
        narrow_flagged_intervals=replayed["narrow_flagged_intervals"],
        all_narrow_intervals_pass=not replayed["narrow_flagged_intervals"],
        screening=replayed["screening"], memory=replayed["memory"],
        gnu_time=checked["gnu_time_companion"]["accounting"],
        native_counts={k: checked[k] for k in ("input_genes", "orthogroups", "root_hogs", "checkpoint_groups", "native_pair_rows") if k in checked},
        observation_points=len(measured["points"]), started_ns=native["started_ns"], finished_ns=native["finished_ns"],
        boot_id=measured["points"][0]["lineage"]["boot_before"], inventory=snapshots,
        checked_evidence=checked_evidence,
        prior_evidence=previous["evidence"], scientific_timings_admitted=False)


def audit(archive, results, scheduler_directory, prior_path):
    # No native archive reads before every assigned job has a detailed terminal record.
    scheduler_paths = [scheduler_directory / f"scheduler_{job}.txt" for job in JOBS]
    scheduler_sources = [record(path) for path in scheduler_paths]
    schedulers = [path.read_text() for path in scheduler_paths]
    if any(terminal_record(raw, job) is None for raw, job in zip(schedulers, JOBS)):
        raise ValueError("All three jobs must be terminal before native inspection")
    context = load_context(results)
    sources = [record(results / name) for name in (
        "dgx_lineage_native_plan_20260919.json", "dgx_lineage_native_recipe_20260919.json",
        "dgx_pressure_overhead_plan_v2_20260919.json", "LINEAGE_NATIVE_PROTOCOL_20260919.md")]
    recipes = recipe_evidence(archive, context["recipe"], "lineage_native_recipe_v1")
    prior, prior_source = prior_results(prior_path)
    rows = []
    for task, scheduler in zip(context["plan"]["runs"], schedulers):
        try:
            rows.append(successful_task(archive, context, task, scheduler, prior))
        except (OSError, ValueError, KeyError, TypeError, NativeOutputFailure) as error:
            directory = archive / Path(task["run"]["measurement_directory"]).parent.relative_to(ROOT)
            rows.append(dict(index=task["index"], method=task["method"], job_id=JOBS[task["index"]],
                status="failed_or_invalid", error_type=type(error).__name__, reason=str(error),
                retained_files=[str(p) for p in sorted(directory.rglob("*")) if p.is_file()],
                scientific_timings_admitted=False))
    observed = [r for r in rows if r["status"] in {"validated", "output_mismatch"}]
    temporal_issues = [dict(left=a["index"], right=b["index"])
        for a, b in zip(observed, observed[1:]) if a["boot_id"] != b["boot_id"] or a["finished_ns"] >= b["started_ns"]]
    for item in [*scheduler_sources, *sources, *recipes, prior_source]:
        check(item)
    return dict(status="lineage_native_diagnostics_audited_not_scientific_admission", runs=rows,
        validated_tasks=sum(r["status"] == "validated" for r in rows), temporal_issues=temporal_issues,
        all_outputs_equivalent=(all(r["output_equivalent"] for r in observed) if len(observed) == 3 else None),
        all_narrow_intervals_pass=(all(r["all_narrow_intervals_pass"] for r in observed) if len(observed) == 3 else None),
        scheduler=scheduler_sources, frozen_sources=sources, archived_recipe=recipes, prior_audit=prior_source,
        source=record(__file__), helpers=[record(Path(__file__).with_name(name)) for name in (
            "verify_lineage_native_provenance.py", "replay_lineage_native_measurement.py",
            "fingerprint_native_overhead_outputs.py", "validate_scaling_outputs.py",
            "measure_native_lineage_step.py", "measure_native_hierarchy_step.py",
            "probe_cgroup_lineage.py", "probe_native_pressure.py", "audit_frontier_overhead.py", "capture_job_scheduler.py")],
        scientific_timings_admitted=False, environmental_validity_established=False, publication_ready=False,
        limitations=["Every prespecified task retained; no selective reruns, overhead subtraction or threshold changes.",
            "Canonical output equivalence is not biological accuracy or identical internal computational work.",
            "Inventory hashing preserves every regular file; semantic validation is limited to checked native products.",
            "Retained runtime checks cannot exclude temporary changes during execution.",
            "This diagnostic does not establish overhead, repeatability, larger-input behavior or absence of interference."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("archive", "results", "scheduler-directory", "prior-audit", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.archive.resolve(), args.results.resolve(), args.scheduler_directory.resolve(), args.prior_audit.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
