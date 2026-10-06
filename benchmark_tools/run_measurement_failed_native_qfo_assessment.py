"""Score recovered native QfO predictions while preserving failed timing provenance."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_qfo_corrected_comparator_assessment import accounting
from benchmark_tools.prepare_measurement_failed_native_qfo_pairs import bind_recovery
from benchmark_tools.prepare_native_factorial_qfo_pairs import ENV_SHA, conversion_kind
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.run_native_factorial_cost import read
from benchmark_tools.run_native_factorial_qfo_assessment import fas_protocol, helper_records
from benchmark_tools.run_qfo_recovered_assessment import command_for, environment_records
from benchmark_tools.validate_native_factorial_outputs import require


def validate_stage(stage, scheduler):
    identity = dict(index=stage.get("native_index"), dataset="qfo_corrected", cell=stage.get("cell"))
    kind = conversion_kind(identity)
    semantics = "native phylogenetically inferred pairs" if kind == "native" else "cross-species group-derived clique pairs"
    require(stage.get("schema") == "measurement_failed_native_qfo_conversion_v1"
        and stage.get("status") == "measurement_failed_native_qfo_pairs_prepared_unscored"
        and stage.get("participant") == "ohmm_qfo_recovered_native_" + identity["cell"]
        and stage.get("conversion_kind") == kind and stage.get("semantics") == semantics
        and all(stage.get(k) is False for k in ("accuracy_evaluated", "native_inference_reexecuted",
            "automatic_retry", "next_identity_authorized", "publication_ready", "original_native_scheduler_success",
            "scientific_timings_admitted", "eligible_for_timing_comparison"))
        and "resources" in stage and stage["resources"] is None,
        "Require a distinct unscored recovered conversion with failed timing retained")
    require((scheduler.get("State"), scheduler.get("ExitCode"), scheduler.get("NodeList"), scheduler.get("AllocCPUS"))
        == ("COMPLETED", "0:0", "bizon", "2") and stage.get("job_id") == scheduler.get("JobIDRaw"),
        "Require successful two-CPU recovered conversion")
    counts = [stage.get(k) for k in ("total_pairs", "retained_pairs", "expected_pairs", "removed_mapping_pairs")]
    require(all(type(v) is int for v in counts) and 0 <= counts[0] == counts[1] == counts[2]
        and counts[3] == 0 and stage.get("empty_predictions") is (counts[0] == 0), "Invalid recovered pair counts")
    require(all(stage["pairs"][k] == stage["filtered_pairs"][k] for k in ("bytes", "sha256")),
        "Recovered reference filtering changed predictions")
    coverage = stage["pair_coverage"]
    require(type(coverage.get("pair_rows")) is int and coverage["pair_rows"] == counts[0]
        and type(coverage.get("input_accessions")) is int and coverage["input_accessions"] > 0
        and type(coverage.get("accessions_in_any_pair")) is int
        and 0 <= coverage["accessions_in_any_pair"] <= coverage["input_accessions"]
        and coverage["fraction_inputs_in_any_pair"] == coverage["accessions_in_any_pair"] / coverage["input_accessions"],
        "Invalid recovered protein relation coverage")


def conversion_binding(root, pairs_ref, conversion_job):
    conversion_text, scheduler = accounting(conversion_job)
    stage = read(pairs_ref)
    validate_stage(stage, scheduler)
    require(stage["source"] == record(Path(__file__).with_name("prepare_measurement_failed_native_qfo_pairs.py")),
        "Recovered conversion source changed")
    request, review, plan, run, kind, outputs, terminal, recovery_records = bind_recovery(
        stage["request"], stage["scientific_recovery"])
    require(stage["plan"] == request["plan"] and stage["native_job_id"] == request["job_id"]
        and stage["native_index"] == run["index"] and stage["cell"] == run["cell"]
        and stage["input_fastas"] == run["inputs"] and stage["conversion_kind"] == kind
        and stage["gene_ownership_sha256"] == outputs["gene_ownership_sha256"]
        and stage["pair_coverage"]["input_accessions"] == run["genes"]
        and stage["native_input"] in outputs["checked_files"], "Recovered conversion/native output binding differs")
    if kind == "native":
        require(stage["retained_pairs"] == outputs["phylogeny"]["native_pair_rows"],
            "Recovered converted count differs from native resolved pairs")
    env_ref = record(root / "benchmark_tools/results/qfo_assessment_environment_20260917.json")
    require(env_ref["sha256"] == ENV_SHA and stage["environment_manifest"] == env_ref,
        "Recovered QfO scoring environment differs")
    manifest = read(env_ref)
    require([r for r in manifest["reference_files"] if Path(r["path"]).name == "mapping.json.gz"] == [stage["mapping"]],
        "Recovered conversion/scoring mapping differs")
    helpers = [record(Path(__file__).with_name(name)) for name in (
        "run_measurement_failed_native_qfo_assessment.py", "admit_measurement_failed_native_qfo_assessment.py",
        "prepare_measurement_failed_native_qfo_pairs.py", "review_native_factorial_measurement_failure.py")]
    records = [pairs_ref, env_ref, *recovery_records, stage["source"], *stage["checked_records"],
        stage["pairs"], stage["filtered_pairs"], *environment_records(manifest), *helper_records(), *helpers]
    unique = {}
    for ref in records:
        require(ref["path"] not in unique or unique[ref["path"]] == ref, "Conflicting recovered scoring evidence pins")
        unique[ref["path"]] = ref
    records = list(unique.values())
    for ref in records:
        check(ref)
    return stage, manifest, dict(conversion_accounting=conversion_text, conversion_scheduler=scheduler,
        native_scheduler=terminal, environment_manifest=env_ref, verified_records=records,
        fas_protocol=fas_protocol(manifest))


def execution_spec(root, pairs_ref, stage, manifest, binding):
    index = stage["native_index"]
    output = root / "benchmarks/results/measurement_failed_native_qfo_assessment_v1" / stage["cell"]
    work = root / "qfo_benchmark/w" / f"mr{index:02d}"
    results = root / "qfo_benchmark/scoring" / f"recovered_native_{index:02d}"
    return dict(schema="measurement_failed_native_qfo_execution_v1", status="prepared_unrun",
        native_index=index, native_job_id=stage["native_job_id"], cell=stage["cell"], stage=stage,
        source=record(__file__), pairs_manifest=pairs_ref, command=command_for(root, stage, manifest, work, results),
        cwd=str(output), work=str(work), results=str(results), environment_overrides=manifest["environment_overrides"],
        accuracy_admitted=False, publication_ready=False, automatic_retry=False, native_inference_reexecuted=False,
        next_identity_authorized=False, original_native_scheduler_success=False, scientific_timings_admitted=False,
        eligible_for_timing_comparison=False, resources=None, **binding)


def prepare(root, pairs_ref, conversion_job):
    stage, manifest, binding = conversion_binding(root, pairs_ref, conversion_job)
    report = execution_spec(root, pairs_ref, stage, manifest, binding)
    for key in ("cwd", "work", "results"):
        path = Path(report[key])
        require(path.is_absolute() and path.resolve() == path and path.is_relative_to(root),
            "Require direct recovered assessment paths")
        if path.exists() or path.is_symlink():
            raise FileExistsError(path)
    return report


def run(root, pairs_ref, conversion_job, check_only=False):
    require(check_only or (os.environ.get("SLURM_CPUS_PER_TASK") == "8"
        and os.environ.get("SLURM_JOB_ID", "").isdigit()), "Require scheduled eight-CPU recovered assessment")
    report = prepare(root, pairs_ref, conversion_job)
    if check_only:
        return report
    output = Path(report["cwd"])
    output.mkdir(parents=True, exist_ok=False)
    report.update(status="running", job_id=os.environ["SLURM_JOB_ID"], started_monotonic_ns=time.monotonic_ns(),
        interval_scope="Nextflow endpoint execution and postflight; separate from conversion and failed inference timing")
    save(output / "preflight.json", report)
    env = {**os.environ, **report["environment_overrides"],
        "NXF_SINGULARITY_CACHEDIR": str(root / "qfo_benchmark/scoring/container_cache")}
    try:
        with (output / "scoring.log").open("x") as log:
            done = subprocess.run(report["command"], cwd=output, env=env, stdout=log, stderr=subprocess.STDOUT)
        report["exit_code"] = done.returncode
        for ref in [report["source"], *report["verified_records"]]:
            check(ref)
        require(done.returncode == 0, f"Recovered native QfO assessment failed: {done.returncode}")
        report["status"] = "process_succeeded_pending_independent_admission"
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report["finished_monotonic_ns"] = time.monotonic_ns()
        report["outputs"] = [record(p) for p in sorted(Path(report["results"]).rglob("*")) if p.is_file()]
        if (output / "scoring.log").exists():
            report["log"] = record(output / "scoring.log")
        save(output / "results.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--pairs-sha256", required=True)
    parser.add_argument("--conversion-job", required=True)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    pairs_ref = record(args.pairs)
    require(pairs_ref["sha256"] == args.pairs_sha256, "Recovered pair-manifest checksum differs")
    result = run(args.root.resolve(), pairs_ref, args.conversion_job, args.check_only)
    print(json.dumps(dict(status=result["status"], command=result["command"]), sort_keys=True))
