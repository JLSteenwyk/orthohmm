"""Assess terminal-reviewed full-native QfO pairs in a separate namespace."""

import argparse
import ast
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_qfo_corrected_comparator_assessment import accounting
from benchmark_tools.prepare_native_factorial_qfo_pairs import admit_conversion, conversion_kind, ENV_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.run_native_factorial_cost import read, validate_plan, validate_request, verify_terminal
from benchmark_tools.run_qfo_recovered_assessment import command_for, environment_records
from benchmark_tools.validate_native_factorial_outputs import require


def validate_stage(stage, scheduler):
    run = dict(index=stage.get("native_index"), dataset="qfo_corrected", cell=stage.get("cell"))
    kind = conversion_kind(run)
    semantics = "native phylogenetically inferred pairs" if kind == "native" else "cross-species group-derived clique pairs"
    require(stage.get("schema") == "full_native_factorial_qfo_conversion_v1"
        and stage.get("status") == "full_native_factorial_qfo_pairs_prepared_unscored"
        and stage.get("participant") == "ohmm_qfo_full_native_" + run["cell"]
        and stage.get("conversion_kind") == kind and stage.get("semantics") == semantics
        and all(stage.get(k) is False for k in ("accuracy_evaluated", "native_inference_reexecuted",
            "automatic_retry", "next_identity_authorized", "publication_ready")), "Wrong full-native conversion identity/semantics")
    require((scheduler.get("State"), scheduler.get("ExitCode"), scheduler.get("NodeList"), scheduler.get("AllocCPUS"))
        == ("COMPLETED", "0:0", "bizon", "2") and stage.get("job_id") == scheduler.get("JobIDRaw"),
        "Require successful two-CPU conversion completion")
    counts = [stage.get(k) for k in ("total_pairs", "retained_pairs", "expected_pairs", "removed_mapping_pairs")]
    require(all(type(v) is int for v in counts) and 0 <= counts[0] == counts[1] == counts[2]
        and counts[3] == 0 and stage.get("empty_predictions") is (counts[0] == 0), "Invalid full-native pair counts")
    require(all(stage["pairs"][k] == stage["filtered_pairs"][k] for k in ("bytes", "sha256")),
        "Reference filtering changed corrected predictions")
    coverage = stage["pair_coverage"]
    require(coverage["pair_rows"] == counts[0] and type(coverage["input_accessions"]) is int
        and type(coverage["accessions_in_any_pair"]) is int
        and 0 <= coverage["accessions_in_any_pair"] <= coverage["input_accessions"]
        and coverage["input_accessions"] > 0
        and coverage["fraction_inputs_in_any_pair"] == coverage["accessions_in_any_pair"] / coverage["input_accessions"],
        "Invalid protein relation coverage")


def helper_records():
    names = ("run_native_factorial_qfo_assessment.py", "prepare_native_factorial_qfo_pairs.py",
        "run_qfo_recovered_assessment.py", "admit_qfo_corrected_comparator_assessment.py",
        "verify_ygob_validation.py", "prepare_ob_candidate_neighborhood.py", "probe_dgx_step_separation.py",
        "admit_native_factorial_qfo_assessment.py", "admit_qfo_recovered_assessment.py",
        "audit_qfo_fas_samples.py", "validate_qfo_native_assessment.py", "qfo_summarize_scores.py")
    return [record(Path(__file__).with_name(name)) for name in names]


def fas_protocol(manifest):
    records = {Path(r["path"]).name: r for r in manifest["pipeline_files"]}
    source, pipeline = records["fas_benchmark.py"], records["main.nf"]
    check(source)
    check(pipeline)
    tree = ast.parse(Path(source["path"]).read_text())
    constants = [n.value.value for n in tree.body if isinstance(n, ast.Assign)
        and any(isinstance(t, ast.Name) and t.id == "MAX_PAIRS_COMPUTE" for t in n.targets)
        and isinstance(n.value, ast.Constant)]
    require(constants == [9000], "FAS sampling cap differs from frozen protocol")
    return dict(source=source, pipeline=pipeline, newly_computed_pair_cap=9000,
        population="precomputed predictions plus missing predictions with both feature annotations; not all submitted pairs",
        sampling="unseeded native shuffle; precomputed sample preserves eligible precomputed/missing ratio when missing pairs exist",
        missing_scores="native scorer omits uncomputed/failed missing-pair scores; no zero imputation",
        species_scope="all supplied species; frozen pipeline does not pass --limited-species",
        uncertainty="native pair-IID SEM is not an independent-family or paired-method confidence interval")


def conversion_binding(root, pairs_ref, conversion_job):
    conversion_text, scheduler = accounting(conversion_job)
    stage = read(pairs_ref)
    validate_stage(stage, scheduler)
    require(stage["source"] == record(Path(__file__).with_name("prepare_native_factorial_qfo_pairs.py")),
        "Conversion source differs")
    request, review = read(stage["request"]), read(stage["terminal_review"])
    validate_request(request, request["plan"], request["job_id"])
    plan = read(request["plan"])
    run = validate_plan(plan)[request["index"]]
    kind = admit_conversion(review, stage["request"], request, run)
    require(stage["plan"] == request["plan"] and stage["native_job_id"] == request["job_id"]
        and stage["native_index"] == run["index"] and stage["cell"] == run["cell"]
        and stage["input_fastas"] == run["inputs"] and stage["conversion_kind"] == kind,
        "Conversion differs from native request/terminal review")
    native_scheduler = verify_terminal(request["job_id"])
    fields = native_scheduler["verified"].get("fields", native_scheduler["verified"])
    require(fields.get("JobState", fields.get("State")) == "COMPLETED" and fields.get("ExitCode") == "0:0",
        "Native scheduler no longer reports successful completion")
    if native_scheduler["source"] == "live_controller":
        require(fields.get("Comment") == stage["request"]["sha256"], "Native scheduler request comment differs")
    env_ref = record(root / "benchmark_tools/results/qfo_assessment_environment_20260917.json")
    require(env_ref["sha256"] == ENV_SHA and stage["environment_manifest"] == env_ref, "QfO scorer manifest differs")
    manifest = read(env_ref)
    require([r for r in manifest["reference_files"] if Path(r["path"]).name == "mapping.json.gz"] == [stage["mapping"]],
        "Conversion/scoring reference mapping differs")
    outputs = read(review["reviews"]["outputs_or_failure"])
    require(outputs.get("native_outputs_validated") is True and outputs.get("request") == stage["request"]
        and outputs.get("gene_ownership_sha256") == stage["gene_ownership_sha256"]
        and stage["pair_coverage"]["input_accessions"] == run["genes"]
        and stage["native_input"] in outputs["checked_files"], "Conversion native output/ownership binding differs")
    if kind == "native":
        require(stage["retained_pairs"] == outputs["phylogeny"]["native_pair_rows"],
            "Converted count differs from independently validated native pairs")
    records = [pairs_ref, env_ref, stage["source"], stage["request"], stage["terminal_review"], request["plan"],
        review["source"], review["scheduler"], *review["reviews"].values(), *stage["checked_records"],
        *outputs["checked_files"], *outputs["evidence"], *plan["helper_sources"],
        *environment_records(manifest), stage["pairs"], stage["filtered_pairs"], *helper_records()]
    for item in records:
        check(item)
    return stage, manifest, dict(conversion_accounting=conversion_text, conversion_scheduler=scheduler,
        native_scheduler=native_scheduler, environment_manifest=env_ref, verified_records=records,
        fas_protocol=fas_protocol(manifest))


def execution_spec(root, pairs_ref, stage, manifest, binding):
    index = stage["native_index"]
    output = root / "benchmarks/results/full_native_qfo_assessment_v1" / stage["cell"]
    work = root / "qfo_benchmark/w" / f"nq{index:02d}"
    results = root / "qfo_benchmark/scoring" / f"full_native_{index:02d}"
    return dict(schema="full_native_factorial_qfo_execution_v1", status="prepared_unrun",
        native_index=index, native_job_id=stage["native_job_id"], cell=stage["cell"], stage=stage,
        source=record(__file__), pairs_manifest=pairs_ref, command=command_for(root, stage, manifest, work, results),
        cwd=str(output), work=str(work), results=str(results), environment_overrides=manifest["environment_overrides"],
        accuracy_admitted=False, publication_ready=False, automatic_retry=False, native_inference_reexecuted=False,
        next_identity_authorized=False, **binding)


def prepare(root, pairs_ref, conversion_job):
    stage, manifest, binding = conversion_binding(root, pairs_ref, conversion_job)
    report = execution_spec(root, pairs_ref, stage, manifest, binding)
    for key in ("cwd", "work", "results"):
        path = Path(report[key])
        require(path.is_absolute() and path.resolve() == path and path.is_relative_to(root), "Require direct assessment paths")
        if path.exists() or path.is_symlink():
            raise FileExistsError(path)
    return report


def run(root, pairs_ref, conversion_job, check_only=False):
    require(check_only or (os.environ.get("SLURM_CPUS_PER_TASK") == "8"
        and os.environ.get("SLURM_JOB_ID", "").isdigit()), "Require scheduled eight-CPU assessment")
    report = prepare(root, pairs_ref, conversion_job)
    if check_only:
        return report
    output = Path(report["cwd"])
    output.mkdir(parents=True, exist_ok=False)
    report.update(status="running", job_id=os.environ["SLURM_JOB_ID"], started_monotonic_ns=time.monotonic_ns(),
        interval_scope="Nextflow endpoint execution and postflight; excludes conversion and preflight")
    save(output / "preflight.json", report)
    env = {**os.environ, **report["environment_overrides"],
        "NXF_SINGULARITY_CACHEDIR": str(root / "qfo_benchmark/scoring/container_cache")}
    try:
        with (output / "scoring.log").open("x") as log:
            done = subprocess.run(report["command"], cwd=output, env=env, stdout=log, stderr=subprocess.STDOUT)
        report["exit_code"] = done.returncode
        for item in [report["source"], *report["verified_records"]]:
            check(item)
        require(done.returncode == 0, f"Native QfO assessment failed: {done.returncode}")
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
    require(pairs_ref["sha256"] == args.pairs_sha256, "Pair preparation checksum differs")
    result = run(args.root.resolve(), pairs_ref, args.conversion_job, args.check_only)
    print(json.dumps(dict(status=result["status"], command=result["command"]), sort_keys=True))
