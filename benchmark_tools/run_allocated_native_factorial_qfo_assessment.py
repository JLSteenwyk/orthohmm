"""Run frozen QfO endpoints for the prospective allocation-aware native route."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_qfo_corrected_comparator_assessment import accounting
from benchmark_tools.native_factorial_allocated_execution import ROOT
from benchmark_tools.prepare_allocated_native_factorial_qfo_pairs import native_binding, ENV_SHA
from benchmark_tools.prepare_native_factorial_qfo_pairs import conversion_kind
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.run_native_factorial_cost import read
from benchmark_tools.run_native_factorial_qfo_assessment import fas_protocol, helper_records as historical_helpers
from benchmark_tools.run_qfo_recovered_assessment import command_for, environment_records
from benchmark_tools.validate_native_factorial_outputs import require


def validate_stage(stage, scheduler):
    run = dict(index=stage.get("native_index"), dataset="qfo_corrected", cell=stage.get("cell"))
    kind = conversion_kind(run)
    semantics = "native phylogenetically inferred pairs" if kind == "native" else "cross-species group-derived clique pairs"
    require(run["index"] in (10, 11, 12)
        and stage.get("schema") == "allocated_native_factorial_qfo_conversion_v1"
        and stage.get("status") == "allocated_native_factorial_qfo_pairs_prepared_unscored"
        and stage.get("participant") == "ohmm_qfo_full_native_" + run["cell"]
        and stage.get("conversion_kind") == kind and stage.get("semantics") == semantics
        and all(stage.get(k) is False for k in ("accuracy_evaluated", "native_inference_reexecuted",
            "automatic_retry", "next_identity_authorized", "publication_ready")),
        "Wrong allocated native conversion identity/semantics")
    require((scheduler.get("State"), scheduler.get("ExitCode"), scheduler.get("NodeList"), scheduler.get("AllocCPUS"))
        == ("COMPLETED", "0:0", "bizon", "2") and stage.get("job_id") == scheduler.get("JobIDRaw"),
        "Require successful two-CPU conversion completion")
    counts = [stage.get(k) for k in ("total_pairs", "retained_pairs", "expected_pairs", "removed_mapping_pairs")]
    require(all(type(v) is int for v in counts) and 0 <= counts[0] == counts[1] == counts[2]
        and counts[3] == 0 and stage.get("empty_predictions") is (counts[0] == 0), "Invalid native pair counts")
    require(all(stage["pairs"][k] == stage["filtered_pairs"][k] for k in ("bytes", "sha256")),
        "Reference filtering changed corrected predictions")
    coverage = stage["pair_coverage"]
    require(type(coverage["pair_rows"]) is int and coverage["pair_rows"] == counts[0]
        and type(coverage["input_accessions"]) is int and type(coverage["accessions_in_any_pair"]) is int
        and 0 <= coverage["accessions_in_any_pair"] <= coverage["input_accessions"]
        and coverage["input_accessions"] > 0
        and type(coverage["fraction_inputs_in_any_pair"]) in (int, float)
        and coverage["fraction_inputs_in_any_pair"] == coverage["accessions_in_any_pair"] / coverage["input_accessions"],
        "Invalid protein relation coverage")
    require(type(stage.get("conversion_started_monotonic_ns")) is int
        and type(stage.get("conversion_finished_monotonic_ns")) is int
        and 0 <= stage["conversion_started_monotonic_ns"] <= stage["conversion_finished_monotonic_ns"],
        "Invalid separate conversion interval")


def helper_records():
    return historical_helpers() + [record(Path(__file__).with_name(name)) for name in (
        "prepare_allocated_native_factorial_qfo_pairs.py", "run_allocated_native_factorial_qfo_assessment.py",
        "admit_allocated_native_factorial_qfo_assessment.py")]


def conversion_binding(root, pairs_ref, conversion_job):
    require(root == ROOT, "Require frozen allocation-aware repository root")
    conversion_text, scheduler = accounting(conversion_job)
    stage = read(pairs_ref)
    validate_stage(stage, scheduler)
    require(stage["source"] == record(Path(__file__).with_name("prepare_allocated_native_factorial_qfo_pairs.py"))
        and stage.get("conversion_kernel_source") == record(Path(__file__).with_name("prepare_native_factorial_qfo_pairs.py")),
        "Allocated conversion source/kernel differs")
    request, execution, plan, run, review, outputs, kind, terminal, native_records = \
        native_binding(stage["request"], stage["terminal_review"])
    require(stage["plan"] == request["plan"] and stage.get("amendment") == request["amendment"]
        and stage["native_job_id"] == request["job_id"] and stage["native_index"] == run["index"]
        and stage["cell"] == run["cell"] and stage["input_fastas"] == run["inputs"]
        and stage["conversion_kind"] == kind and stage.get("native_cpu_ids") == outputs["native_cpu_ids"]
        and stage.get("allocated_ready") == outputs["allocated_ready"], "Conversion differs from allocated native request")
    env_ref = record(root / "benchmark_tools/results/qfo_assessment_environment_20260917.json")
    require(env_ref["sha256"] == ENV_SHA and stage["environment_manifest"] == env_ref, "QfO scorer manifest differs")
    manifest = read(env_ref)
    require([r for r in manifest["reference_files"] if Path(r["path"]).name == "mapping.json.gz"] == [stage["mapping"]],
        "Conversion/scoring reference mapping differs")
    require(outputs["gene_ownership_sha256"] == stage["gene_ownership_sha256"]
        and stage["pair_coverage"]["input_accessions"] == run["genes"]
        and stage["native_input"] in outputs["checked_files"], "Conversion native output/ownership binding differs")
    suffix = "native/orthohmm_phylogeny/orthohmm_pairwise_orthologs.tsv" if kind == "native" else \
        "native/orthohmm_working_res/orthohmm_edges_clustered.txt"
    require(stage["native_input"]["path"] == str(Path(run["output_root"]) / suffix), "Converted prediction path differs")
    if kind == "native":
        require(stage["retained_pairs"] == outputs["phylogeny"]["native_pair_rows"],
            "Converted count differs from independently validated native pairs")
    records = [pairs_ref, env_ref, stage["source"], stage["conversion_kernel_source"], *native_records,
        *stage["checked_records"], *environment_records(manifest), stage["pairs"], stage["filtered_pairs"], *helper_records()]
    for item in records:
        check(item)
    return stage, manifest, dict(conversion_accounting=conversion_text, conversion_scheduler=scheduler,
        native_scheduler=terminal, environment_manifest=env_ref, verified_records=records, fas_protocol=fas_protocol(manifest))


def execution_spec(root, pairs_ref, stage, manifest, binding):
    index = stage["native_index"]
    output = root / "benchmarks/results/allocated_native_qfo_assessment_v1" / stage["cell"]
    work = root / "qfo_benchmark/w" / f"aq{index:02d}"
    results = root / "qfo_benchmark/scoring" / f"allocated_native_{index:02d}"
    return dict(schema="allocated_native_factorial_qfo_execution_v1", status="prepared_unrun",
        native_index=index, native_job_id=stage["native_job_id"], cell=stage["cell"], stage=stage,
        amendment=stage["amendment"], source=record(__file__), pairs_manifest=pairs_ref,
        command=command_for(root, stage, manifest, work, results), cwd=str(output), work=str(work), results=str(results),
        environment_overrides=manifest["environment_overrides"], accuracy_admitted=False, publication_ready=False,
        automatic_retry=False, native_inference_reexecuted=False, next_identity_authorized=False, **binding)


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
