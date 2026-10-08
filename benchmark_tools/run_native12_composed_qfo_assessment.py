"""Run unchanged six QfO endpoints with explicit composed native12 lineage."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import time

from benchmark_tools import prepare_native12_composed_qfo_pairs as converter
from benchmark_tools.admit_qfo_corrected_comparator_assessment import accounting
from benchmark_tools.native12_composed_review_binding import native_binding
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.run_native_factorial_cost import ROOT, read
from benchmark_tools.run_native_factorial_qfo_assessment import fas_protocol, helper_records as base_helpers
from benchmark_tools.run_qfo_recovered_assessment import command_for, environment_records
from benchmark_tools.validate_native_factorial_outputs import require


SCHEMA = "native12_composed_qfo_execution_v1"


def validate_stage(stage, producer):
    require(stage.get("schema") == converter.SCHEMA and stage.get("status") == converter.STATUS
        and type(stage.get("native_index")) is int and stage["native_index"] == 12
        and stage.get("native_job_id") == 24036 and stage.get("cell") == "p1_c1_r1"
        and stage.get("participant") == converter.PARTICIPANT and stage.get("conversion_kind") == "native"
        and stage.get("semantics") == converter.SEMANTICS
        and all(stage.get(k) is False for k in ("accuracy_evaluated", "native_inference_reexecuted",
            "automatic_retry", "next_identity_authorized", "publication_ready", "original_review_translated")),
        "Require explicit unscored composed native12 conversion")
    require((producer.get("State"), producer.get("ExitCode"), producer.get("NodeList"), producer.get("AllocCPUS"))
        == ("COMPLETED", "0:0", "bizon", "2") and stage.get("job_id") == producer.get("JobIDRaw")
        and producer.get("ReqMem") in {"32G", "32Gn", "32768M", "32768Mn"},
        "Require successful two-CPU/32GiB conversion producer")
    counts = [stage.get(k) for k in ("total_pairs", "retained_pairs", "expected_pairs", "removed_mapping_pairs")]
    require(all(type(v) is int for v in counts) and 0 <= counts[0] == counts[1] == counts[2]
        and counts[3] == 0 and stage.get("empty_predictions") is (counts[0] == 0), "Invalid conversion pair counts")
    require(isinstance(stage.get("pairs"), dict) and isinstance(stage.get("filtered_pairs"), dict)
        and all(stage["pairs"][k] == stage["filtered_pairs"][k] for k in ("bytes", "sha256")),
        "Reference mapping altered predictions")
    coverage = stage["pair_coverage"]
    require(isinstance(coverage, dict) and type(coverage.get("pair_rows")) is int and coverage["pair_rows"] == counts[0]
        and type(coverage.get("input_accessions")) is int and coverage["input_accessions"] > 0
        and type(coverage.get("accessions_in_any_pair")) is int
        and 0 <= coverage["accessions_in_any_pair"] <= coverage["input_accessions"]
        and type(coverage.get("fraction_inputs_in_any_pair")) in (int, float)
        and coverage["fraction_inputs_in_any_pair"] == coverage["accessions_in_any_pair"] / coverage["input_accessions"],
        "Invalid full-input relation coverage")
    require(type(stage.get("conversion_started_monotonic_ns")) is int
        and type(stage.get("conversion_finished_monotonic_ns")) is int
        and 0 <= stage["conversion_started_monotonic_ns"] <= stage["conversion_finished_monotonic_ns"],
        "Invalid separate conversion interval")
    provenance = stage["composed_binding"]
    require(isinstance(provenance, dict) and provenance.get("composed_schema_preserved") is True
        and provenance.get("original_review_translated") is False
        and provenance.get("next_identity_authorized") is False
        and type(provenance.get("review_producer_job_id")) is int
        and provenance["review_producer_job_id"] > 24036
        and str(producer.get("JobIDRaw", "")).isdigit()
        and int(producer["JobIDRaw"]) > provenance["review_producer_job_id"]
        and isinstance(provenance.get("review_held"), dict)
        and isinstance(provenance.get("review_release"), dict),
        "Final binding was translated, broadened or lacks distinct review/conversion producers")


def helper_records():
    return base_helpers() + [record(ROOT / "benchmark_tools" / name) for name in (
        "native12_composed_review_binding.py", "prepare_native12_composed_qfo_pairs.py",
        "run_native12_composed_qfo_assessment.py", "admit_native12_composed_qfo_assessment.py")]


def conversion_binding(root, pairs_ref, conversion_job):
    require(root == ROOT and pairs_ref["path"] == str(converter.DESTINATION / "results.json"),
        "Require explicit composed conversion result namespace")
    raw, producer = accounting(conversion_job, include_memory=True)
    stage = read(pairs_ref)
    validate_stage(stage, producer)
    require(stage.get("source") == record(converter.__file__)
        and stage.get("binding_source") == record(ROOT / "benchmark_tools/native12_composed_review_binding.py")
        and stage.get("conversion_kernel_source") == record(ROOT / "benchmark_tools/prepare_native_factorial_qfo_pairs.py"),
        "Composed converter/binding or original kernel source differs")
    provenance = stage["composed_binding"]
    request, execution, plan, run, review, outputs, kind, terminal, native_records, binding = \
        native_binding(stage["request"], stage["terminal_review"], provenance["review_producer_job_id"],
            provenance["review_held"], provenance["review_release"])
    require(stage["plan"] == request["plan"] and stage["amendment"] == request["amendment"]
        and stage["native_job_id"] == request["job_id"] and stage["native_index"] == run["index"]
        and stage["cell"] == run["cell"] and stage["input_fastas"] == run["inputs"]
        and stage["conversion_kind"] == kind and kind == "native"
        and stage["native_cpu_ids"] == outputs["native_cpu_ids"]
        and stage["allocated_ready"] == outputs["allocated_ready"] and stage["composed_binding"] == binding
        and stage["expected_pairs"] == outputs["phylogeny"]["native_pair_rows"],
        "Conversion differs from admitted native identity")
    env_ref = record(root / "benchmark_tools/results/qfo_assessment_environment_20260917.json")
    require(env_ref["sha256"] == converter.ENV_SHA and stage["environment_manifest"] == env_ref,
        "Frozen QfO scoring environment differs")
    manifest = read(env_ref)
    require([r for r in manifest["reference_files"] if Path(r["path"]).name == "mapping.json.gz"] == [stage["mapping"]],
        "Conversion/scoring mapping differs")
    require(outputs["gene_ownership_sha256"] == stage["gene_ownership_sha256"]
        and stage["pair_coverage"]["input_accessions"] == run["genes"]
        and stage["native_input"] in outputs["checked_files"]
        and stage["native_input"]["path"] == str(Path(run["output_root"]) /
            "native/orthohmm_phylogeny/orthohmm_pairwise_orthologs.tsv"),
        "Converted prediction or all-input ownership differs")
    records = [pairs_ref, env_ref, stage["source"], stage["binding_source"], stage["conversion_kernel_source"],
        *native_records, *stage["checked_records"], *environment_records(manifest),
        stage["pairs"], stage["filtered_pairs"], *helper_records()]
    for ref in records:
        check(ref)
    return stage, manifest, dict(conversion_accounting=raw, conversion_scheduler=producer,
        native_scheduler=terminal, environment_manifest=env_ref, composed_binding=binding,
        verified_records=records, fas_protocol=fas_protocol(manifest))


def execution_spec(root, pairs_ref, stage, manifest, binding):
    output = root / "benchmarks/results/native12_composed_qfo_assessment_v1"
    work = root / "qfo_benchmark/w/cn12"
    results = root / "qfo_benchmark/scoring/composed_native12_v1"
    return dict(schema=SCHEMA, status="prepared_unrun", native_index=12, native_job_id=24036,
        cell=stage["cell"], stage=stage, amendment=stage["amendment"], source=record(__file__),
        pairs_manifest=pairs_ref, command=command_for(root, stage, manifest, work, results),
        cwd=str(output), work=str(work), results=str(results), environment_overrides=manifest["environment_overrides"],
        assessment_resource_limits=dict(cpu_slots=8, memory_bytes=128 * 1024 ** 3, scheduler_limit_s=93600),
        resource_scope="QfO assessment only; separate from inference, conversion and admission",
        accuracy_admitted=False, publication_ready=False, automatic_retry=False,
        native_inference_reexecuted=False, next_identity_authorized=False, original_review_translated=False, **binding)


def prepare(root, pairs_ref, conversion_job):
    stage, manifest, binding = conversion_binding(root, pairs_ref, conversion_job)
    report = execution_spec(root, pairs_ref, stage, manifest, binding)
    for key in ("cwd", "work", "results"):
        path = Path(report[key])
        require(path.is_absolute() and path.resolve() == path and path.is_relative_to(root), "Require direct fresh scoring paths")
        if path.exists() or path.is_symlink():
            raise FileExistsError(path)
    return report


def run(root, pairs_ref, conversion_job, expected_source_sha256, check_only=False):
    require(record(__file__)["sha256"] == expected_source_sha256, "Prospective scoring worker changed")
    require(check_only or (os.environ.get("SLURM_CPUS_PER_TASK") == "8"
        and os.environ.get("SLURM_JOB_ID", "").isdigit()
        and os.environ.get("SLURM_MEM_PER_NODE") == "131072"), "Require scheduled eight-CPU/128GiB assessment")
    report = prepare(root, pairs_ref, conversion_job)
    if check_only:
        return report
    output = Path(report["cwd"])
    output.mkdir(parents=True, exist_ok=False)
    report.update(status="running", job_id=os.environ["SLURM_JOB_ID"], started_monotonic_ns=time.monotonic_ns(),
        interval_scope="Nextflow execution/postflight; excludes inference/conversion/preflight")
    save(output / "preflight.json", report)
    env = {**os.environ, **report["environment_overrides"],
        "NXF_SINGULARITY_CACHEDIR": str(root / "qfo_benchmark/scoring/container_cache")}
    try:
        with (output / "scoring.log").open("x") as log:
            done = subprocess.run(report["command"], cwd=output, env=env, stdout=log, stderr=subprocess.STDOUT)
        report["exit_code"] = done.returncode
        for ref in [report["source"], *report["verified_records"]]:
            check(ref)
        require(done.returncode == 0, f"Native12 composed QfO assessment failed: {done.returncode}")
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
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--pairs-sha256", required=True)
    parser.add_argument("--conversion-job", required=True)
    parser.add_argument("--source-sha256", required=True)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    ref = record(args.pairs)
    require(ref["sha256"] == args.pairs_sha256, "Conversion manifest checksum differs")
    report = run(ROOT, ref, args.conversion_job, args.source_sha256, args.check_only)
    print(json.dumps(dict(status=report["status"], command=report["command"]), sort_keys=True))
