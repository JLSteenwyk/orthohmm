"""Convert scientifically recovered QfO predictions without admitting failed timing."""

import argparse
import json
import os
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_native_factorial_qfo_pairs import (
    ENV_SHA, FIXED_INPUTS, Evidence, check, conversion_kind, convert, filter_pairs,
    load_mapping, normalize_owners, pair_coverage, record, require, universe,
)
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.review_native_factorial_attempt import scheduler_fields
from benchmark_tools.run_native_factorial_cost import (
    ROOT, SCOPE, read, validate_plan, validate_request, verify_terminal,
)
from benchmark_tools.validate_native_factorial_outputs import expected_factors


def admit_recovery(review, request_ref, request, run):
    kind = conversion_kind(run)
    require(review.get("schema") == "native_factorial_measurement_failure_review_v1"
        and review.get("status") == "native_scientific_outputs_recovered_measurement_failure_retained"
        and review.get("request") == request_ref and review.get("plan") == request["plan"]
        and all(review.get(k) == run[k] for k in ("index", "cell", "dataset", "repeat"))
        and review.get("job_id") == request["job_id"]
        and review.get("scheduler_state") == "FAILED" and review.get("scheduler_exit_code") == "1:0"
        and all(review.get(k) is True for k in ("terminal_reviewed", "native_outputs_validated",
            "native_command_success"))
        and all(review.get(k) is False for k in ("primary_resources_replayed",
            "shared_host_resources_reviewed", "scientific_timings_admitted", "eligible_for_timing_comparison",
            "scheduler_success", "original_receipts_rewritten", "inference_reexecuted", "automatic_retry",
            "accuracy_evaluated", "uncontended_timing"))
        and "resources" in review and review["resources"] is None
        and review.get("execution_scope") == SCOPE,
        "Require an explicit bound scientific recovery with failed timing retained")
    return kind


def bind_recovery(request_ref, review_ref):
    request, review = read(request_ref), read(review_ref)
    plan = read(request["plan"])
    validate_request(request, request["plan"], request["job_id"])
    run = validate_plan(plan)[request["index"]]
    kind = admit_recovery(review, request_ref, request, run)
    require(review["source"] == record(ROOT / "benchmark_tools/review_native_factorial_measurement_failure.py"),
        "Scientific recovery source changed")
    terminal = verify_terminal(request["job_id"])
    fields = scheduler_fields(terminal)
    require(fields.get("JobState", fields.get("State")) == "FAILED" and fields.get("ExitCode") == "1:0",
        "Original failed native allocation outcome changed")
    if terminal["source"] == "live_controller":
        require(fields.get("Comment") == request_ref["sha256"], "Original native request comment changed")
    outputs = read(review["outputs"])
    require(outputs.get("schema") == "native_factorial_output_review_v1"
        and outputs.get("status") == "native_outputs_validated" and outputs.get("native_outputs_validated") is True
        and outputs.get("cell") == run["cell"] and outputs.get("input_genes") == run["genes"]
        and outputs.get("factors") == expected_factors(run["cell"])
        and outputs.get("source") == record(ROOT / "benchmark_tools/validate_native_factorial_outputs.py")
        and all(outputs.get(k) is False for k in ("accuracy_evaluated", "resource_measurements_admitted",
            "next_identity_authorized")),
        "Recovered scientific output binding differs")
    records = [request_ref, review_ref, request["plan"], review["source"], review["original_failed_review"],
        review["cadence_diagnosis"], review["environment_report"], review["runtime_readback"], review["outputs"],
        *review["evidence"], outputs["source"], *outputs["checked_files"], plan["baseline"],
        *plan["helper_sources"]]
    # Recheck retained bindings, not the multi-gigabyte raw cadence census or semantic inference.
    unique = {}
    for ref in records:
        require(ref["path"] not in unique or unique[ref["path"]] == ref, "Conflicting recovery evidence pins")
        unique[ref["path"]] = ref
    records = list(unique.values())
    for ref in records:
        check(ref)
    return request, review, plan, run, kind, outputs, terminal, records


def materialize(kind, prediction, destination, inputs, owners, mapping, normalized, expected):
    pairs, filtered = destination / "pairs.partial.tsv", destination / "pairs.qfo.partial.tsv"
    count, command = convert(kind, prediction, pairs, inputs, owners, expected)
    total, retained = filter_pairs(pairs, filtered, mapping)
    require(count == total == retained, "Pair count differs or corrected mapping loses predictions")
    coverage = pair_coverage(filtered, normalized)
    require(coverage["pair_rows"] == count, "Independent converted coverage count differs")
    a, b = record(pairs), record(filtered)
    require(all(a[k] == b[k] for k in ("bytes", "sha256")),
        "Mapping changed corrected predictions")
    return pairs, filtered, count, command, coverage


def prepare(request_ref, review_ref, destination):
    source = record(__file__)
    request, review, plan, run, kind, outputs, terminal, records = bind_recovery(request_ref, review_ref)
    require(os.environ.get("SLURM_CPUS_PER_TASK") == "2"
        and os.environ.get("SLURM_JOB_ID", "").isdigit(), "Require scheduled two-CPU conversion")
    destination = Path(destination)
    require(destination.is_absolute() and destination.resolve() == destination
        and destination.is_relative_to(ROOT) and not destination.is_relative_to(Path(plan["panel_root"])),
        "Require separate direct recovery conversion destination")
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    evidence = Evidence()
    for ref in records:
        require(record(evidence.bind(ref["path"])) == ref, "Recovery evidence changed")
    name, digest = FIXED_INPUTS["qfo_preparation"]
    prepared_ref = record(ROOT / name)
    require(prepared_ref["sha256"] == digest and read(prepared_ref)["input_fastas"] == run["inputs"],
        "Recovered inputs differ from frozen corrected-QfO release")
    evidence.bind(prepared_ref["path"])
    env_ref = record(ROOT / "benchmark_tools/results/qfo_assessment_environment_20260917.json")
    require(env_ref["sha256"] == ENV_SHA, "Frozen QfO reference environment differs")
    evidence.bind(env_ref["path"])
    mappings = [r for r in read(env_ref)["reference_files"] if Path(r["path"]).name == "mapping.json.gz"]
    require(len(mappings) == 1, "Require one frozen QfO accession mapping")
    mapping_ref = mappings[0]
    check(mapping_ref)
    mapping = load_mapping(evidence.bind(mapping_ref["path"]))
    inputs = Path(run["output_root"]) / "input"
    owners, owner_digest, per_species = universe(dict(run, input_directory=str(inputs)), evidence)
    require(owner_digest == outputs["gene_ownership_sha256"] and per_species == outputs["per_species_counts"],
        "Recovered conversion ownership differs")
    normalized = normalize_owners(owners, mapping)
    suffix = ("native/orthohmm_phylogeny/orthohmm_pairwise_orthologs.tsv" if kind == "native"
        else "native/orthohmm_working_res/orthohmm_edges_clustered.txt")
    prediction = Path(run["output_root"]) / suffix
    checked = {r["path"]: r for r in outputs["checked_files"]}
    require(str(prediction) in checked, "Selected prediction was not scientifically recovered")
    prediction_ref = checked[str(prediction)]
    check(prediction_ref)
    require({p.name for p in inputs.glob("*.fasta")} == set(run["native_order"]), "Recovered FASTA inventory differs")
    for path in (Path(__file__), ROOT / "benchmark_tools/prepare_native_factorial_qfo_pairs.py",
        ROOT / "benchmark_tools/prepare_qfo_factorial_pairs.py", ROOT / "benchmark_tools/qfo_filter_pairs.py",
        ROOT / "qfo_benchmark/og_to_pairwise.py"):
        evidence.bind(path)
    destination.mkdir(parents=True, exist_ok=False)
    report = dict(schema="measurement_failed_native_qfo_conversion_v1", status="preparing_unscored",
        source=source, job_id=os.environ["SLURM_JOB_ID"], native_job_id=request["job_id"],
        native_index=run["index"], cell=run["cell"], participant="ohmm_qfo_recovered_native_" + run["cell"],
        request=request_ref, plan=request["plan"], scientific_recovery=review_ref, native_scheduler=terminal,
        prepared=prepared_ref, environment_manifest=env_ref, input_fastas=run["inputs"],
        gene_ownership_sha256=owner_digest, mapping=mapping_ref, native_input=prediction_ref,
        semantics="native phylogenetically inferred pairs" if kind == "native" else "cross-species group-derived clique pairs",
        conversion_kind=kind, accuracy_evaluated=False, native_inference_reexecuted=False,
        automatic_retry=False, next_identity_authorized=False, publication_ready=False,
        original_native_scheduler_success=False, scientific_timings_admitted=False,
        eligible_for_timing_comparison=False, resources=None, conversion_started_monotonic_ns=time.monotonic_ns(),
        conversion_interval_scope="pair materialization/filter/coverage/postflight, outside failed inference measurement")
    save(destination / "preflight.json", report)
    try:
        expected = outputs["phylogeny"]["native_pair_rows"] if kind == "native" else None
        pairs, filtered, count, command, coverage = materialize(kind, prediction, destination, inputs,
            owners, mapping, normalized, expected)
        checked_files = evidence.finish()
        pairs.rename(destination / "pairs.tsv")
        filtered.rename(destination / "pairs.qfo.tsv")
        report.update(status="measurement_failed_native_qfo_pairs_prepared_unscored",
            pairs=record(destination / "pairs.tsv"), filtered_pairs=record(destination / "pairs.qfo.tsv"),
            total_pairs=count, retained_pairs=count, expected_pairs=count, removed_mapping_pairs=0,
            empty_predictions=count == 0, pair_coverage=coverage, command=command, checked_records=checked_files,
            limitations=["Distinct recovered-science conversion; the original allocation and timing remain failed/ineligible.",
                "Only materialization, mapping and coverage; separate native assessment/admission remain required.",
                "Relation coverage is not benchmark accuracy and includes all input accessions in its denominator.",
                "R-on resolved pairs and R-off group cliques have different prediction semantics.",
                "No native inference, semantic review, raw cadence census or cached score is rerun or substituted.",
                "Shared-host effects are unknown and potentially tool-dependent; no isolated ranking or corrected time."])
        if (destination / "conversion.log").exists():
            report["conversion_log"] = record(destination / "conversion.log")
    except BaseException as error:
        report.update(status="measurement_failed_native_qfo_conversion_failed_retained",
            error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report["conversion_finished_monotonic_ns"] = time.monotonic_ns()
        save(destination / "results.json", report)
    return record(destination / "results.json")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("request", "scientific-recovery", "output-directory"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("request-sha256", "scientific-recovery-sha256"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    refs = [record(p) for p in (args.request, args.scientific_recovery)]
    require([r["sha256"] for r in refs] == [args.request_sha256, args.scientific_recovery_sha256],
        "Explicit recovered conversion request/review checksum differs")
    print(json.dumps(prepare(*refs, args.output_directory.absolute()), sort_keys=True))
