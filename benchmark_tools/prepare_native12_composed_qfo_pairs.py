"""Convert the final ablation's reviewed native ortholog pairs, never cliques."""

import argparse
import json
import os
from pathlib import Path
import time

from benchmark_tools.native12_composed_review_binding import native_binding
from benchmark_tools.prepare_native_factorial_qfo_pairs import (
    normalize_owners, pair_coverage, convert, FIXED_INPUTS, ENV_SHA)
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.qfo_filter_pairs import filter_pairs, load_mapping
from benchmark_tools.run_native_factorial_cost import ROOT, read
from benchmark_tools.validate_native_factorial_outputs import Evidence, require, universe


SCHEMA = "native12_composed_qfo_conversion_v1"
STATUS = "native12_composed_qfo_pairs_prepared_unscored"
DESTINATION = ROOT / "benchmarks/results/native12_composed_qfo_pairs_v1"
PARTICIPANT = "ohmm_qfo_full_native_p1_c1_r1"
SEMANTICS = "native phylogenetically inferred pairs"


def materialize(prediction, inputs, owners, mapping, destination, expected):
    normalized = normalize_owners(owners, mapping)
    pairs, filtered = destination / "pairs.partial.tsv", destination / "pairs.qfo.partial.tsv"
    require(not pairs.exists() and not pairs.is_symlink()
        and not filtered.exists() and not filtered.is_symlink(), "Require fresh native conversion outputs")
    count, command = convert("native", prediction, pairs, inputs, owners, expected)
    total, retained = filter_pairs(pairs, filtered, mapping)
    require(count == total == retained == expected, "Native pair count differs or mapping loses predictions")
    coverage = pair_coverage(filtered, normalized)
    require(coverage["pair_rows"] == count, "Independent native coverage row count differs")
    require(record(pairs)["sha256"] == record(filtered)["sha256"]
        and pairs.stat().st_size == filtered.stat().st_size, "Mapping changed native predictions")
    return dict(total_pairs=total, retained_pairs=retained, expected_pairs=count,
        removed_mapping_pairs=0, empty_predictions=count == 0, pair_coverage=coverage, command=command)


def prepare(request_ref, review_ref, producer_job, held_ref, release_ref, destination, expected_source_sha256):
    source_ref = record(__file__)
    require(source_ref["sha256"] == expected_source_sha256, "Prospective native converter changed")
    require(os.environ.get("SLURM_CPUS_PER_TASK") == "2"
        and os.environ.get("SLURM_JOB_ID", "").isdigit()
        and os.environ.get("SLURM_MEM_PER_NODE") == "32768", "Require scheduled two-CPU/32GiB conversion")
    destination = Path(destination)
    require(destination == DESTINATION and destination.is_absolute()
        and destination.resolve() == destination and not destination.exists() and not destination.is_symlink(),
        "Require separate fresh final conversion destination")
    request, execution, plan, run, review, outputs, kind, terminal, records, binding = native_binding(
        request_ref, review_ref, producer_job, held_ref, release_ref)
    require(kind == "native" and run["index"] == 12 and run["cell"] == "p1_c1_r1",
        "Final native conversion must use R-on native ortholog pairs")
    expected = outputs.get("phylogeny", {}).get("native_pair_rows")
    require(type(expected) is int and expected >= 0, "Require independently reviewed native pair count")
    evidence = Evidence()
    for ref in records:
        evidence.bind(ref["path"])
    name, digest = FIXED_INPUTS["qfo_preparation"]
    prepared_ref = record(ROOT / name)
    require(prepared_ref["sha256"] == digest and read(prepared_ref)["input_fastas"] == run["inputs"],
        "Final native inputs differ from frozen corrected QfO release")
    evidence.bind(prepared_ref["path"])
    env_ref = record(ROOT / "benchmark_tools/results/qfo_assessment_environment_20260917.json")
    require(env_ref["sha256"] == ENV_SHA, "Frozen QfO environment differs")
    environment = read(env_ref)
    evidence.bind(env_ref["path"])
    mappings = [ref for ref in environment["reference_files"] if Path(ref["path"]).name == "mapping.json.gz"]
    require(len(mappings) == 1, "Require one frozen accession mapping")
    mapping_ref = mappings[0]
    check(mapping_ref)
    mapping = load_mapping(evidence.bind(mapping_ref["path"]))
    inputs = Path(run["output_root"]) / "input"
    owners, owner_digest, per_species = universe(dict(run, input_directory=str(inputs)), evidence)
    require(owner_digest == outputs["gene_ownership_sha256"] and per_species == outputs["per_species_counts"],
        "Final conversion gene ownership differs from semantic review")
    prediction = Path(run["output_root"]) / "native/orthohmm_phylogeny/orthohmm_pairwise_orthologs.tsv"
    checked = {ref["path"]: ref for ref in outputs["checked_files"]}
    require(str(prediction) in checked, "Native ortholog prediction was not independently reviewed")
    prediction_ref = checked[str(prediction)]
    check(prediction_ref)
    evidence.bind(prediction_ref["path"])
    require({path.name for path in inputs.glob("*.fasta")} == set(run["native_order"]),
        "Native converter input glob differs from native enumeration")
    for path in [Path(__file__), *[ROOT / "benchmark_tools" / name for name in (
        "native12_composed_review_binding.py", "prepare_native_factorial_qfo_pairs.py",
        "prepare_qfo_factorial_pairs.py", "qfo_filter_pairs.py", "simulation_method_outputs.py",
        "validate_native_factorial_outputs.py")]]:
        evidence.bind(path)
    require(record(__file__) == source_ref, "Native converter source changed")
    destination.mkdir(parents=True, exist_ok=False)
    report = dict(schema=SCHEMA, status="preparing_unscored", source=source_ref,
        conversion_kernel_source=record(ROOT / "benchmark_tools/prepare_native_factorial_qfo_pairs.py"),
        binding_source=record(ROOT / "benchmark_tools/native12_composed_review_binding.py"),
        job_id=os.environ["SLURM_JOB_ID"], native_job_id=request["job_id"], native_index=12, cell=run["cell"],
        participant=PARTICIPANT, request=request_ref, plan=request["plan"], amendment=request["amendment"],
        terminal_review=review_ref, native_scheduler=terminal, native_cpu_ids=outputs["native_cpu_ids"],
        allocated_ready=outputs["allocated_ready"], prepared=prepared_ref, environment_manifest=env_ref,
        input_fastas=run["inputs"], gene_ownership_sha256=owner_digest, mapping=mapping_ref,
        native_input=prediction_ref, semantics=SEMANTICS, conversion_kind="native", composed_binding=binding,
        accuracy_evaluated=False, native_inference_reexecuted=False, automatic_retry=False,
        next_identity_authorized=False, publication_ready=False, original_review_translated=False,
        conversion_started_monotonic_ns=time.monotonic_ns(),
        conversion_interval_scope="pair materialization/filter/coverage/postflight; excludes review and ownership indexing")
    save(destination / "preflight.json", report)
    try:
        converted = materialize(prediction, inputs, owners, mapping, destination, expected)
        checked_files = evidence.finish()
        (destination / "pairs.partial.tsv").rename(destination / "pairs.tsv")
        (destination / "pairs.qfo.partial.tsv").rename(destination / "pairs.qfo.tsv")
        report.update(converted, status=STATUS, pairs=record(destination / "pairs.tsv"),
            filtered_pairs=record(destination / "pairs.qfo.tsv"), checked_records=checked_files,
            limitations=[
                "Explicit final-identity lineage, not ordinary-review or replay schema substitution.",
                "Only native reconciled ortholog pairs are submitted; no clique expansion or graph-edge substitution.",
                "All inference accessions remain in the coverage denominator; no mapping loss accepted.",
                "Empty predictions remain unscored; no fabricated endpoint values.",
                "Conversion is separate from shared-host inference; no isolated speed claim."])
    except BaseException as error:
        report.update(status="native12_composed_qfo_conversion_failed_retained",
            error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report["conversion_finished_monotonic_ns"] = time.monotonic_ns()
        save(destination / "results.json", report)
    return record(destination / "results.json")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("request", "terminal-review", "review-held", "review-release"):
        parser.add_argument("--" + name, type=Path, required=True)
        parser.add_argument("--" + name + "-sha256", required=True)
    parser.add_argument("--review-job", type=int, required=True)
    parser.add_argument("--source-sha256", required=True)
    args = parser.parse_args()
    refs = []
    for name in ("request", "terminal_review", "review_held", "review_release"):
        ref = record(getattr(args, name))
        require(ref["sha256"] == getattr(args, name + "_sha256"), "Native conversion argument digest differs")
        refs.append(ref)
    print(json.dumps(prepare(refs[0], refs[1], args.review_job, refs[2], refs[3],
        DESTINATION, args.source_sha256), sort_keys=True))
