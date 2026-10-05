"""Convert terminal-reviewed full-native QfO outputs, never cached admissions.

R-off uses cross-species cluster cliques; R-on uses resolved native pairs.
Preparation does not execute endpoint assessment or authorize another inference.
"""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.collect_factorial_resources import FIXED_INPUTS
from benchmark_tools.derive_threadripper_resources import SCOPES
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.prepare_qfo_corrected_comparator_pairs import ENV_SHA
from benchmark_tools.prepare_qfo_corrected_group_pairs import expected_pairs
from benchmark_tools.prepare_qfo_factorial_pairs import write_native_pairs
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.qfo_filter_pairs import filter_pairs, load_mapping
from benchmark_tools.run_native_factorial_cost import ROOT, SCOPE, IDENTITIES, read, validate_plan, validate_request, verify_terminal
from benchmark_tools.validate_native_factorial_outputs import Evidence, require, universe
from qfo_benchmark.og_to_pairwise import _strip_to_uniprot


def conversion_kind(run):
    require(run.get("dataset") == "qfo_corrected" and type(run.get("index")) is int
        and 6 <= run["index"] < 13, "Require one of the seven full-native corrected-QfO identities")
    require((run["dataset"], run["cell"]) == IDENTITIES[run["index"]], "Native QfO index/cell binding differs")
    return "native" if run["cell"].endswith("r1") else "group"


def admit_conversion(review, request_ref, request, run):
    kind = conversion_kind(run)
    require(review.get("schema") == "native_factorial_terminal_review_v1"
        and review.get("request") == request_ref and review.get("plan") == request["plan"]
        and review.get("job_id") == request["job_id"] and review.get("index") == run["index"]
        and review.get("dataset") == run["dataset"] and review.get("cell") == run["cell"]
        and review.get("repeat") == run["repeat"] and review.get("status") == "native_success"
        and review.get("scheduler_state") == "COMPLETED" and review.get("scheduler_exit_code") == "0:0"
        and all(review.get(k) is True for k in ("terminal_reviewed", "native_outputs_validated",
            "primary_resources_replayed", "shared_host_resources_reviewed"))
        and review.get("execution_scope") == SCOPE and review.get("resource_scopes") == SCOPES
        and review.get("uncontended_timing") is False
        and review.get("automatic_retry") is False
        and set(review.get("reviews", {})) == {"runtime", "resources", "environment", "outputs_or_failure"},
        "Require a bound successful full-native terminal/runtime/resource/environment/output review")
    return kind


def normalize_owners(owners, mapping):
    normalized = {}
    for gene, species in owners.items():
        accession = _strip_to_uniprot(gene)
        require(accession and accession not in normalized, "Noninjective accession normalization")
        require(accession in mapping, "Corrected inference input absent from frozen QfO mapping")
        normalized[accession] = species
    return normalized


def pair_coverage(path, normalized):
    covered, count = set(), 0
    with path.open() as stream:
        for line in stream:
            pair = line.rstrip("\n").split("\t")
            require(len(pair) == 2 and pair[0] < pair[1]
                and all(g in normalized for g in pair), "Invalid normalized prediction pair")
            require(normalized[pair[0]] != normalized[pair[1]], "Converted pair is not cross-species")
            covered.update(pair)
            count += 1
    return dict(pair_rows=count, input_accessions=len(normalized), accessions_in_any_pair=len(covered),
        fraction_inputs_in_any_pair=len(covered) / len(normalized) if normalized else 0.,
        definition="fraction of all inference input accessions present in at least one submitted pair")


def convert(kind, source, destination, input_directory, owners, expected_native_count=None, runner=subprocess.run):
    require(kind in {"group", "native"}, "Unknown conversion semantics")
    if kind == "native":
        require(type(expected_native_count) is int and expected_native_count >= 0, "Invalid admitted native pair count")
        count = write_native_pairs(source, destination, owners, expected_native_count)
        return count, None
    count = expected_pairs(source, owners)
    command = [sys.executable, str(ROOT / "qfo_benchmark/og_to_pairwise.py"), str(source), str(input_directory)]
    with destination.open("x") as stream, destination.with_name("conversion.log").open("x") as log:
        runner(command, stdout=stream, stderr=log, check=True)
    return count, command


def prepare(request_ref, review_ref, destination):
    source_ref = record(__file__)
    request, review = read(request_ref), read(review_ref)
    plan = read(request["plan"])
    validate_request(request, request["plan"], request["job_id"])
    run = validate_plan(plan)[request["index"]]
    kind = admit_conversion(review, request_ref, request, run)
    require(os.environ.get("SLURM_CPUS_PER_TASK") == "2"
        and os.environ.get("SLURM_JOB_ID", "").isdigit(), "Require a scheduled two-CPU conversion")
    destination = Path(destination)
    require(destination.is_absolute() and destination.resolve() == destination
        and destination.is_relative_to(ROOT) and not destination.is_relative_to(Path(plan["panel_root"])),
        "Require separate direct conversion destination")
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    terminal = verify_terminal(request["job_id"])
    fields = terminal["verified"].get("fields", terminal["verified"])
    require(fields.get("JobState", fields.get("State")) == "COMPLETED" and fields.get("ExitCode") == "0:0",
            "Native inference no longer has a successful terminal scheduler outcome")
    if terminal["source"] == "live_controller":
        require(fields.get("Comment") == request_ref["sha256"], "Native scheduler comment differs")
    evidence = Evidence()
    for ref in [request_ref, request["plan"], review_ref, review["source"], review["scheduler"],
                *review["reviews"].values(), plan["baseline"], *plan["helper_sources"]]:
        check(ref)
        evidence.bind(ref["path"])
    require(review["source"] == record(ROOT / "benchmark_tools/review_native_factorial_attempt.py"),
            "Terminal reviewer source differs")
    outputs = read(review["reviews"]["outputs_or_failure"])
    require(outputs.get("native_outputs_validated") is True and outputs.get("request") == request_ref
        and outputs.get("plan") == request["plan"] and outputs.get("index") == run["index"]
        and outputs.get("cell") == run["cell"] and outputs.get("job_id") == request["job_id"],
        "Independent native output binding differs")
    for ref in [*outputs["checked_files"], *outputs["evidence"]]:
        check(ref)
        evidence.bind(ref["path"])
    name, digest = FIXED_INPUTS["qfo_preparation"]
    prepared_ref = record(ROOT / name)
    require(prepared_ref["sha256"] == digest and read(prepared_ref)["input_fastas"] == run["inputs"],
            "Full-native inputs differ from frozen corrected-QfO release")
    evidence.bind(prepared_ref["path"])
    env_ref = record(ROOT / "benchmark_tools/results/qfo_assessment_environment_20260917.json")
    require(env_ref["sha256"] == ENV_SHA, "Frozen QfO reference environment differs")
    environment = read(env_ref)
    evidence.bind(env_ref["path"])
    mappings = [r for r in environment["reference_files"] if Path(r["path"]).name == "mapping.json.gz"]
    require(len(mappings) == 1, "Require one frozen QfO accession mapping")
    mapping_ref = mappings[0]
    check(mapping_ref)
    mapping = load_mapping(evidence.bind(mapping_ref["path"]))
    input_directory = Path(run["output_root"]) / "input"
    owners, owner_digest, per_species = universe(dict(run, input_directory=str(input_directory)), evidence)
    require(owner_digest == outputs["gene_ownership_sha256"]
        and per_species == outputs["per_species_counts"], "Conversion gene ownership differs from independent admission")
    normalized = normalize_owners(owners, mapping)
    suffix = "native/orthohmm_phylogeny/orthohmm_pairwise_orthologs.tsv" if kind == "native" else "native/orthohmm_working_res/orthohmm_edges_clustered.txt"
    prediction = Path(run["output_root"]) / suffix
    checked = {r["path"]: r for r in outputs["checked_files"]}
    require(str(prediction) in checked, "Selected prediction was not independently admitted")
    prediction_ref = checked[str(prediction)]
    check(prediction_ref)
    require({p.name for p in input_directory.glob("*.fasta")} == set(run["native_order"]),
            "QfO converter FASTA glob differs from admitted native input inventory")
    for path in [Path(__file__), ROOT / "qfo_benchmark/og_to_pairwise.py", *[ROOT / "benchmark_tools" / name
        for name in ("prepare_qfo_corrected_group_pairs.py", "prepare_qfo_factorial_pairs.py", "qfo_filter_pairs.py",
                     "simulation_method_outputs.py", "validate_native_factorial_outputs.py")]]:
        evidence.bind(path)
    require(record(__file__) == source_ref, "Conversion source changed")
    destination.mkdir(parents=True, exist_ok=False)
    report = dict(schema="full_native_factorial_qfo_conversion_v1", status="preparing_unscored", source=source_ref,
        job_id=os.environ["SLURM_JOB_ID"], native_job_id=request["job_id"], native_index=run["index"], cell=run["cell"],
        participant="ohmm_qfo_full_native_" + run["cell"], request=request_ref, plan=request["plan"],
        terminal_review=review_ref, native_scheduler=terminal, prepared=prepared_ref, environment_manifest=env_ref,
        input_fastas=run["inputs"], gene_ownership_sha256=owner_digest, mapping=mapping_ref, native_input=prediction_ref,
        semantics="native phylogenetically inferred pairs" if kind == "native" else "cross-species group-derived clique pairs",
        conversion_kind=kind, accuracy_evaluated=False, native_inference_reexecuted=False, automatic_retry=False,
        next_identity_authorized=False, publication_ready=False, conversion_started_monotonic_ns=time.monotonic_ns())
    report["conversion_interval_scope"] = "pair materialization/filter/coverage/postflight; excludes prior admission and ownership indexing"
    save(destination / "preflight.json", report)
    try:
        pairs, filtered = destination / "pairs.partial.tsv", destination / "pairs.qfo.partial.tsv"
        expected = outputs["phylogeny"]["native_pair_rows"] if kind == "native" else None
        count, command = convert(kind, prediction, pairs, input_directory, owners, expected)
        total, retained = filter_pairs(pairs, filtered, mapping)
        require(count == total == retained, "Pair count differs or corrected mapping loses predictions")
        coverage = pair_coverage(filtered, normalized)
        require(coverage["pair_rows"] == count, "Independent converted coverage row count differs")
        require(record(pairs)["sha256"] == record(filtered)["sha256"]
            and pairs.stat().st_size == filtered.stat().st_size, "Mapping changed corrected predictions")
        checked_files = evidence.finish()
        pairs.rename(destination / "pairs.tsv")
        filtered.rename(destination / "pairs.qfo.tsv")
        report.update(status="full_native_factorial_qfo_pairs_prepared_unscored", pairs=record(destination / "pairs.tsv"),
            filtered_pairs=record(destination / "pairs.qfo.tsv"), total_pairs=total, retained_pairs=retained,
            removed_mapping_pairs=0, expected_pairs=count, empty_predictions=count == 0, pair_coverage=coverage,
            command=command, checked_records=checked_files,
            limitations=["Conversion only; native QfO endpoints and independent assessment admission remain required.",
                "R-on resolved native pairs and R-off group cliques have different prediction semantics.",
                "Protein relation coverage is reported separately from pair volume; singleton/no-relation inputs remain in the denominator.",
                "Empty predictions are retained, not assigned fabricated benchmark scores.",
                "Full-native and older cached-replay provenance/namespaces are not interchangeable.",
                "Conversion is outside the measured inference interval; no isolated efficiency or causal resource claim."])
        if (destination / "conversion.log").exists():
            report["conversion_log"] = record(destination / "conversion.log")
    except BaseException as error:
        report.update(status="full_native_factorial_qfo_conversion_failed_retained",
            error_type=type(error).__name__, error=str(error), accuracy_evaluated=False)
        raise
    finally:
        report["conversion_finished_monotonic_ns"] = time.monotonic_ns()
        save(destination / "results.json", report)
    return record(destination / "results.json")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("request", "terminal-review", "output-directory"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("request-sha256", "terminal-review-sha256"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    request_ref, review_ref = record(args.request), record(args.terminal_review)
    require(request_ref["sha256"] == args.request_sha256 and review_ref["sha256"] == args.terminal_review_sha256,
            "Conversion request/review checksum differs")
    print(json.dumps(prepare(request_ref, review_ref, args.output_directory.absolute()), sort_keys=True))
