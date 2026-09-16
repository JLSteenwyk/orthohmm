"""Gated OrthoBench scoring of two sequence-search controls against p0_c0_r0."""

import argparse
import json
import math
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
from benchmark_tools.assemble_orthobench_factorial import compare_official, coverage_and_resources, load_reference_snapshot
from benchmark_tools.audit_historical_profile_ablation import read_partition, verify_file
from benchmark_tools.bootstrap_orthobench import paired_bootstrap, render_report
from benchmark_tools.orthobench_stage_diagnostics import file_provenance, run_official_benchmark
from benchmark_tools.recover_factorial_postflight import PREPARED_HASH
from benchmark_tools.run_orthobench_factorial_cell import select_cell, verify_prepared
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.score_orthobench_partition import score_partition
from benchmark_tools.validate_sequence_graph_control import validate, VARIANTS
from benchmark_tools.verify_ygob_validation import require_completed_job

BASELINE = "p0_c0_r0"


def bootstrap(scores):
    if set(scores) != {BASELINE, *VARIANTS}:
        raise ValueError("Require the frozen HMM baseline and both sequence controls")
    return paired_bootstrap({name: value["refog_records"] for name, value in scores.items()},
                            BASELINE, replicates=20000, seed=20260918)


def graph_resources(native):
    resources = {BASELINE: {"wall_s": None, "scope": "Historical cached baseline; no separately matched graph timing"}}
    for label in VARIANTS:
        record = native["variants"][label]["native_records"]["replay.json"]
        verify_file(Path(record["path"]), record)
        metrics = json.loads(Path(record["path"]).read_text())
        values = {key: metrics[key] for key in ("wall_s", "peak_process_rss_gib")}
        if any(isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) or v < 0 for v in values.values()):
            raise ValueError("Invalid graph resource measurement")
        resources[label] = {**values, "timings": metrics["timings"], "source": record,
            "scope": "Incremental shared-node graph replay; search/conversion excluded; peak RSS is replay process only, not process tree"}
    return resources


def admit_coverage(root, native):
    accounting = subprocess.check_output(["sacct", "-j", "21294", "--parsable2",
                                         "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = require_completed_job(accounting, 21294)
    path = root / "benchmarks/results/ob_search_hit_coverage_v1.json"
    record = file_provenance(path)
    report = json.loads(path.read_text())
    if report["accuracy_evaluated"] is not False or set(report["searches"]) != {"hmm", *VARIANTS}:
        raise ValueError("Wrong coverage inventory or prior accuracy use")
    if report["conversion"] != native["conversion"] or report["scheduler"]["JobIDRaw"] != "21292":
        raise ValueError("Coverage used a different conversion")
    source = root / "benchmarks/work/publication_ob_hit_coverage_v1/benchmark_tools/compare_search_hit_coverage.py"
    if report["source"] != file_provenance(source):
        raise ValueError("Coverage source changed")
    for item in (report["conversion"], report["hmm_cache"]):
        verify_file(Path(item["path"]), item)
    for label in VARIANTS:
        for key in ("manifest", "summary"):
            if report["numeric_admissions"][label][key] != native["variants"][label]["numeric_admission"][key]:
                raise ValueError("Coverage and graph checkpoint evidence differ")
    return report, {"report": record, "scheduler": scheduler}


def assemble(root, output):
    root, output = root.resolve(), output.resolve()
    if output.exists():
        raise FileExistsError(output)
    native = validate(root)
    hit_coverage, coverage_gate = admit_coverage(root, native)
    resources = graph_resources(native)
    results = root / "benchmark_tools/results"
    prepared = read_frozen(results / "orthobench_factorial_prepared_20260916.json", PREPARED_HASH)
    cell, _, launcher = select_cell(prepared, 0)
    verify_prepared(prepared, cell, launcher, 21161)
    if hit_coverage["hmm_cache"] != prepared["cache"]:
        raise ValueError("Coverage HMM cache differs from baseline cache")
    gene_species = {}
    for item in prepared["fasta_inputs"]:
        for record in SeqIO.parse(item["path"], "fasta"):
            if record.id in gene_species:
                raise ValueError("Duplicate FASTA gene ID")
            gene_species[record.id] = Path(item["path"]).name
    sources = {BASELINE: prepared["candidate_arms"]["p0_c0"]["candidate_partition"],
               **{label: native["variants"][label]["prediction"] for label in VARIANTS}}
    predictions, coverage = {}, {}
    for label, source in sources.items():
        path = Path(source["path"])
        verify_file(path, source)
        predictions[label] = read_partition(path, set(gene_species))
        coverage[label] = coverage_and_resources(predictions[label], gene_species)
    references, uncertain, official, records = load_reference_snapshot(
        results / "orthobench_paired_uncertainty_20260916.json")
    if not set().union(*references.values()).issubset(gene_species):
        raise ValueError("Reference genes absent from input universe")
    output.mkdir(parents=True)
    scores, official_scores = {}, {}
    for label, groups in predictions.items():
        converted = output / f"{label}.txt"
        with converted.open("x") as handle:
            for group in groups:
                handle.write(" ".join(sorted(group)) + "\n")
        scores[label] = score_partition(groups, references, uncertain)
        official_scores[label] = run_official_benchmark(official, converted)
        compare_official(scores[label], official_scores[label])
    for record in [*records, *sources.values(), *prepared["fasta_inputs"], coverage_gate["report"]]:
        verify_file(Path(record["path"]), record)
    result = bootstrap(scores)
    result.update(schema_version=1, scores=scores, official_scores=official_scores,
                  native_validation=native, hit_coverage=hit_coverage, coverage_admission=coverage_gate,
                  predictions=sources, references=records, gene_coverage=coverage,
                  graph_resources=resources,
                  official_scorer=file_provenance(official), assembler=file_provenance(Path(__file__)),
                  publication_ready=False, analysis_scope="Exploratory initial-search replacement; profiles, candidate expansion and reconciliation off")
    result["limitations"].extend([
        "All-hit DIAMOND is the primary control; post-search top100 is a reporting diagnostic, not HMM prefilter emulation.",
        "Equal E-value cutoffs and length divisors do not establish matched sensitivity or score calibration.",
        "Costs are shared-node incremental graph replays plus separately recorded search; not matched efficiency evidence.",
        "Control specified after development and YGOB outcomes; not independent confirmation or a new default."])
    with (output / "results.json").open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    with (output / "results.md").open("x") as handle:
        handle.write("# Sequence-Search Control\n\n" + result["analysis_scope"] + ".\n\n" + render_report(result))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    assemble(args.root, args.output)


if __name__ == "__main__":
    main()
