"""Score the terminal, native-validated unconstrained OrthoBench diagnostic."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
from benchmark_tools.assemble_orthobench_factorial import (
    compare_official, coverage_and_resources, load_reference_snapshot,
)
from benchmark_tools.bootstrap_orthobench import paired_bootstrap, render_report
from benchmark_tools.orthobench_stage_diagnostics import file_provenance, run_official_benchmark
from benchmark_tools.recover_factorial_postflight import PREPARED_HASH
from benchmark_tools.run_orthobench_factorial_cell import select_cell, unconstrained_cell
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.score_orthobench_partition import score_partition
from benchmark_tools.score_ygob_groups import read_predictions
from benchmark_tools.validate_factorial_native import validate as validate_constrained
from benchmark_tools.validate_unconstrained_control import validate as validate_unconstrained

BASELINE = "p1_c1_r1"
DIAGNOSTIC = "p1_c1_r1_unconstrained_v2"


def bootstrap(scores):
    if set(scores) != {BASELINE, DIAGNOSTIC}:
        raise ValueError("Require exactly the constrained and unconstrained diagnostic pair")
    return paired_bootstrap({name: value["refog_records"] for name, value in scores.items()},
                            BASELINE, replicates=20000, seed=20260918)


def assemble(root, output):
    root, output = root.resolve(), output.resolve()
    if output.exists():
        raise FileExistsError(output)
    # Gate execution before loading predictions or benchmark reference labels.
    native = {DIAGNOSTIC: validate_unconstrained(root), BASELINE: validate_constrained(root, 3)}
    results = root / "benchmark_tools/results"
    prepared_path = results / "orthobench_factorial_prepared_20260916.json"
    prepared = read_frozen(prepared_path, PREPARED_HASH)
    original, directory, _ = select_cell(prepared, 3)
    diagnostic = unconstrained_cell(original, directory)
    gene_species = {}
    for item in prepared["fasta_inputs"]:
        for record in SeqIO.parse(item["path"], "fasta"):
            if record.id in gene_species:
                raise ValueError("Duplicate FASTA ID")
            gene_species[record.id] = Path(item["path"]).name
    predictions, sources, coverage = {}, {}, {}
    for cell in (original, diagnostic):
        label, path = cell["label"], Path(cell["prediction"])
        sources[label] = file_provenance(path)
        groups = read_predictions(path, "root_hogs")
        predictions[label] = [set(group) for group in groups.values()]
        metrics = json.loads(Path(native[label]["native_metrics"]["path"]).read_text())
        coverage[label] = coverage_and_resources(predictions[label], gene_species, metrics)
    references, uncertain, official, records = load_reference_snapshot(
        results / "orthobench_paired_uncertainty_20260916.json")
    if not set().union(*references.values()).issubset(gene_species):
        raise ValueError("Reference genes absent from inference universe")
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
    for record in [*records, *sources.values(), *prepared["fasta_inputs"]]:
        if file_provenance(Path(record["path"])) != record:
            raise ValueError("Scoring input changed during evaluation")
    result = bootstrap(scores)
    result.update(schema_version=1, scores=scores, official_scores=official_scores,
                  native_validation=native, predictions=sources, references=records,
                  coverage_resources=coverage, official_scorer=file_provenance(official),
                  assembler=file_provenance(Path(__file__)), publication_ready=False,
                  analysis_scope="Exploratory unconstrained diagnostic; not a ninth factorial cell",
                  timing_scope="Incremental cached reconciliation on a shared node")
    result["limitations"].append(
        "Detailed execution specified after factorial outcomes; not independent confirmation or a new default.")
    with (output / "results.json").open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    with (output / "results.md").open("x") as handle:
        handle.write("# Unconstrained Reconciliation Diagnostic\n\n" + result["analysis_scope"] + ".\n\n")
        handle.write(render_report(result))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    assemble(args.root, args.output)


if __name__ == "__main__":
    main()
