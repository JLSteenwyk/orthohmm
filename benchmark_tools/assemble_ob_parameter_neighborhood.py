"""Score all six prespecified parameter variants after fresh native admission."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
from benchmark_tools.assemble_orthobench_factorial import compare_official, coverage_and_resources, load_reference_snapshot
from benchmark_tools.audit_historical_profile_ablation import verify_file
from benchmark_tools.bootstrap_orthobench import render_report
from benchmark_tools.orthobench_stage_diagnostics import file_provenance, run_official_benchmark
from benchmark_tools.parameter_neighborhood_statistics import BASELINE, VARIANTS, summarize
from benchmark_tools.recover_factorial_postflight import PREPARED_HASH
from benchmark_tools.run_ob_candidate_neighborhood import VARIANTS as THRESHOLD_VARIANTS, CPM_VARIANTS
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.score_orthobench_partition import score_partition
from benchmark_tools.score_ygob_groups import read_predictions
from benchmark_tools.validate_factorial_native import validate as validate_baseline
from benchmark_tools.validate_ob_candidate_phylogeny import validate as validate_variant


def admit_panel(root):
    baseline = validate_baseline(root, 3)
    variants = {}
    for labels, cpm in ((CPM_VARIANTS, True), (THRESHOLD_VARIANTS, False)):
        for index, label in enumerate(labels):
            variants[label] = validate_variant(root, index, cpm=cpm)
    return {"control": baseline, "variants": variants}


def native_methods(panel):
    if set(panel["variants"]) != set(VARIANTS):
        raise ValueError("Require exactly all six parameter variants")
    for label, row in panel["variants"].items():
        if (row["status"] != "candidate_variant_native_validated_unscored"
                or row["accuracy_evaluated"] is not False or row["variant"] != label):
            raise ValueError("Invalid or mismatched variant admission")
    methods = {BASELINE: panel["control"],
               **{label: panel["variants"][label]["native_validation"] for label in VARIANTS}}
    if any(row["status"] != "native_group_output_verified" or row["native_outputs_validated"] is not True
           or row["accuracy_evaluated"] is not False for row in methods.values()):
        raise ValueError("Require successful unscored native admission for every output")
    return methods


def assemble(root, output):
    root, output = root.resolve(), output.resolve()
    if output.exists():
        raise FileExistsError(output)
    helpers = [file_provenance(Path(__file__).with_name(name)) for name in (
        "parameter_neighborhood_statistics.py", "bootstrap_orthobench.py", "score_orthobench_partition.py",
        "assemble_orthobench_factorial.py", "validate_ob_candidate_phylogeny.py", "validate_factorial_native.py")]
    # Finish label-blind admission for the entire panel before reading references.
    panel = admit_panel(root)
    methods = native_methods(panel)
    results = root / "benchmark_tools/results"
    prepared_path = results / "orthobench_factorial_prepared_20260916.json"
    prepared = read_frozen(prepared_path, PREPARED_HASH)
    gene_species = {}
    for item in prepared["fasta_inputs"]:
        verify_file(Path(item["path"]), item)
        for sequence in SeqIO.parse(item["path"], "fasta"):
            if sequence.id in gene_species:
                raise ValueError("Duplicate FASTA ID")
            gene_species[sequence.id] = Path(item["path"]).name
    predictions, sources, coverage, tracked = {}, {}, {}, []
    for label, admission in methods.items():
        metrics_record = admission["native_metrics"]
        verify_file(Path(metrics_record["path"]), metrics_record)
        metrics = json.loads(Path(metrics_record["path"]).read_text())
        source = metrics["outputs"]["root_hogs"]
        path = Path(source["path"])
        verify_file(path, source)
        groups = [set(group) for group in read_predictions(path, "root_hogs").values()]
        predictions[label], sources[label] = groups, source
        coverage[label] = coverage_and_resources(groups, gene_species, metrics)
        tracked.extend([metrics_record, source])
    references, uncertain, official, records = load_reference_snapshot(
        results / "orthobench_paired_uncertainty_20260916.json")
    if not set().union(*references.values()).issubset(gene_species):
        raise ValueError("Reference genes absent from inference universe")
    output.mkdir(parents=True)
    scores, official_scores, conversions = {}, {}, {}
    for label, groups in predictions.items():
        converted = output / f"{label}.txt"
        with converted.open("x") as handle:
            for group in groups:
                handle.write(" ".join(sorted(group)) + "\n")
        conversions[label] = file_provenance(converted)
        scores[label] = score_partition(groups, references, uncertain)
        official_scores[label] = run_official_benchmark(official, converted)
        compare_official(scores[label], official_scores[label])
    result = summarize({label: score["refog_records"] for label, score in scores.items()}, {})
    for item in [*records, *tracked, *prepared["fasta_inputs"], *conversions.values(), *helpers]:
        verify_file(Path(item["path"]), item)
    result.update(schema_version=1, scores=scores, official_scores=official_scores,
                  native_validation=panel, predictions=sources, references=records,
                  conversions=conversions, coverage_resources=coverage, helpers=helpers,
                  official_scorer=file_provenance(official), assembler=file_provenance(Path(__file__)),
                  analysis_scope="Post-development six-variant parameter sensitivity; no default selection",
                  timing_scope="Incremental inferred-phylogeny runs with exact-input checkpoint reuse on a shared node; not end-to-end")
    result["limitations"].append(
        "CPM variants recompute their own HMM seeds and candidates; species trees are inferred under unchanged native rules. "
        "Reported phylogeny costs omit upstream replay and candidate construction and are not matched timing comparisons.")
    with (output / "results.json").open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    with (output / "results.md").open("x") as handle:
        handle.write("# OrthoBench Parameter Neighborhood\n\n" + result["analysis_scope"] + ".\n\n"
                     + render_report(result) + "\n\n" + result["timing_scope"] + ".\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    assemble(args.root, args.output)
