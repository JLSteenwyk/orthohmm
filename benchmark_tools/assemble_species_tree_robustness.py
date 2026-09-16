"""Score the complete, native-validated species-tree perturbation panel."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
from benchmark_tools.assemble_orthobench_factorial import compare_official, coverage_and_resources, load_reference_snapshot
from benchmark_tools.audit_historical_profile_ablation import verify_file
from benchmark_tools.bootstrap_orthobench import paired_bootstrap, render_report
from benchmark_tools.orthobench_stage_diagnostics import file_provenance, run_official_benchmark
from benchmark_tools.recover_factorial_postflight import PREPARED_HASH
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.run_species_tree_control import PERTURBATIONS
from benchmark_tools.score_orthobench_partition import score_partition
from benchmark_tools.score_ygob_groups import read_predictions
from benchmark_tools.validate_species_tree_perturbations import validate

BASELINE = "supplied_control"


def bootstrap(scores):
    if set(scores) != {BASELINE, *PERTURBATIONS}:
        raise ValueError("Require exactly the supplied control and all six fixed perturbations")
    return paired_bootstrap({name: value["refog_records"] for name, value in scores.items()},
                            BASELINE, replicates=20000, seed=20260918)


def native_methods(panel):
    if (panel["status"] != "native_panel_validated_unscored" or panel["accuracy_evaluated"] is not False
            or set(panel["variants"]) != set(PERTURBATIONS)
            or panel["control_validation"]["status"] != "equivalent"):
        raise ValueError("Incomplete panel or failed supplied-control equivalence")
    methods = {BASELINE: panel["control_validation"]["native_validation"],
               **{label: panel["variants"][label]["native_validation"] for label in PERTURBATIONS}}
    if any(item["status"] != "native_group_output_verified" or item["native_outputs_validated"] is not True
           or item["accuracy_evaluated"] is not False for item in methods.values()):
        raise ValueError("Every output must have native admission without accuracy evaluation")
    return methods


def render_resources(result):
    lines = ["", "## Coverage And Incremental Resources", "",
             "All inputs include singletons. Costs reuse raw gene trees on a shared node; they are not end-to-end timings.", "",
             "| Tree | Rooted Clade Distance | Groups | Singleton Groups | Genes In Multispecies Groups | Wall (s) | Peak Tree RSS (GiB) |",
             "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for label in (BASELINE, *PERTURBATIONS):
        row = result["coverage_resources"][label]
        distance = 0 if label == BASELINE else result["native_validation"]["variants"][label]["rooted_rf_clade_distance"]
        lines.append(f"| {label} | {distance} | {row['groups']} | {row['singleton_groups']} | "
                     f"{row['genes_in_multispecies_groups']} | {row['wall_s']:.2f} | {row['peak_process_tree_rss_bytes'] / 2**30:.3f} |")
    return "\n".join(lines) + "\n"


def assemble(root, output):
    root, output = root.resolve(), output.resolve()
    if output.exists():
        raise FileExistsError(output)
    # Complete label-blind admission before accessing reference outcomes.
    panel = validate(root)
    native = native_methods(panel)
    results = root / "benchmark_tools/results"
    prepared_path = results / "orthobench_factorial_prepared_20260916.json"
    prepared = read_frozen(prepared_path, PREPARED_HASH)
    gene_species = {}
    for item in prepared["fasta_inputs"]:
        verify_file(Path(item["path"]), item)
        for record in SeqIO.parse(item["path"], "fasta"):
            if record.id in gene_species:
                raise ValueError("Duplicate FASTA ID")
            gene_species[record.id] = Path(item["path"]).name
    predictions, sources, coverage, tracked = {}, {}, {}, []
    for label, admission in native.items():
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
    scores, official_scores = {}, {}
    for label, groups in predictions.items():
        converted = output / f"{label}.txt"
        with converted.open("x") as handle:
            for group in groups:
                handle.write(" ".join(sorted(group)) + "\n")
        scores[label] = score_partition(groups, references, uncertain)
        official_scores[label] = run_official_benchmark(official, converted)
        compare_official(scores[label], official_scores[label])
    for record in [*records, *tracked, *prepared["fasta_inputs"]]:
        verify_file(Path(record["path"]), record)
    result = bootstrap(scores)
    result.update(schema_version=1, scores=scores, official_scores=official_scores,
                  native_validation=panel, predictions=sources, references=records,
                  coverage_resources=coverage, official_scorer=file_provenance(official),
                  assembler=file_provenance(Path(__file__)), publication_ready=False,
                  analysis_scope="Exploratory prespecified topology stress test; unchanged frozen method and reused raw gene trees",
                  timing_scope="Incremental cached supplied-tree reconciliation on a shared node")
    result["limitations"].extend([
        "All six fixed variants and 18 F1/precision/recall endpoints are retained; no best-tree selection.",
        "The unchanged supplied-tree control exactly reproduces the inferred baseline; this is not independent tree reconstruction.",
        "Rooted NNI variants are controlled perturbations, not posterior draws or an empirical tree-error distribution.",
        "Lengths travel with subtrees; these topology edits do not preserve all evolutionary distances.",
        "Specified after development and YGOB outcomes; not independent confirmation, parameter robustness, or a new default."])
    with (output / "results.json").open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    with (output / "results.md").open("x") as handle:
        handle.write("# Species-Tree Robustness\n\n" + result["analysis_scope"] + ".\n\n" + render_report(result) + render_resources(result))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    assemble(args.root, args.output)


if __name__ == "__main__":
    main()
