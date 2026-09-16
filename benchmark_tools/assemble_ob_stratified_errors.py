"""Evaluate the frozen 84-endpoint, development-exposed OrthoBench strata."""

import argparse
import json
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
import numpy as np
from benchmark_tools.assemble_orthobench_factorial import REFERENCE_SNAPSHOT, compare_official, load_reference_snapshot
from benchmark_tools.audit_historical_profile_ablation import verify_file
from benchmark_tools.bootstrap_orthobench import METRICS, paired_bootstrap, statistics, weighted_records
from benchmark_tools.orthobench_stage_diagnostics import file_provenance, run_official_benchmark
from benchmark_tools.prepare_ob_error_strata import CATEGORIES, assign_strata, prepare as recheck_strata
from benchmark_tools.recover_factorial_postflight import PREPARED_HASH
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.score_orthobench_partition import score_partition
from benchmark_tools.score_ygob_groups import membership, read_predictions

STRATA_SHA = "60208446dca4e51e720f4be0dac5c1ee20a49e72920a13389dd9fbb5d76516a9"
BASELINE = "orthofinder_3_1_5_full"
COMPARATORS = ("orthohmm_high_sensitivity", "orthohmm_phylogeny_satellite_v2")
METHODS = (BASELINE, *COMPARATORS)


def validate_strata(report):
    if (report["status"] != "family_strata_prepared_unscored" or report["accuracy_evaluated"] is not False
            or len(report["families"]) != 70 or report["multiplicity_endpoints"] != 84):
        raise ValueError("Incomplete or changed frozen strata")
    recomputed = assign_strata(report["families"])
    if any(report[key] != value for key, value in recomputed.items()):
        raise ValueError("Strata or illustrative selections differ from input features")


def stratified_statistics(scores, strata):
    if set(scores) != set(METHODS) or set(strata) != {d + ":" + c for d, cats in CATEGORIES.items() for c in cats}:
        raise ValueError("Require all three methods and fourteen planned strata")
    expected_names, expected_sizes, _ = weighted_records(scores[BASELINE]["refog_records"])
    indexed = {}
    for method in METHODS:
        records = scores[method]["refog_records"]
        names, sizes, _ = weighted_records(records)
        if names != expected_names or not np.array_equal(sizes, expected_sizes):
            raise ValueError("Method reference families or sizes differ")
        indexed[method] = {row["refog"]: row for row in records}
    for dimension, categories in CATEGORIES.items():
        combined = [name for category in categories for name in strata[dimension + ":" + category]]
        if sorted(combined) != expected_names:
            raise ValueError("Each feature dimension must partition every reference family exactly once")
    results = {}
    for label, names in sorted(strata.items()):
        selected = {method: [indexed[method][name] for name in sorted(names)] for method in METHODS}
        row = {"families": sorted(names), "family_count": len(names), "multiplicity_endpoints": 84}
        if len(names) >= 5:
            result = paired_bootstrap(selected, BASELINE, replicates=20000, seed=20260918, multiplicity_endpoints=84)
            row.update(status="paired_bootstrap", point_estimates_percent=result["point_estimates_percent"],
                       comparisons=result["comparisons"], bootstrap=result)
        else:
            weights = {method: weighted_records(records)[2] for method, records in selected.items()} if names else {}
            observed = {method: statistics(value.sum(axis=0)) for method, value in weights.items()}
            comparisons = {}
            for method in COMPARATORS:
                individual = statistics(weights[method])[:, 0] - statistics(weights[BASELINE])[:, 0] if names else np.array([])
                comparisons[method] = {"versus": BASELINE, "family_f1_wins": int(np.sum(individual > 1e-10)),
                    "family_f1_ties": int(np.sum(np.abs(individual) <= 1e-10)), "family_f1_losses": int(np.sum(individual < -1e-10)),
                    "metrics": {metric: {"difference_percentage_points": float(observed[method][i] - observed[BASELINE][i]) if names else None,
                        "paired_percentile_ci": None, "bonferroni_percentile_ci": None} for i, metric in enumerate(METRICS)}}
            row.update(status="descriptive_only_lt5" if names else "empty_nonestimable", comparisons=comparisons,
                       point_estimates_percent={method: dict(zip(METRICS, observed[method].tolist())) if names else None for method in METHODS})
        row["weighted_counts"] = {method: dict(zip(("tp", "fp", "fn"), weighted_records(records)[2].sum(axis=0).tolist())) if names else None
                                   for method, records in selected.items()}
        results[label] = row
    return results


def parse_groups(path, method):
    if method not in METHODS:
        raise ValueError("Unknown method format")
    if method == BASELINE:
        groups = read_predictions(path, "named_groups")
    elif method == COMPARATORS[1]:
        groups = read_predictions(path, "root_hogs")
    else:
        with path.open() as handle:
            groups = {str(i): line.split() for i, line in enumerate(handle) if line.strip()}
    membership(groups)
    return groups


def render(result):
    lines = ["# Stratified OrthoBench Error Analysis", "",
             "Exploratory, development-exposed comparisons; descriptors do not establish biological mechanisms.", "",
             "Both OrthoHMM configurations minus full OrthoFinder3.1.5. Differences in percentage points.",
             "20,000 paired RefOG resamples; seed20260918; adjustment retains all84 planned endpoints.",
             "Bins with fewer than five families have no bootstrap intervals. Empty bins are not scored as zero.", "",
             "| Stratum | Families | Method | F1 (%) | F1 Difference | Nominal 95% CI | Adjusted CI | Status |",
             "| --- | ---: | --- | ---: | ---: | --- | --- | --- |"]
    def number(value):
        return "NA" if value is None else f"{value:.3f}"
    def interval(value):
        return "NA" if value is None else "[" + ", ".join(number(v) for v in value) + "]"
    for label, row in result["strata"].items():
        for method in COMPARATORS:
            estimate = row["point_estimates_percent"][method]
            value = row["comparisons"][method]["metrics"]["f_score"]
            lines.append(f"| {label} | {row['family_count']} | {method} | {number(estimate['f_score'] if estimate else None)} | "
                         f"{number(value['difference_percentage_points'])} | {interval(value['paired_percentile_ci'])} | "
                         f"{interval(value['bonferroni_percentile_ci'])} | {row['status']} |")
    lines.extend(["", "## All Endpoint Effects", "", "| Stratum | Method | Metric | Difference | Nominal 95% CI | Adjusted CI |",
                  "| --- | --- | --- | ---: | --- | --- |"])
    for label, row in result["strata"].items():
        for method, comparison in row["comparisons"].items():
            for metric, value in comparison["metrics"].items():
                lines.append(f"| {label} | {method} | {metric} | {number(value['difference_percentage_points'])} | "
                             f"{interval(value['paired_percentile_ci'])} | {interval(value['bonferroni_percentile_ci'])} |")
    lines.extend(["", "## Limitations", "", *["- " + note for note in result["limitations"]]])
    return "\n".join(lines) + "\n"


def assemble(root, output):
    root, output = root.resolve(), output.resolve()
    if output.exists():
        raise FileExistsError(output)
    results = root / "benchmark_tools/results"
    strata_path = results / "ob_error_strata_prepared_20260916.json"
    frozen = read_frozen(strata_path, STRATA_SHA)
    validate_strata(frozen)
    fresh = recheck_strata(root, output / "strata_admission")
    validate_strata(fresh)
    for key in ("families", "strata", "assignments", "identity_median", "illustrative_families_by_stratum"):
        if frozen[key] != fresh[key]:
            raise ValueError("Fresh feature admission differs from frozen strata")
    snapshot_path = results / "orthobench_paired_uncertainty_20260916.json"
    snapshot = read_frozen(snapshot_path, REFERENCE_SNAPSHOT)
    if set(snapshot["scores"]) != set(METHODS) or set(snapshot["inputs"]["predictions"]) != set(METHODS):
        raise ValueError("Historical score snapshot has different methods")
    prepared = read_frozen(results / "orthobench_factorial_prepared_20260916.json", PREPARED_HASH)
    universe = set()
    for item in prepared["fasta_inputs"]:
        verify_file(Path(item["path"]), item)
        for protein in SeqIO.parse(item["path"], "fasta"):
            if protein.id in universe:
                raise ValueError("Duplicate input protein")
            universe.add(protein.id)
    references, uncertain, official, records = load_reference_snapshot(snapshot_path)
    predictions, coverage, scores, official_scores = {}, {}, {}, {}
    sources = snapshot["inputs"]["predictions"]
    for method in METHODS:
        path = Path(sources[method]["path"])
        verify_file(path, sources[method])
        groups = parse_groups(path, method)
        seen = set(membership(groups))
        if not seen <= universe:
            raise ValueError("Native prediction has genes outside the frozen input universe")
        predictions[method] = seen
        coverage[method] = {"input_genes": len(universe), "observed_prediction_genes": len(seen),
                            "missing_input_genes": len(universe - seen), "native_groups": len(groups),
                            "note": "Native group coverage; no unassigned genes silently appended to predictions."}
        clusters = [set(group) for group in groups.values()]
        scores[method] = score_partition(clusters, references, uncertain)
        if scores[method]["refog_records"] != snapshot["scores"][method]["refog_records"]:
            raise ValueError("Fresh full-reference sufficient statistics differ from frozen benchmark")
        converted = output / (method + ".txt")
        with converted.open("x") as handle:
            for group in clusters:
                handle.write(" ".join(sorted(group)) + "\n")
        official_scores[method] = run_official_benchmark(official, converted)
        compare_official(scores[method], official_scores[method])
    analysis = stratified_statistics(scores, frozen["strata"])
    for label, row in analysis.items():
        genes = set().union(*(references[name] for name in row["families"]))
        row["coverage"] = {method: {"reference_genes": len(genes), "observed_prediction_genes": len(genes & predictions[method])}
                           for method in METHODS}
    for item in [*sources.values(), *records, *prepared["fasta_inputs"]]:
        verify_file(Path(item["path"]), item)
    result = {"status": "stratified_analysis_complete", "publication_ready": False, "baseline": BASELINE,
              "job_id": os.environ.get("SLURM_JOB_ID"), "command": [sys.executable, *sys.argv],
              "multiplicity_endpoints": 84, "replicates": 20000, "seed": 20260918, "alpha": .05,
              "strata": analysis, "scores": scores, "official_scores": official_scores, "native_coverage": coverage,
              "feature_admission": fresh, "frozen_strata": file_provenance(strata_path),
              "frozen_scores": file_provenance(snapshot_path), "predictions": sources, "references": records,
              "source": file_provenance(Path(__file__)), "bootstrap_source": file_provenance(Path(__file__).with_name("bootstrap_orthobench.py")),
              "official_scorer": file_provenance(official),
              "limitations": [
                  "Post-development exploratory analysis; no independent confirmation, equivalence claim or new default.",
                  "Strata overlap and share histories; within-stratum RefOG resampling assumes exchangeable families, not independent gene pairs.",
                  "All84 endpoints remain in multiplicity adjustment, including empty and descriptive-only bins.",
                  "Adjusted percentile tails have about six draws per tail with20,000 replicates; uncertainty is approximate.",
                  "Copy number, identity, relative length and global composition do not establish duplication history, fragments, domains or causal mechanisms.",
                  "Sufficient statistics retain full-reference scoring and low-certainty conventions before restriction to each stratum.",
                  "This analysis uses previously audited retained group outputs, not new end-to-end inference or QfO evidence."]}
    (output / "results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    (output / "results.md").write_text(render(result))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    assemble(args.root, args.output)


if __name__ == "__main__":
    main()
