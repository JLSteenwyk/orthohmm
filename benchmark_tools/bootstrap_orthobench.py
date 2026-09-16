#!/usr/bin/env python3
"""Paired RefOG bootstrap, recomputing the weighted OrthoBench statistic."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.orthobench_stage_diagnostics import file_provenance, read_clusters
from benchmark_tools.score_orthobench_partition import read_named_groups, score_partition


METRICS = ("f_score", "precision", "recall")


def statistics(counts):
    """Last axis is weighted TP, FP, FN; outputs are percentages."""
    tp, fp, fn = np.moveaxis(np.asarray(counts, dtype=float), -1, 0)
    def ratio(numerator, denominator):
        return np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator > 0)
    return np.stack((100 * ratio(2 * tp, 2 * tp + fp + fn),
                     100 * ratio(tp, tp + fp), 100 * ratio(tp, tp + fn)), axis=-1)


def weighted_records(records):
    names = [r["refog"] for r in records]
    if not records or len(names) != len(set(names)):
        raise ValueError("RefOG records must be nonempty and unique")
    by_name = {r["refog"]: r for r in records}
    ordered = [by_name[name] for name in sorted(names)]
    sizes = np.array([r["genes"] for r in ordered], dtype=float)
    counts = np.array([[r[k] for k in ("true_positive", "false_positive", "false_negative")]
                       for r in ordered], dtype=float)
    if not np.isfinite(counts).all() or not np.isfinite(sizes).all():
        raise ValueError("Nonfinite sufficient statistics")
    if (sizes < 2).any() or (sizes != np.floor(sizes)).any() or (counts < 0).any():
        raise ValueError("Invalid family sizes or counts")
    if not np.allclose(counts[:, 0] + counts[:, 2], sizes * (sizes - 1) / 2):
        raise ValueError("TP + FN must equal all reference pairs")
    return sorted(names), sizes, counts / (sizes[:, None] - 1)


def paired_bootstrap(records_by_method, baseline, replicates=20000, seed=20260916, alpha=0.05, *, multiplicity_endpoints=None):
    if baseline not in records_by_method or len(records_by_method) < 2:
        raise ValueError("Need a baseline and at least one comparator")
    if replicates < 100 or not 0 < alpha < 1:
        raise ValueError("Need at least 100 replicates and 0 < alpha < 1")
    names, sizes, _ = weighted_records(records_by_method[baseline])
    weights = {}
    for method, records in sorted(records_by_method.items()):
        other_names, other_sizes, weights[method] = weighted_records(records)
        if names != other_names or not np.array_equal(sizes, other_sizes):
            raise ValueError("Methods must have identical RefOG names and sizes")
    rng = np.random.default_rng(seed)
    # All methods share the same multiplicities, preserving paired comparisons.
    multiplicities = rng.multinomial(len(names), np.full(len(names), 1 / len(names)), size=replicates)
    draws = {method: statistics(multiplicities @ counts) for method, counts in weights.items()}
    observed = {method: statistics(counts.sum(axis=0)) for method, counts in weights.items()}
    contrasts = len(weights) - 1
    endpoints = contrasts * len(METRICS) if multiplicity_endpoints is None else multiplicity_endpoints
    if type(endpoints) is not int or endpoints < contrasts * len(METRICS):
        raise ValueError("Multiplicity endpoints must cover every reported contrast/metric")
    comparisons = {}
    for method in weights:
        if method == baseline:
            continue
        differences = draws[method] - draws[baseline]
        individual = statistics(weights[method])[:, 0] - statistics(weights[baseline])[:, 0]
        metrics = {}
        for index, metric in enumerate(METRICS):
            metrics[metric] = {
                "difference_percentage_points": float(observed[method][index] - observed[baseline][index]),
                "paired_percentile_ci": np.quantile(differences[:, index], [alpha / 2, 1 - alpha / 2]).tolist(),
                "bonferroni_percentile_ci": np.quantile(
                    differences[:, index], [alpha / (2 * endpoints), 1 - alpha / (2 * endpoints)]
                ).tolist(),
            }
        comparisons[method] = {
            "versus": baseline, "metrics": metrics,
            "family_f1_wins": int(np.sum(individual > 1e-10)),
            "family_f1_ties": int(np.sum(np.abs(individual) <= 1e-10)),
            "family_f1_losses": int(np.sum(individual < -1e-10)),
        }
    return {
        "baseline": baseline, "families": names, "replicates": replicates,
        "seed": seed, "alpha": alpha, "rng": "numpy.default_rng PCG64 multinomial",
        "numpy_version": np.__version__,
        "point_estimates_percent": {m: dict(zip(METRICS, v.tolist())) for m, v in observed.items()},
        "comparisons": comparisons,
        "multiplicity": f"Bonferroni tail adjustment over {endpoints} " +
                        ("reported contrasts/metrics" if multiplicity_endpoints is None else "planned endpoints"),
        "limitations": [
            "Development-exposed benchmark; intervals do not correct for previous method selection.",
            "RefOG resampling assumes exchangeable families; shared histories and fused predictions can violate independence.",
            "Percentile intervals are approximate, not a guarantee of simultaneous coverage.",
            "Family F1 wins are descriptive; the benchmark is not a mean of family F1 values.",
            "No gene-pair independence assumption and no bootstrap-derived p-values are used.",
        ],
    }


def render_report(result):
    rows = [
        "# OrthoBench Paired Uncertainty", "",
        "Development-exposed analysis; not independent confirmation or a superiority claim.", "",
        f"Baseline: `{result['baseline']}`. {len(result['families'])} RefOGs; "
        f"{result['replicates']:,} paired bootstrap replicates; seed {result['seed']}.", "",
        "| Method | F1 (%) | Precision (%) | Recall (%) |",
        "| --- | ---: | ---: | ---: |",
    ]
    for method, metrics in result["point_estimates_percent"].items():
        rows.append("| " + method + " | " + " | ".join(f"{metrics[k]:.6f}" for k in METRICS) + " |")
    rows.extend(["", "All differences below are method minus baseline, in percentage points.", "",
                 "| Method | Metric | Difference | Paired 95% CI | Multiplicity-adjusted CI |",
                 "| --- | --- | ---: | --- | --- |"])
    for method, comparison in result["comparisons"].items():
        for metric, values in comparison["metrics"].items():
            intervals = [", ".join(f"{v:.3f}" for v in values[key]) for key in
                         ("paired_percentile_ci", "bonferroni_percentile_ci")]
            rows.append(f"| {method} | {metric} | {values['difference_percentage_points']:.3f} | "
                        f"[{intervals[0]}] | [{intervals[1]}] |")
    rows.extend(["", result["multiplicity"] + ".", "",
                 "| Method | Family F1 wins | Ties | Losses |", "| --- | ---: | ---: | ---: |"])
    for method, values in result["comparisons"].items():
        rows.append(f"| {method} | {values['family_f1_wins']} | {values['family_f1_ties']} | {values['family_f1_losses']} |")
    rows.extend(["", "## Limitations", "", *["- " + note for note in result["limitations"]]])
    return "\n".join(rows) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", action="append", required=True, help="NAME=PATH")
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--refogs", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=20260916)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args()
    paths = {}
    for item in args.predictions:
        name, path = item.split("=", 1)
        if not name or name in paths:
            raise ValueError("Prediction names must be nonempty and unique")
        paths[name] = Path(path)
    names = sorted(p.name for p in args.refogs.glob("RefOG*.txt"))
    references = read_named_groups(args.refogs, names)
    uncertain_paths = [args.refogs / "low_certainty_assignments" / n for n in names]
    uncertain = {p.name: set(p.read_text().split()) if p.exists() else set() for p in uncertain_paths}
    scores = {}
    for method, path in paths.items():
        predictions = read_clusters(path)
        seen = set()
        for group in predictions:
            if seen.intersection(group):
                raise ValueError(f"Overlapping predicted groups for {method}")
            seen.update(group)
        scores[method] = score_partition(predictions, references, uncertain)
    result = paired_bootstrap({m: s["refog_records"] for m, s in scores.items()}, args.baseline,
                              args.replicates, args.seed)
    result.update(
        schema_version=1, generated_at=datetime.now(timezone.utc).isoformat(),
        command=[sys.executable, *sys.argv], source=file_provenance(Path(__file__)),
        scorer_source=file_provenance(Path(__file__).with_name("score_orthobench_partition.py")),
        scores=scores,
        inputs={"predictions": {m: file_provenance(p) for m, p in paths.items()},
                "references": [file_provenance(args.refogs / n) for n in names],
                "uncertain": [file_provenance(p) for p in uncertain_paths if p.exists()]},
    )
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(render_report(result))
    print(json.dumps(result["comparisons"], indent=2))


if __name__ == "__main__":
    main()
