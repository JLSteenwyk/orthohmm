"""Prespecified paired seed-level analysis; never replace failures with zeros."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import sys

import numpy as np


SEEDS = tuple(range(20261001, 20261011))
CONDITIONS = ("baseline", "divergent", "turnover", "divergent_turnover",
              "missing20", "uneven_taxa", "taxon_count_control")
METHODS = ("orthohmm_high_sensitivity", "orthohmm_satellite_v2", "orthofinder_full", "orthofinder_sequence_only")
METRICS = ("f1", "precision", "recall")


def validate_score(score):
    counts = [score[k] for k in ("tp", "fp", "fn", "input_genes", "eligible_true_pairs")]
    if any(type(n) is not int or n < 0 for n in counts):
        raise ValueError("Nonnegative integer counts required")
    tp, fp, fn = counts[:3]
    if tp + fn != score["eligible_true_pairs"]:
        raise ValueError("True-pair count mismatch")
    denominators = (2 * tp + fp + fn, tp + fp, tp + fn)
    for metric, numerator, denominator in zip(METRICS, (2 * tp, tp, tp), denominators):
        expected = numerator / denominator if denominator else 0.0
        actual = score[metric]
        if not math.isfinite(actual) or not math.isclose(actual, expected, rel_tol=0, abs_tol=1e-12):
            raise ValueError("Metric disagrees with pair counts")
    if score.get("undefined_ratios") != [m for m, d in zip(METRICS, denominators) if not d]:
        raise ValueError("Undefined ratio flags disagree with counts")


def paired_seed_summary(comparator, baseline):
    """Resample complete paired seeds; average their metrics, not pooled pairs."""
    if set(comparator) != set(baseline):
        raise ValueError("Paired seed sets differ")
    seeds = sorted(comparator)
    if not seeds:
        return {"status": "no_complete_pairs", "included_seeds": []}
    a = np.array([[comparator[s][m] for m in METRICS] for s in seeds])
    b = np.array([[baseline[s][m] for m in METRICS] for s in seeds])
    differences = 100 * (a - b)
    result = {"status": "estimated" if len(seeds) >= 2 else "insufficient_seeds",
              "included_seeds": seeds, "metrics": {}}
    intervals = None
    if len(seeds) >= 2:
        # Reset per contrast: identical included seeds use identical draws.
        rng = np.random.Generator(np.random.PCG64(20261031))
        weights = rng.multinomial(len(seeds), np.full(len(seeds), 1 / len(seeds)), size=20000)
        draws = weights @ differences / len(seeds)
        intervals = (np.quantile(draws, [0.025, 0.975], axis=0),
                     np.quantile(draws[:, 0], [0.025 / 14, 1 - 0.025 / 14]))
    for i, metric in enumerate(METRICS):
        record = {"comparator_mean": float(a[:, i].mean()), "baseline_mean": float(b[:, i].mean()),
                  "difference_percentage_points": float(differences[:, i].mean()),
                  "seed_differences_percentage_points": differences[:, i].tolist(),
                  "inference_role": "primary" if metric == "f1" else "exploratory"}
        if intervals is not None:
            record["paired_95_percent_ci"] = intervals[0][:, i].tolist()
            if metric == "f1":
                record["bonferroni_14_ci"] = intervals[1].tolist()
        result["metrics"][metric] = record
    return result


def summarize(records):
    index = {}
    for row in records:
        key = (row["condition"], row["seed"], row["method"])
        if key in index or key[0] not in CONDITIONS or type(key[1]) is not int or key[1] not in SEEDS or key[2] not in METHODS:
            raise ValueError("Duplicate or out-of-protocol row")
        if row["status"] not in {"complete", "failed", "inapplicable"}:
            raise ValueError("Every row needs an explicit terminal outcome")
        if row["status"] == "complete":
            validate_score(row["score"])
            if not re.fullmatch(r"[0-9a-f]{64}", row.get("truth_sha256", "")):
                raise ValueError("Complete score needs truth provenance")
        elif "score" in row or not row.get("reason"):
            raise ValueError("Failures need reasons, not imputed accuracy scores")
        index[key] = row
    if len(index) != len(CONDITIONS) * len(SEEDS) * len(METHODS):
        raise ValueError("Missing protocol rows; pending work is not a failed seed")
    result = {"schema_version": 1, "publication_ready": False, "statistic": "mean of seed-level metrics",
              "bootstrap": {"replicates": 20000, "seed": 20261031, "rng": "PCG64 multinomial",
                            "reset_per_contrast": True, "f1_multiplicity_count": 14, "numpy_version": np.__version__},
              "conditions": {}, "records": records,
              "limitations": ["Ten planned seeds provide limited tail resolution.",
                              "Intervals assume exchangeable independent simulation seeds, not independent gene pairs.",
                              "Complete-case contrasts are conditional and may be biased by failures.",
                              "Bonferroni percentile intervals are approximate, not exact simultaneous coverage."]}
    for condition in CONDITIONS:
        for seed in SEEDS:
            complete = [index[(condition, seed, m)] for m in METHODS if index[(condition, seed, m)]["status"] == "complete"]
            signatures = {(r["truth_sha256"], r["score"]["input_genes"], r["score"]["eligible_true_pairs"]) for r in complete}
            if len(signatures) > 1:
                raise ValueError("Methods have different truth or input universes")
        block = {"methods": {}, "contrasts": {}}
        for method in METHODS:
            rows = [index[(condition, s, method)] for s in SEEDS]
            successes = [r for r in rows if r["status"] == "complete"]
            block["methods"][method] = {"complete_seeds": [r["seed"] for r in successes],
                "failed_seeds": [r["seed"] for r in rows if r["status"] == "failed"],
                "inapplicable_seeds": [r["seed"] for r in rows if r["status"] == "inapplicable"],
                "failure_fraction_of_planned": sum(r["status"] == "failed" for r in rows) / len(SEEDS),
                "available_case_means": {m: float(np.mean([r["score"][m] for r in successes])) if successes else None for m in METRICS}}
        for method in METHODS[:2]:
            included, excluded = [], []
            for seed in SEEDS:
                pair = [index[(condition, seed, m)] for m in (method, "orthofinder_full")]
                if all(r["status"] == "complete" for r in pair):
                    included.append(seed)
                else:
                    excluded.append({"seed": seed, "outcomes": {r["method"]: {k: r[k] for k in ("status", "reason") if k in r} for r in pair}})
            summary = paired_seed_summary({s: index[(condition, s, method)]["score"] for s in included},
                                          {s: index[(condition, s, "orthofinder_full")]["score"] for s in included})
            summary["excluded_seeds"] = excluded
            summary["conditional_on_success"] = bool(excluded)
            block["contrasts"][method] = summary
        result["conditions"][condition] = block
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--records", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("Refusing to overwrite analysis")
    data = args.records.read_bytes()
    result = summarize(json.loads(data)["records"])
    result["input"] = {"path": str(args.records.resolve()), "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
    source = Path(__file__).read_bytes()
    result["source"] = {"path": str(Path(__file__).resolve()), "sha256": hashlib.sha256(source).hexdigest()}
    result["python"] = sys.version
    result["command"] = [sys.executable, *sys.argv]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
