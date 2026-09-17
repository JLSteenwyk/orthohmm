"""Prespecified 126-endpoint tree robustness analysis with paired seed resampling."""

import argparse
import json
from pathlib import Path
import re
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.summarize_simulation_panel import CONDITIONS, METRICS, validate_score

METHODS = ("orthohmm_satellite_v2", "orthofinder_full")
SEEDS = tuple(range(20261101, 20261111))
ARMS = ("inferred", "generating", "nni1", "nni2")
CONTRASTS = (("generating", "inferred"), ("nni1", "generating"), ("nni2", "generating"))
BOOTSTRAP_SEED = 20260918
REPLICATES = 20000
MULTIPLICITY = 126


def index_records(records):
    expected = {(c, s, m, a) for c in CONDITIONS for s in SEEDS for m in METHODS for a in ARMS}
    indexed = {}
    for row in records:
        key = row["condition"], row["seed"], row["method"], row["arm"]
        if key not in expected or key in indexed or type(row["seed"]) is not int:
            raise ValueError("Duplicate or out-of-protocol tree score")
        if row["status"] == "complete":
            validate_score(row["score"])
            if not re.fullmatch(r"[0-9a-f]{64}", row.get("truth_sha256", "")):
                raise ValueError("Completed score lacks truth identity")
        elif row["status"] != "failed" or "score" in row or not row.get("reason"):
            raise ValueError("Every planned arm needs a terminal score or explicit failure without imputation")
        indexed[key] = row
    if set(indexed) != expected:
        raise ValueError("Incomplete 560-row panel; do not analyze partial results")
    for condition in CONDITIONS:
        for seed in SEEDS:
            available = [indexed[condition, seed, m, a] for m in METHODS for a in ARMS
                         if indexed[condition, seed, m, a]["status"] == "complete"]
            if len({(r["truth_sha256"], r["score"]["input_genes"], r["score"]["eligible_true_pairs"])
                    for r in available}) > 1:
                raise ValueError("Tree arms use different truth or input universes")
    return indexed


def paired_summary(target, reference):
    if set(target) != set(reference):
        raise ValueError("Paired seed sets differ")
    seeds = sorted(target)
    result = {"included_seeds": seeds, "paired_seed_count": len(seeds), "metrics": {}}
    if not seeds:
        result["status"] = "no_complete_pairs"
        result["metrics"] = {m: {"difference_percentage_points": None, "paired_95_percent_ci": None,
                                  "bonferroni_126_ci": None} for m in METRICS}
        return result
    a = np.array([[target[s][m] for m in METRICS] for s in seeds])
    b = np.array([[reference[s][m] for m in METRICS] for s in seeds])
    differences = 100 * (a - b)
    intervals = None
    if len(seeds) >= 2:
        rng = np.random.Generator(np.random.PCG64(BOOTSTRAP_SEED))
        weights = rng.multinomial(len(seeds), np.full(len(seeds), 1 / len(seeds)), size=REPLICATES)
        draws = weights @ differences / len(seeds)
        intervals = np.quantile(draws, [.025, .975, .025 / MULTIPLICITY, 1 - .025 / MULTIPLICITY], axis=0)
    result["status"] = "estimated" if intervals is not None else "insufficient_seeds"
    for i, metric in enumerate(METRICS):
        result["metrics"][metric] = {"target_mean": float(a[:, i].mean()), "reference_mean": float(b[:, i].mean()),
            "difference_percentage_points": float(differences[:, i].mean()),
            "seed_differences_percentage_points": differences[:, i].tolist(),
            "paired_95_percent_ci": intervals[:2, i].tolist() if intervals is not None else None,
            "bonferroni_126_ci": intervals[2:, i].tolist() if intervals is not None else None,
            "target_undefined_ratio_seeds": [s for s in seeds if metric in target[s]["undefined_ratios"]],
            "reference_undefined_ratio_seeds": [s for s in seeds if metric in reference[s]["undefined_ratios"]]}
    return result


def summarize(records):
    indexed = index_records(records)
    contrasts = []
    for condition in CONDITIONS:
        for method in METHODS:
            for target, reference in CONTRASTS:
                included, excluded = [], []
                for seed in SEEDS:
                    arms = {a: indexed[condition, seed, method, a] for a in (target, reference)}
                    if all(r["status"] == "complete" for r in arms.values()):
                        included.append(seed)
                    else:
                        excluded.append({"seed": seed, "outcomes": {a: {k: r[k] for k in ("status", "reason") if k in r}
                                                                     for a, r in arms.items()}})
                result = paired_summary({s: indexed[condition, s, method, target]["score"] for s in included},
                                        {s: indexed[condition, s, method, reference]["score"] for s in included})
                contrasts.append({"condition": condition, "method": method, "target": target, "reference": reference,
                                  "excluded_seeds": excluded, "conditional_on_success": bool(excluded), **result})
    if len(contrasts) != 42 or sum(len(c["metrics"]) for c in contrasts) != MULTIPLICITY:
        raise ValueError("Missing prespecified exploratory endpoints")
    return {"status": "tree_robustness_summarized", "publication_ready": False, "records": records, "contrasts": contrasts,
            "bootstrap": {"replicates": REPLICATES, "seed": BOOTSTRAP_SEED, "rng": "PCG64 multinomial",
                          "reset_per_contrast": True, "numpy_version": np.__version__, "multiplicity": MULTIPLICITY},
            "statistic": "Arithmetic mean of paired per-seed metric differences; percentage-point effects",
            "inference_role": "All 126 endpoints exploratory; fixed multiplicity even with unavailable contrasts",
            "limitations": ["Ten planned seeds provide limited bootstrap-tail resolution.",
                            "Bonferroni percentile intervals are approximate, not exact simultaneous coverage.",
                            "Complete-case contrasts may be biased by failures; no failed seed is assigned zero accuracy.",
                            "Undefined ratios follow frozen scorer zero conventions and remain flagged explicitly.",
                            "Development-exposed simplified simulations, not independent biological generalization."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    data = json.loads(args.scores.read_text())
    if data["status"] != "complete_tree_panel_scored":
        raise ValueError("Require independently admitted complete-panel scores")
    result = summarize(data["records"])
    result.update(source=record(__file__), scores=record(args.scores))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
