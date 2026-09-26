"""Standalone count-level reproduction; only this script, result JSON and NumPy are needed."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np


RESULT_SHA = "d31dbc0f4e3f02dfb534c8b59b559db00260227f23f2640b4fed13d5aee8ec89"
CONDITIONS = ("baseline", "divergent", "divergent_turnover", "missing20", "taxon_count_control", "turnover", "uneven_taxa")
SEEDS = tuple(range(20261106, 20261111))
METRICS = ("f1", "precision", "recall")


def record(path):
    path = Path(path).resolve()
    data = path.read_bytes()
    return dict(path=str(path), bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def reproduce(result):
    values = {}
    for row in result["records"]:
        key = row["condition"], row["seed"], row["arm"]
        if key in values:
            raise ValueError("Duplicate result row")
        score = row["score"]
        tp, fp, fn = [score[k] for k in ("tp", "fp", "fn")]
        if any(type(n) is not int or n < 0 for n in (tp, fp, fn)) or tp + fn != score["eligible_true_pairs"]:
            raise ValueError("Invalid sufficient counts")
        denominators = [2*tp+fp+fn, tp+fp, tp+fn]
        point = np.array([n/d if d else 0 for n, d in zip((2*tp, tp, tp), denominators)])
        if not np.allclose(point, [score[m] for m in METRICS], atol=1e-12, rtol=0):
            raise ValueError("Scores disagree with counts")
        if score["undefined_ratios"] != [m for m, d in zip(METRICS, denominators) if not d]:
            raise ValueError("Undefined ratio flags differ")
        values[key] = point
    expected = {(c, s, a) for c in CONDITIONS for s in SEEDS for a in ("hmm", "diamond")}
    if set(values) != expected:
        raise ValueError("Incomplete reporting panel")
    rng = np.random.Generator(np.random.PCG64(20260927))
    weights = rng.multinomial(5, np.full(5, .2), size=20000)
    condition_arrays = {}
    for c in CONDITIONS:
        condition_arrays[c] = {arm: np.stack([values[c, seed, arm] for seed in SEEDS]) for arm in ("hmm", "diamond")}
    condition_arrays["overall"] = {arm: sum(condition_arrays[c][arm] for c in CONDITIONS)/7 for arm in ("hmm", "diamond")}
    checked = 0
    for condition, arms in condition_arrays.items():
        differences = 100 * (arms["hmm"] - arms["diamond"])
        boot = weights @ differences / 5
        for i, metric in enumerate(METRICS):
            row = result["contrasts"][condition][metric]
            checks = [(row["hmm_mean"], arms["hmm"][:, i].mean()),
                      (row["diamond_mean"], arms["diamond"][:, i].mean()),
                      (row["difference_percentage_points"], differences[:, i].mean()),
                      (row["seed_differences_percentage_points"], differences[:, i]),
                      (row["marginal_95_percent_ci"], np.quantile(boot[:, i], [.025, .975]))]
            if metric == "f1":
                checks.append((row["bonferroni_8_ci"], np.quantile(boot[:, i], [.025/8, 1-.025/8])))
            if any(not np.allclose(a, b, atol=1e-12, rtol=0) for a, b in checks):
                raise ValueError("Paired effect/interval reproduction differs")
            signs = [int((differences[:, i] > 0).sum()), int((differences[:, i] == 0).sum()), int((differences[:, i] < 0).sum())]
            if signs != [row[k] for k in ("wins", "ties", "losses")]:
                raise ValueError("Seed win/tie/loss counts differ")
            checked += 1
    return dict(status="all_24_metric_contrasts_reproduced_from_counts", records=len(values),
                metric_contrasts=checked, f1_adjusted_intervals=8, absolute_tolerance=1e-12,
                numpy_version=np.__version__, python=sys.version,
                source=record(__file__), original_data_paths_opened=False,
                limitations=["Only count-level statistics; does not revalidate predictions, simulation truth or inference.",
                             "Five seed blocks and development exposure remain limitations.",
                             "NumPy and Python are runtime prerequisites, not bundled dependencies."],
                publication_ready=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    identity = record(args.results)
    if identity["sha256"] != RESULT_SHA:
        raise ValueError("Unexpected retained result identity")
    result = reproduce(json.loads(args.results.read_text()))
    result["input"] = identity
    if record(args.results) != identity:
        raise ValueError("Result changed during reproduction")
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
