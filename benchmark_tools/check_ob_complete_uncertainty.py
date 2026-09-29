"""Reproduce complete OB bootstrap via explicit per-family accumulation."""

import argparse
from fractions import Fraction
import json
from pathlib import Path

import numpy as np

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def ratios(counts):
    tp, fp, fn = counts
    return dict(f_score=100 * 2 * tp / (2 * tp + fp + fn),
                precision=100 * tp / (tp + fp), recall=100 * tp / (tp + fn))


def verify(base, output):
    if output.exists():
        raise FileExistsError(output)
    source = base / "ob_complete_uncertainty_20260928.json"
    report = json.loads(source.read_text())
    for ref in report["checked_records"]:
        check(ref)
    main = json.loads((base / "orthobench_paired_uncertainty_20260916.json").read_text())
    other = json.loads((base / "retained_ob_comparator_readback_20260926.json").read_text())
    scores = dict(main["scores"])
    scores.update({r["key"]: r["score"] for r in other["rows"]})
    names = sorted(report["families"])
    if (len(names) != 70 or len(scores) != 8 or report["replicates"] != 100000
            or report["seed"] != 20260928 or report["alpha"] != .05):
        raise ValueError("Changed frozen analysis")
    draws = np.random.Generator(np.random.PCG64(20260928)).multinomial(70, np.full(70, 1 / 70), size=100000)
    statistics, maximum_error = {}, 0.
    for method, score in scores.items():
        indexed = {r["refog"]: r for r in score["refog_records"]}
        totals = np.zeros((3, 100000))
        exact = [Fraction(), Fraction(), Fraction()]
        for i, name in enumerate(names):
            row = indexed[name]
            for j, key in enumerate(("true_positive", "false_positive", "false_negative")):
                weight = Fraction(str(row[key])) / (row["genes"] - 1)
                totals[j] += draws[:, i] * float(weight)
                exact[j] += weight
        point = ratios(exact)
        for metric in point:
            error = abs(float(point[metric]) - report["point_estimates_percent"][method][metric])
            maximum_error = max(maximum_error, error)
            if error > 1e-10:
                raise ValueError("Point estimate differs")
        statistics[method] = ratios(totals)
    rows = []
    for method, comparison in report["comparisons"].items():
        for metric, retained in comparison["metrics"].items():
            point_difference = (report["point_estimates_percent"][method][metric]
                                - report["point_estimates_percent"][report["baseline"]][metric])
            if abs(point_difference - retained["difference_percentage_points"]) > 1e-10:
                raise ValueError("Point difference differs")
            differences = statistics[method][metric] - statistics[report["baseline"]][metric]
            for key, tails in (("paired_percentile_ci", [.025, .975]),
                               ("bonferroni_percentile_ci", [.05 / 42, 1 - .05 / 42])):
                calculated = np.quantile(differences, tails)
                error = float(np.max(np.abs(calculated - retained[key])))
                maximum_error = max(maximum_error, error)
                if error > 1e-10:
                    raise ValueError("Interval differs")
            rows.append(dict(method=method, metric=metric, **retained))
    if len(rows) != 21:
        raise ValueError("Incomplete endpoint inventory")
    result = dict(status="complete_ob_intervals_independently_reproduced", endpoints=21,
        maximum_absolute_error_percentage_points=maximum_error,
        sources=[record(source.resolve()), record(Path(__file__).resolve())],
        rows=rows, publication_ready=False,
        limitations="Alternative count accumulation and rational points; same NumPy RNG/quantile implementation. Not statistical coverage validation.")
    output.mkdir(parents=True)
    (output / "crosscheck.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    lines = ["# Complete OrthoBench Paired Comparisons", "",
             "Exploratory, development-exposed; conditional on exchangeable RefOGs. Differences in percentage points versus full OrthoFinder 3.1.5.", "",
             "| Method | Metric | Difference | Nominal 95% CI | Adjusted CI (21 endpoints) |",
             "|---|---|---:|---|---|"]
    for row in rows:
        intervals = [f"[{row[k][0]:.3f}, {row[k][1]:.3f}]" for k in ("paired_percentile_ci", "bonferroni_percentile_ci")]
        lines.append(f"| {row['method']} | {row['metric']} | {row['difference_percentage_points']:.3f} | {intervals[0]} | {intervals[1]} |")
    lines += ["", "100,000 paired draws; PCG64 seed 20260928. Approximate percentile intervals, not guarantees of coverage. No tuning or new independent confirmation."]
    (output / "TABLE.md").write_text("\n".join(lines) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = verify(args.base, args.output)
    print(json.dumps({k: result[k] for k in ("endpoints", "maximum_absolute_error_percentage_points")}))
