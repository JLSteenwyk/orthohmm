"""Reproduce complete OB bootstrap via explicit per-family accumulation."""

import argparse
from fractions import Fraction
import json
from pathlib import Path

import numpy as np

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def numerical_error(calculated, retained):
    left, right = np.asarray(calculated, dtype=float), np.asarray(retained, dtype=float)
    if left.shape != right.shape or not np.isfinite(left).all() or not np.isfinite(right).all():
        raise ValueError("Nonfinite or wrong-shaped numerical result")
    error = float(np.max(np.abs(left - right)))
    if error > 1e-10:
        raise ValueError("Numerical result differs")
    return error


def family_outcomes(method_records, baseline_records):
    def index(records):
        indexed = {row["refog"]: row for row in records}
        if not indexed or len(indexed) != len(records):
            raise ValueError("Duplicate or missing families")
        return indexed

    method, baseline = index(method_records), index(baseline_records)
    if set(method) != set(baseline):
        raise ValueError("Family panels differ")
    outcomes = dict(family_f1_wins=0, family_f1_ties=0, family_f1_losses=0)
    for name in method:
        values = []
        for row in (method[name], baseline[name]):
            tp, fp, fn = [Fraction(str(row[k])) for k in
                          ("true_positive", "false_positive", "false_negative")]
            genes = row["genes"]
            if (type(genes) is not int or genes < 2 or min(tp, fp, fn) < 0
                    or tp + fn != genes * (genes - 1) // 2):
                raise ValueError("Invalid family sufficient statistics")
            values.append(200 * tp / (2 * tp + fp + fn))
        if method[name]["genes"] != baseline[name]["genes"]:
            raise ValueError("Family sizes differ")
        difference = values[0] - values[1]
        key = ("family_f1_ties" if abs(difference) <= Fraction("1e-10") else
               "family_f1_wins" if difference > 0 else "family_f1_losses")
        outcomes[key] += 1
    return outcomes


def check_family_outcomes(calculated, retained):
    for key, value in calculated.items():
        if type(retained.get(key)) is not int or retained[key] != value:
            raise ValueError("Family outcome count differs")


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
    baseline = "orthofinder_3_1_5_full"
    metrics = {"f_score", "precision", "recall"}
    if (len(names) != 70 or len(set(names)) != 70 or len(scores) != 8 or report["replicates"] != 100000
            or report["seed"] != 20260928 or report["alpha"] != .05):
        raise ValueError("Changed frozen analysis")
    if (report["baseline"] != baseline or set(report["comparisons"]) != set(scores) - {baseline}
            or set(report["point_estimates_percent"]) != set(scores)):
        raise ValueError("Changed method inventory or baseline")
    draws = np.random.Generator(np.random.PCG64(20260928)).multinomial(70, np.full(70, 1 / 70), size=100000)
    statistics, maximum_error = {}, 0.
    for method, score in scores.items():
        indexed = {r["refog"]: r for r in score["refog_records"]}
        if len(indexed) != len(score["refog_records"]) or set(indexed) != set(names):
            raise ValueError("Changed family inventory")
        if set(report["point_estimates_percent"][method]) != metrics:
            raise ValueError("Changed point metric inventory")
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
            error = numerical_error(float(point[metric]), report["point_estimates_percent"][method][metric])
            maximum_error = max(maximum_error, error)
        statistics[method] = ratios(totals)
    rows = []
    family_checks = {}
    for method, comparison in report["comparisons"].items():
        if comparison["versus"] != baseline or set(comparison["metrics"]) != metrics:
            raise ValueError("Changed comparison inventory")
        outcomes = family_outcomes(scores[method]["refog_records"], scores[baseline]["refog_records"])
        check_family_outcomes(outcomes, comparison)
        family_checks[method] = outcomes
        for metric, retained in comparison["metrics"].items():
            point_difference = (report["point_estimates_percent"][method][metric]
                                - report["point_estimates_percent"][report["baseline"]][metric])
            maximum_error = max(maximum_error, numerical_error(point_difference, retained["difference_percentage_points"]))
            differences = statistics[method][metric] - statistics[report["baseline"]][metric]
            for key, tails in (("paired_percentile_ci", [.025, .975]),
                               ("bonferroni_percentile_ci", [.05 / 42, 1 - .05 / 42])):
                calculated = np.quantile(differences, tails)
                error = numerical_error(calculated, retained[key])
                maximum_error = max(maximum_error, error)
            rows.append(dict(method=method, metric=metric, **retained))
    if len(rows) != 21:
        raise ValueError("Incomplete endpoint inventory")
    result = dict(status="complete_ob_intervals_independently_reproduced", endpoints=21,
        maximum_absolute_error_percentage_points=maximum_error,
        sources=[record(source.resolve()), record(Path(__file__).resolve())],
        rows=rows, family_outcomes=family_checks, publication_ready=False,
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
    lines += ["", "Family F1 comparisons use rational counts and the retained 1e-10 percentage-point tie tolerance.", "",
              "| Method | Family F1 wins | Ties | Losses |", "|---|---:|---:|---:|"]
    for method, counts in family_checks.items():
        lines.append(f"| {method} | {counts['family_f1_wins']} | {counts['family_f1_ties']} | {counts['family_f1_losses']} |")
    (output / "TABLE.md").write_text("\n".join(lines) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = verify(args.base, args.output)
    print(json.dumps({k: result[k] for k in ("endpoints", "maximum_absolute_error_percentage_points")}))
