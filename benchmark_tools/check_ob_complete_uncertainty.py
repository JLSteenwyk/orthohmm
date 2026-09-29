"""Reproduce complete OB bootstrap via explicit per-family accumulation."""

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

def file_record(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return dict(path=str(path.resolve()), bytes=path.stat().st_size, sha256=digest.hexdigest())


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
    from benchmark_tools.prepare_ob_candidate_neighborhood import check

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
    return verify_data(report, scores, output, [file_record(source), file_record(Path(__file__))],
                       "local_raw_input_records_checked")


def export_portable(base, output):
    from benchmark_tools.bootstrap_ob_complete import load_scores
    from benchmark_tools.prepare_ob_candidate_neighborhood import check

    if output.exists():
        raise FileExistsError(output)
    scores, refs = load_scores(base)
    source = base / "ob_complete_uncertainty_20260928.json"
    report = json.loads(source.read_text())
    for ref in report["checked_records"]:
        check(ref)
    columns = ("refog", "genes", "true_positive", "false_positive", "false_negative")
    compact = {method: {"refog_records": [{key: row[key] for key in columns}
               for row in score["refog_records"]]} for method, score in scores.items()}
    report_keys = ("baseline", "families", "replicates", "seed", "alpha", "rng",
                   "numpy_version", "point_estimates_percent", "comparisons", "multiplicity", "limitations")
    payload = dict(schema="orthohmm_ob_complete_statistics_v1", scores=compact,
        report={key: report[key] for key in report_keys},
        provenance=[dict(name=Path(ref["path"]).name, bytes=ref["bytes"], sha256=ref["sha256"])
                    for ref in [file_record(source), *refs[:2]]],
        limitations="Derived family sufficient statistics only. No sequences, gene memberships or raw-input verification on replay; not an inference package or independent biological confirmation.")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    return file_record(output)


def verify_portable(source, expected_sha256, output):
    if output.exists():
        raise FileExistsError(output)
    data = source.read_bytes()
    if hashlib.sha256(data).hexdigest() != expected_sha256:
        raise ValueError("Portable input checksum differs")
    payload = json.loads(data)
    if payload["schema"] != "orthohmm_ob_complete_statistics_v1":
        raise ValueError("Unknown portable schema")
    return verify_data(payload["report"], payload["scores"], output,
        [dict(path=str(source.resolve()), bytes=len(data), sha256=expected_sha256), file_record(Path(__file__))],
        "portable_derived_counts_only_no_raw_input_verification")


def verify_data(report, scores, output, sources, verification_scope):
    if output.exists():
        raise FileExistsError(output)
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
        sources=sources, verification_scope=verification_scope,
        runtime=dict(python=sys.version, numpy=np.__version__),
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
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--base", type=Path)
    inputs.add_argument("--export-from", type=Path)
    inputs.add_argument("--portable", type=Path)
    parser.add_argument("--sha256")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if bool(args.portable) != bool(args.sha256):
        parser.error("--sha256 is required only with --portable")
    if args.export_from:
        result = export_portable(args.export_from, args.output)
        print(json.dumps(result))
    else:
        result = (verify_portable(args.portable, args.sha256, args.output) if args.portable
                  else verify(args.base, args.output))
        print(json.dumps({k: result[k] for k in ("endpoints", "maximum_absolute_error_percentage_points")}))
