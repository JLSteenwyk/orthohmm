"""Independently check frozen seed-level simulation arithmetic without native reads."""

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import re
import sys

import numpy as np


METHODS = ("orthohmm_high_sensitivity", "orthohmm_satellite_v2", "orthofinder_full", "orthofinder_sequence_only")
CONDITIONS = ("baseline", "divergent", "turnover", "divergent_turnover", "missing20", "uneven_taxa", "taxon_count_control")
METRICS = ("f1", "precision", "recall")
PANELS = {
    "fixed_length_v1": dict(seeds=list(range(20261001, 20261011)), bootstrap_seed=20261031,
                            bytes=795597, sha256="305dcb1dde0c0f57d8b390b0f00cc6148d95103a7efb98bd9f6dce37b96e64a8"),
    "variable_length_v2": dict(seeds=list(range(20261101, 20261111)), bootstrap_seed=20261130,
                               bytes=2217675, sha256="cc99fc31c3433809098212d3c8dc12f4ad829b4cfeb66dabcb13aede4e95258f"),
}


def close(actual, expected):
    actual, expected = np.asarray(actual), np.asarray(expected)
    if (actual.shape != expected.shape or actual.dtype.kind not in "iuf"
            or not np.isfinite(actual).all()
            or not np.allclose(actual, expected, rtol=0, atol=1e-12)):
        raise ValueError("Simulation arithmetic differs from retained statistics")


def ratio(n, d):
    return Fraction(n, d) if d else Fraction(0)


def score_values(score):
    keys = ("tp", "fp", "fn", "input_genes", "eligible_true_pairs", "predicted_pairs",
            "genes_in_predicted_pairs", "duplicate_prediction_rows")
    if any(type(score[k]) is not int or score[k] < 0 for k in keys):
        raise ValueError("Require nonnegative integer simulation counts")
    tp, fp, fn, genes = (score[k] for k in ("tp", "fp", "fn", "input_genes"))
    if (tp + fn != score["eligible_true_pairs"] or tp + fp != score["predicted_pairs"]
            or max(tp + fn, tp + fp) > genes * (genes - 1) // 2
            or score["genes_in_predicted_pairs"] > genes
            or score["coverage_definition"] != "input genes occurring in at least one predicted cross-species pair"):
        raise ValueError("Pair count/universe or endpoint coverage metadata differs")
    close(score["pair_endpoint_coverage"], float(ratio(score["genes_in_predicted_pairs"], genes)))
    denominators = (2 * tp + fp + fn, tp + fp, tp + fn)
    values = [float(ratio(n, d)) for n, d in zip((2 * tp, tp, tp), denominators)]
    if (any(type(score[m]) not in (int, float) or not math.isfinite(score[m]) for m in METRICS)
            or score["undefined_ratios"] != [m for m, d in zip(METRICS, denominators) if not d]):
        raise ValueError("Metric types or undefined-ratio flags differ")
    close([score[m] for m in METRICS], values)
    return values


def checked_records(report):
    panel = report["panel"]
    if panel not in PANELS:
        raise ValueError("Unknown frozen simulation panel")
    controls = PANELS[panel]
    expected_bootstrap = dict(replicates=20000, seed=controls["bootstrap_seed"], rng="PCG64 multinomial",
                              reset_per_contrast=True, f1_multiplicity_count=14, numpy_version="2.2.6")
    if (report["schema_version"] != 1 or report["publication_ready"] is not False
            or report["statistic"] != "mean of seed-level metrics" or report["bootstrap"] != expected_bootstrap
            or set(report["conditions"]) != set(CONDITIONS)):
        raise ValueError("Frozen seed-level design or bootstrap controls differ")
    rows, scores = {}, {}
    for row in report["records"]:
        key = (row["condition"], row["seed"], row["method"])
        if (key in rows or key[0] not in CONDITIONS or type(key[1]) is not int
                or key[1] not in controls["seeds"] or key[2] not in METHODS):
            raise ValueError("Duplicate, missing or out-of-protocol simulation identity")
        if row["status"] == "complete":
            if not re.fullmatch("[0-9a-f]{64}", row.get("truth_sha256", "")):
                raise ValueError("Missing scored truth identity")
            scores[key] = score_values(row["score"])
        elif (row["status"] not in ("failed", "inapplicable") or "score" in row
              or not isinstance(row.get("reason"), str) or not row["reason"]):
            raise ValueError("Explicit terminal failures need reasons, not imputed scores")
        rows[key] = row
    if len(rows) != 280:
        raise ValueError("Incomplete planned simulation inventory")
    for condition in CONDITIONS:
        for seed in controls["seeds"]:
            admitted = [rows[condition, seed, m] for m in METHODS if (condition, seed, m) in scores]
            universes = {(r["truth_sha256"], r["score"]["input_genes"], r["score"]["eligible_true_pairs"]) for r in admitted}
            if len(universes) > 1:
                raise ValueError("Methods do not share truth/input universes")
            if (condition, seed, METHODS[3]) in scores and (condition, seed, METHODS[2]) not in scores:
                raise ValueError("Diagnostic checkpoint has no admitted full parent")
    return controls, rows, scores


def paired_intervals(differences, bootstrap_seed):
    differences = np.asarray(differences, dtype=float)
    n = len(differences)
    weights = np.random.Generator(np.random.PCG64(bootstrap_seed)).multinomial(n, np.full(n, 1 / n), size=20000)
    # Sum each seed's contribution explicitly, independently of the producer's matrix product.
    draws = np.zeros((20000, 3))
    for i in range(n):
        draws += weights[:, i, None] * differences[i] / n
    nominal = np.quantile(draws, [.025, .975], axis=0, method="linear")
    adjusted = np.quantile(draws[:, 0], [.025 / 14, 1 - .025 / 14], method="linear")
    return nominal, adjusted


def verify(report):
    controls, records, scores = checked_records(report)
    if np.__version__ != report["bootstrap"]["numpy_version"]:
        raise ValueError("Use the retained NumPy 2.2.6 numerical implementation")
    estimated = absent = insufficient = method_cells = bound_values = paired_effects = 0
    for condition in CONDITIONS:
        block = report["conditions"][condition]
        if set(block["methods"]) != set(METHODS) or set(block["contrasts"]) != set(METHODS[:2]):
            raise ValueError("Wrong condition method/contrast inventory")
        for method in METHODS:
            row = block["methods"][method]
            expected_fields = {"complete_seeds", "failed_seeds", "inapplicable_seeds", "failure_fraction_of_planned", "available_case_means"}
            if set(row) != expected_fields or set(row["available_case_means"]) != set(METRICS):
                raise ValueError("Wrong method summary fields")
            for status in ("complete", "failed", "inapplicable"):
                seeds = [s for s in controls["seeds"] if records[condition, s, method]["status"] == status]
                if row[status + "_seeds"] != seeds:
                    raise ValueError("Terminal outcome inventory differs")
            close(row["failure_fraction_of_planned"], len(row["failed_seeds"]) / 10)
            values = [scores[condition, s, method] for s in row["complete_seeds"]]
            for i, metric in enumerate(METRICS):
                expected = math.fsum(v[i] for v in values) / len(values) if values else None
                actual = row["available_case_means"][metric]
                if expected is None:
                    if actual is not None:
                        raise ValueError("Unavailable method mean was imputed")
                else:
                    close(actual, expected)
                method_cells += 1
        for method in METHODS[:2]:
            summary = block["contrasts"][method]
            paired = [s for s in controls["seeds"] if all((condition, s, m) in scores for m in (method, METHODS[2]))]
            excluded = [dict(seed=s, outcomes={m: {k: records[condition, s, m][k] for k in ("status", "reason") if k in records[condition, s, m]}
                                             for m in (method, METHODS[2])}) for s in controls["seeds"] if s not in paired]
            if (summary["included_seeds"] != paired or summary["excluded_seeds"] != excluded
                    or summary["conditional_on_success"] is not bool(excluded)):
                raise ValueError("Paired seed inclusion, exclusions or conditioning differ")
            status = "no_complete_pairs" if not paired else "insufficient_seeds" if len(paired) == 1 else "estimated"
            expected_fields = {"status", "included_seeds", "excluded_seeds", "conditional_on_success"} | ({"metrics"} if paired else set())
            if summary["status"] != status or set(summary) != expected_fields:
                raise ValueError("Wrong complete-case inferential status")
            if not paired:
                absent += 1
                continue
            a = np.array([scores[condition, s, method] for s in paired])
            b = np.array([scores[condition, s, METHODS[2]] for s in paired])
            differences = 100 * (a - b)
            intervals = paired_intervals(differences, controls["bootstrap_seed"]) if len(paired) >= 2 else None
            estimated += intervals is not None
            insufficient += intervals is None
            if set(summary["metrics"]) != set(METRICS):
                raise ValueError("Wrong paired metric inventory")
            for i, metric in enumerate(METRICS):
                item = summary["metrics"][metric]
                fields = {"comparator_mean", "baseline_mean", "difference_percentage_points", "seed_differences_percentage_points", "inference_role"}
                if intervals is not None:
                    fields.add("paired_95_percent_ci")
                    if i == 0:
                        fields.add("bonferroni_14_ci")
                if set(item) != fields or item["inference_role"] != ("primary" if i == 0 else "exploratory"):
                    raise ValueError("Wrong metric role, adjustment or interval availability")
                close(item["comparator_mean"], math.fsum(a[:, i]) / len(paired))
                close(item["baseline_mean"], math.fsum(b[:, i]) / len(paired))
                close(item["difference_percentage_points"], math.fsum(differences[:, i]) / len(paired))
                close(item["seed_differences_percentage_points"], differences[:, i])
                paired_effects += 1
                if intervals is not None:
                    close(item["paired_95_percent_ci"], intervals[0][:, i])
                    bound_values += 2
                    if i == 0:
                        close(item["bonferroni_14_ci"], intervals[1])
                        bound_values += 2
    return dict(panel=report["panel"], outcome_records=280, scored_records=len(scores),
                failed_records=sum(r["status"] == "failed" for r in records.values()),
                inapplicable_records=sum(r["status"] == "inapplicable" for r in records.values()),
                method_mean_cells=method_cells, planned_contrasts=14, estimated_contrasts=estimated,
                unavailable_contrasts=absent, insufficient_seed_contrasts=insufficient,
                paired_metric_effects=paired_effects, interval_bound_values=bound_values)


def reproduce(results, output):
    output = Path(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    result = dict(status="validation_failed", publication_ready=False, native_inference_repeated=False,
                  raw_scoring_repeated=False, historical_evidence_paths_accessed=False)
    try:
        reports, contents = [], []
        if len(results) != 2:
            raise ValueError("Supply the two distinct frozen panels")
        for path in map(Path, results):
            if path.is_symlink() or not path.is_file() or path.stat().st_size > 3 * 1024 ** 2:
                raise ValueError("Require bounded direct report files")
            raw = path.read_bytes()
            report = json.loads(raw)
            expected = PANELS.get(report["panel"])
            if expected is None or (len(raw), hashlib.sha256(raw).hexdigest()) != (expected["bytes"], expected["sha256"]):
                raise ValueError("Frozen simulation report identity differs")
            reports.append(report)
            contents.append(raw)
        if {r["panel"] for r in reports} != set(PANELS):
            raise ValueError("Duplicate or mixed panel inputs")
        summaries = [verify(report) for report in reports]
        if any(Path(path).read_bytes() != raw for path, raw in zip(results, contents)):
            raise ValueError("Input changed during numerical reproduction")
        source = Path(__file__).read_bytes()
        result.update(status="frozen_simulation_panels_arithmetically_reproduced", panels=summaries,
                      inputs=[dict(bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest()) for raw in contents],
                      source=dict(bytes=len(source), sha256=hashlib.sha256(source).hexdigest()),
                      numpy_version=np.__version__, python=sys.version, absolute_tolerance=1e-12,
                      panels_pooled=False, failed_scores_imputed=False, replicates_per_estimated_contrast=20000,
                      limitations=["Retained-count arithmetic only, not native output, truth history or generation re-admission.",
                          "Independent rational per-seed metrics and explicit seed sums; same NumPy PCG64/linear quantiles, not a new statistical engine.",
                          "Complete-case contrasts remain conditional on success; ten planned seeds limit interval resolution.",
                          "Fixed-length unavailable contrasts remain unavailable, not zero or tool wins.",
                          "No controlled resources, arbitrary-dataset validation, complete release or rights clearance."])
    except Exception as error:
        result.update(error_type=type(error).__name__, error=str(error))
        raise
    finally:
        with output.open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, action="append", required=True, help="Supply each frozen panel once")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(reproduce(args.results, args.output), indent=2, sort_keys=True))
