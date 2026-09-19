"""Numerical engine for the prespecified corrected-QfO parameter contrasts.

This module does not admit source counts or write publication evidence. A
separate provenance-bound count audit is required before applying it to runs.
"""

import numpy as np

from benchmark_tools.audit_qfo_swiss_counts import LABELS, REFERENCE_SHA, statistics
from benchmark_tools.bootstrap_qfo_swiss_stages import METRICS, aggregate

ARMS = ("control", "cpm_low", "cpm_high", "norm_low", "norm_high", "margin_low", "margin_high")
REPLICATES = 100000
SEED = 20260925
MULTIPLICITY = 18


def validated_values(report):
    if (report["status"] != "corrected_qfo_parameter_swiss_counts_verified"
            or report["publication_ready"] is not False or report["uncertainty_admitted"] is not False
            or report["shared_represented_genes"] or report["reference"]["sha256"] != REFERENCE_SHA):
        raise ValueError("Require audited disjoint corrected SwissTrees counts")
    families = report["families"]
    if (len(families) != 18 or len(set(families)) != 18
            or [r["arm"] for r in report["arms"]] != list(ARMS)):
        raise ValueError("Require 18 families and all seven planned arms in order")
    if type(report["reference_relation_count"]) is not int or report["reference_relation_count"] <= 0:
        raise ValueError("Invalid reference relation count")
    values, universe = {}, None
    for arm in report["arms"]:
        if arm["status"] == "not_admitted":
            if (not isinstance(arm.get("reason"), str) or not arm["reason"].strip()
                    or "families" in arm or "aggregate" in arm):
                raise ValueError("Unavailable arm must state reason without imputed scores")
            continue
        if arm["status"] != "counts_verified" or [r["family"] for r in arm["families"]] != families:
            raise ValueError("Wrong arm status or family order")
        seen, current, reference, total = set(), [], [], 0
        for row in arm["families"]:
            counts, genes = row["counts_without_prior"], row["represented_genes"]
            if set(counts) != set(LABELS) or any(type(v) is not int or v < 0 for v in counts.values()):
                raise ValueError("Invalid raw confusion counts")
            if not sum(counts.values()) or len(genes) <= 5 or len(set(genes)) != len(genes) or seen.intersection(genes):
                raise ValueError("Empty, duplicate or overlapping family")
            seen.update(genes)
            reference.append((tuple(sorted(genes)), counts["TP"] + counts["FN"], counts["FP"] + counts["TN"]))
            calculated = statistics(counts)
            stored = row["statistics_with_prior"]
            if set(stored) != set(METRICS) or any(
                    type(stored[m]) not in (int, float) or not np.isfinite(stored[m])
                    or abs(calculated[m] - stored[m]) > 1e-12 for m in METRICS):
                raise ValueError("Stored family statistic differs from raw counts")
            current.append([calculated["PPV"], calculated["TPR"]])
            total += sum(counts.values())
        if universe is None:
            universe = reference
        if reference != universe or total != report["reference_relation_count"]:
            raise ValueError("Reference genes or truth totals differ across arms")
        array = np.asarray(current)
        point = aggregate(array.mean(axis=0))
        stored = arm["aggregate"]
        if set(stored) != set(METRICS) or any(
                type(stored[m]) not in (int, float) or not np.isfinite(stored[m])
                or abs(point[j] - stored[m]) > 1e-12 for j, m in enumerate(METRICS)):
            raise ValueError("Aggregate is not the native macro statistic")
        values[arm["arm"]] = array
    return values


def calculate(report, *, replicates=REPLICATES, seed=SEED):
    if type(replicates) is not int or replicates < 100 or type(seed) is not int or seed < 0:
        raise ValueError("Invalid bootstrap controls")
    values = validated_values(report)
    n = len(report["families"])
    weights = np.random.Generator(np.random.PCG64(seed)).multinomial(n, np.full(n, 1 / n), size=replicates)
    draws = {name: aggregate(weights @ cell / n) for name, cell in values.items()}
    points = {name: aggregate(cell.mean(axis=0)) for name, cell in values.items()}
    comparisons = []
    for name in ARMS[1:]:
        row = {"candidate": name, "reference": "control"}
        if name not in values or "control" not in values:
            row.update(status="not_estimable", reason="baseline_not_admitted" if "control" not in values else "variant_not_admitted",
                       metrics=None, family_differences=None)
        else:
            differences = draws[name] - draws["control"]
            family_differences = aggregate(values[name]) - aggregate(values["control"])
            metrics = {}
            for j, metric in enumerate(METRICS):
                metrics[metric] = {"difference": float(points[name][j] - points["control"][j]),
                    "paired_percentile_ci": np.quantile(differences[:, j], [.025, .975], method="linear").tolist(),
                    "bonferroni_percentile_ci": np.quantile(differences[:, j], [.05 / (2 * MULTIPLICITY), 1 - .05 / (2 * MULTIPLICITY)], method="linear").tolist(),
                    "family_wins": int(np.sum(family_differences[:, j] > 1e-10)),
                    "family_ties": int(np.sum(np.abs(family_differences[:, j]) <= 1e-10)),
                    "family_losses": int(np.sum(family_differences[:, j] < -1e-10))}
            row.update(status="estimated", metrics=metrics, family_differences=[
                dict(family=family, **dict(zip(METRICS, delta.tolist())))
                for family, delta in zip(report["families"], family_differences)])
        comparisons.append(row)
    return {"status": "parameter_swiss_intervals_computed_pending_provenance_admission",
        "publication_ready": False, "uncertainty_admitted": False,
        "protocol_controls_match": replicates == REPLICATES and seed == SEED,
        "replicates": replicates, "seed": seed, "numpy_version": np.__version__,
        "rng": "numpy.PCG64 multinomial; shared family draws across all available arms",
        "alpha": .05, "multiplicity_endpoints": MULTIPLICITY, "quantile_method": "linear",
        "families": report["families"], "units": "raw 0-to-1 metric units",
        "arms": [{k: row[k] for k in ("arm", "status", "reason") if k in row} for row in report["arms"]],
        "point_estimates": {name: dict(zip(METRICS, points[name].tolist())) if name in points else None for name in ARMS},
        "comparisons": comparisons, "limitations": [
            "Development-exposed conditional sensitivity analysis, not independent or selection-adjusted confirmation.",
            "Only 18 families; shared history and merged predictions can violate exchangeability.",
            "Adjustment retains all 18 planned endpoints regardless of missing variants.",
            "Missing arms are not zeros; unavailable control prevents every paired contrast.",
            "No TreeFam family intervals or uncertainty on other endpoints or the secondary mean.",
            "Intervals containing zero do not establish equivalence or general robustness.",
            "Numerical calculation only; source-count and protocol provenance admission is required."]}
