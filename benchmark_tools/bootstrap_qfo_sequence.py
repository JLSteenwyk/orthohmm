"""Numerical engine for the prespecified corrected SwissTrees search control.

This module does not admit source files or expose a production CLI. An
independent source-bound count audit is required before using real results.
"""

import numpy as np

from benchmark_tools.audit_qfo_swiss_counts import LABELS, REFERENCE_SHA, statistics
from benchmark_tools.bootstrap_qfo_swiss_stages import aggregate, METRICS

VARIANTS = ("p0_c0_r0", "all_hits", "top100")


def validated_values(report):
    if (report.get("status") != "corrected_qfo_sequence_swiss_counts_verified"
            or report.get("publication_ready") is not False
            or report.get("uncertainty_admitted") is not False
            or report["shared_represented_genes"]):
        raise ValueError("Require verified unresampled disjoint corrected counts")
    if (report["reference"]["sha256"] != REFERENCE_SHA
            or type(report["reference_relation_count"]) is not int
            or report["reference_relation_count"] != 10765):
        raise ValueError("Changed reference or relation universe")
    families = report["families"]
    if (len(families) != 18 or len(set(families)) != 18
            or [v["variant"] for v in report["variants"]] != list(VARIANTS)):
        raise ValueError("Require fixed family and variant inventory")
    values, reference = [], None
    for variant in report["variants"]:
        if [r["family"] for r in variant["families"]] != families:
            raise ValueError("Changed family order")
        seen, universe, current, total = set(), [], [], 0
        for row in variant["families"]:
            counts, genes = row["counts_without_prior"], row["represented_genes"]
            if set(counts) != set(LABELS) or any(type(v) is not int or v < 0 for v in counts.values()):
                raise ValueError("Invalid confusion counts")
            if not sum(counts.values()) or len(genes) <= 5 or len(set(genes)) != len(genes) or seen.intersection(genes):
                raise ValueError("Empty, duplicate or overlapping family")
            seen.update(genes)
            universe.append((tuple(sorted(genes)), counts["TP"] + counts["FN"], counts["FP"] + counts["TN"]))
            calculated = statistics(counts)
            stored = row["statistics_with_prior"]
            if set(stored) != set(METRICS) or any(not np.isfinite(stored[m]) or abs(stored[m] - calculated[m]) > 1e-12 for m in METRICS):
                raise ValueError("Family statistics disagree with counts")
            current.append([calculated["PPV"], calculated["TPR"]])
            total += sum(counts.values())
        if reference is None:
            reference = universe
        if universe != reference or total != 10765:
            raise ValueError("Reference membership or truth totals differ")
        point = aggregate(np.mean(current, axis=0))
        if set(variant["aggregate"]) != set(METRICS) or any(
                not np.isfinite(variant["aggregate"][m]) or abs(variant["aggregate"][m] - point[j]) > 1e-12
                for j, m in enumerate(METRICS)):
            raise ValueError("Aggregate differs from native macro statistic")
        values.append(current)
    return np.asarray(values)


def bootstrap(report, replicates=100000, seed=20260923):
    if type(replicates) is not int or replicates < 100 or type(seed) is not int or seed < 0:
        raise ValueError("Invalid bootstrap controls")
    values = validated_values(report)
    weights = np.random.Generator(np.random.PCG64(seed)).multinomial(18, np.full(18, 1 / 18), size=replicates)
    draws = np.asarray([aggregate(weights @ variant / 18) for variant in values])
    observed, per_family = aggregate(values.mean(axis=1)), aggregate(values)
    comparisons = []
    for index in (1, 2):
        delta, family_delta = draws[index] - draws[0], per_family[index] - per_family[0]
        metrics = {}
        for j, metric in enumerate(METRICS):
            metrics[metric] = {
                "difference": float(observed[index, j] - observed[0, j]),
                "paired_percentile_ci": np.quantile(delta[:, j], [.025, .975], method="linear").tolist(),
                "bonferroni_percentile_ci": np.quantile(delta[:, j], [.05 / 12, 1 - .05 / 12], method="linear").tolist(),
                "family_wins": int(np.sum(family_delta[:, j] > 1e-10)),
                "family_ties": int(np.sum(np.abs(family_delta[:, j]) <= 1e-10)),
                "family_losses": int(np.sum(family_delta[:, j] < -1e-10)),
            }
        comparisons.append({"candidate": VARIANTS[index], "reference": VARIANTS[0], "metrics": metrics,
            "family_differences": [dict(family=name, **dict(zip(METRICS, row.tolist())))
                                   for name, row in zip(report["families"], family_delta)]})
    return {"status": "corrected_sequence_swiss_intervals_pending_source_admission",
            "publication_ready": False, "uncertainty_admitted": False,
            "replicates": replicates, "seed": seed, "numpy_version": np.__version__,
            "rng": "numpy.PCG64 multinomial; shared family draws across three variants",
            "alpha": .05, "multiplicity_endpoints": 6, "quantile_method": "linear",
            "units": "raw 0-to-1 metric units", "families": report["families"],
            "point_estimates": {v: dict(zip(METRICS, row.tolist())) for v, row in zip(VARIANTS, observed)},
            "comparisons": comparisons,
            "limitations": ["Development-exposed exploratory analysis, not independent confirmation.",
                "Only 18 families; disjoint genes do not establish exchangeability.",
                "Adjustment covers these six endpoints only, not other QfO scores or their custom mean.",
                "Equal E-values do not establish matched sensitivity or cost.",
                "Numerical output alone does not verify source provenance."]}
