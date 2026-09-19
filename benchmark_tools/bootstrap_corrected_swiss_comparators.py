"""Corrected-release comparator bootstrap, retaining unavailable contrasts."""

import numpy as np

from benchmark_tools.audit_qfo_swiss_counts import LABELS, REFERENCE_SHA, statistics
from benchmark_tools.bootstrap_qfo_swiss_comparators import METHODS, CONTRASTS, METRICS
from benchmark_tools.bootstrap_qfo_swiss_stages import aggregate


def validated_values(report):
    if (report["status"] != "corrected_comparison_swiss_counts_verified"
            or report["shared_represented_genes"] or report["reference"]["sha256"] != REFERENCE_SHA
            or report["publication_ready"] is not False or report["uncertainty_admitted"] is not False):
        raise ValueError("Require audited disjoint corrected comparator counts")
    families = report["families"]
    if len(families) != 18 or len(set(families)) != 18 or tuple(r["method"] for r in report["methods"]) != METHODS:
        raise ValueError("Changed method or family inventory")
    if type(report["reference_relation_count"]) is not int or report["reference_relation_count"] <= 0:
        raise ValueError("Invalid reference relation count")
    values, universe = {}, None
    for method in report["methods"]:
        if method["status"] == "not_admitted":
            if set(method) != {"method", "status", "reason"} or not isinstance(method["reason"], str) or not method["reason"].strip():
                raise ValueError("Missing method must state a reason without imputed evidence")
            continue
        if method["status"] != "counts_verified" or [r["family"] for r in method["families"]] != families:
            raise ValueError("Wrong method status or family order")
        current, reference, seen = [], [], set()
        for row in method["families"]:
            counts, genes = row["counts_without_prior"], row["represented_genes"]
            if set(counts) != set(LABELS) or any(type(v) is not int or v < 0 for v in counts.values()):
                raise ValueError("Invalid raw confusion counts")
            if not sum(counts.values()) or len(genes) <= 5 or len(set(genes)) != len(genes) or seen.intersection(genes):
                raise ValueError("Empty or overlapping family")
            seen.update(genes)
            reference.append((tuple(sorted(genes)), counts["TP"] + counts["FN"], counts["FP"] + counts["TN"]))
            scores, stored = statistics(counts), row["statistics_with_prior"]
            if set(stored) != set(METRICS) or any(type(stored[m]) not in (int, float)
                    or not np.isfinite(stored[m]) or abs(scores[m] - stored[m]) > 1e-12 for m in METRICS):
                raise ValueError("Stored family statistic differs from raw counts")
            current.append([scores["PPV"], scores["TPR"]])
        if universe is None:
            universe = reference
        if reference != universe or sum(r[1] + r[2] for r in reference) != report["reference_relation_count"]:
            raise ValueError("Reference genes or truth totals differ")
        array = np.asarray(current)
        point, stored = aggregate(array.mean(axis=0)), method["aggregate"]
        if set(stored) != set(METRICS) or any(type(stored[m]) not in (int, float)
                or not np.isfinite(stored[m]) or abs(point[j] - stored[m]) > 1e-12 for j, m in enumerate(METRICS)):
            raise ValueError("Aggregate is not the native macro statistic")
        values[method["method"]] = array
    return values


def bootstrap(report, replicates=100000, seed=20260920):
    if type(replicates) is not int or replicates < 100 or type(seed) is not int or seed < 0:
        raise ValueError("Invalid bootstrap controls")
    values = validated_values(report)
    weights = np.random.Generator(np.random.PCG64(seed)).multinomial(18, np.full(18, 1 / 18), size=replicates)
    draws = {name: aggregate(weights @ cell / 18) for name, cell in values.items()}
    points = {name: aggregate(cell.mean(axis=0)) for name, cell in values.items()}
    comparisons = []
    for candidate_index, reference_index in CONTRASTS:
        candidate, reference = METHODS[candidate_index], METHODS[reference_index]
        row = {"candidate": candidate, "reference": reference}
        missing = [name for name in (candidate, reference) if name not in values]
        if missing:
            row.update(status="not_estimable", unavailable_methods=missing, metrics=None, family_differences=None)
        else:
            differences = draws[candidate] - draws[reference]
            family_differences = aggregate(values[candidate]) - aggregate(values[reference])
            metrics = {}
            for j, metric in enumerate(METRICS):
                metrics[metric] = {"difference": float(points[candidate][j] - points[reference][j]),
                    "paired_percentile_ci": np.quantile(differences[:, j], [.025, .975], method="linear").tolist(),
                    "bonferroni_percentile_ci": np.quantile(differences[:, j], [.05 / 48, 1 - .05 / 48], method="linear").tolist(),
                    "family_wins": int(np.sum(family_differences[:, j] > 1e-10)),
                    "family_ties": int(np.sum(np.abs(family_differences[:, j]) <= 1e-10)),
                    "family_losses": int(np.sum(family_differences[:, j] < -1e-10))}
            row.update(status="estimated", metrics=metrics, family_differences=[
                dict(family=name, **dict(zip(METRICS, delta.tolist())))
                for name, delta in zip(report["families"], family_differences)])
        comparisons.append(row)
    return {"status": "corrected_comparator_intervals_pending_provenance_admission",
        "publication_ready": False, "uncertainty_admitted": False,
        "replicates": replicates, "seed": seed, "alpha": .05, "multiplicity_endpoints": 24,
        "protocol_controls_match": replicates == 100000 and seed == 20260920,
        "numpy_version": np.__version__, "quantile_method": "linear",
        "rng": "numpy.PCG64 multinomial; shared family multiplicities across all available methods",
        "units": "raw 0-to-1 differences", "families": report["families"],
        "point_estimates": {name: dict(zip(METRICS, points[name].tolist())) if name in points else None for name in METHODS},
        "comparisons": comparisons, "limitations": [
            "Retrospective, development-exposed conditional family analysis; not independent or selection-adjusted confirmation.",
            "Only 18 families; shared evolutionary history and merged predictions can violate exchangeability.",
            "All eight contrasts and 24 planned endpoints are retained even when methods are unavailable.",
            "Unavailable results are not zeros. Intervals containing zero do not establish equivalence.",
            "Sequence-only OrthoFinder is an MCL-checkpoint diagnostic; FastOMA uses a supplied OrthoFinder tree.",
            "Phylogenetic versus sensitive OrthoHMM also changes candidate expansion, not only reconciliation.",
            "No transfer to historical inputs, other QfO challenges, the secondary mean or global superiority."]}
