"""Numerical kernel for the prespecified corrected SwissTrees strata analysis.

This module does not admit scientific inputs. A caller must independently bind
corrected prediction/count evidence and frozen sequence-stratum membership.
"""

import numpy as np

from benchmark_tools.audit_qfo_swiss_counts import LABELS, statistics
from benchmark_tools.bootstrap_qfo_swiss_stages import METRICS, aggregate

METHODS = ("high_sensitivity", "phylogenetic", "orthofinder_full")
CONTRASTS = ((0, 2), (1, 2), (1, 0))
BINS = ("lower", "higher", "missing")
ENDPOINTS = 27


def validated_values(counts, membership):
    if set(counts) != set(METHODS):
        raise ValueError("Require the three prespecified methods")
    families = list(membership)
    if len(families) != 18 or any(not isinstance(f, str) or not f for f in families):
        raise ValueError("Require 18 named reference families")
    if any(b not in BINS for b in membership.values()):
        raise ValueError("Unknown primary bin")
    values, reference = [], None
    for method in METHODS:
        rows = counts[method]
        if set(rows) != set(families):
            raise ValueError("Changed family inventory")
        current, universe, seen = [], [], set()
        for family in families:
            row = rows[family]
            raw, genes = row["counts_without_prior"], row["represented_genes"]
            if set(raw) != set(LABELS) or any(type(v) is not int or v < 0 for v in raw.values()):
                raise ValueError("Invalid raw confusion counts")
            if (not sum(raw.values()) or len(genes) <= 5
                    or any(not isinstance(g, str) or not g for g in genes)
                    or len(set(genes)) != len(genes) or seen.intersection(genes)):
                raise ValueError("Empty, duplicate or overlapping reference members")
            seen.update(genes)
            universe.append((tuple(sorted(genes)), raw["TP"] + raw["FN"], raw["FP"] + raw["TN"]))
            scores = statistics(raw)
            current.append([scores["PPV"], scores["TPR"]])
        if reference is None:
            reference = universe
        if universe != reference:
            raise ValueError("Reference members or truth totals differ")
        values.append(current)
    return families, np.asarray(values)


def _intervals(draws):
    if draws is None:
        return {"paired_percentile_ci": None, "bonferroni_percentile_ci": None}
    return {
        "paired_percentile_ci": np.quantile(draws, [.025, .975], method="linear").tolist(),
        "bonferroni_percentile_ci": np.quantile(
            draws, [.05 / (2 * ENDPOINTS), 1 - .05 / (2 * ENDPOINTS)], method="linear").tolist(),
    }


def bootstrap(counts, membership, replicates=100000, seed=20260924):
    """Calculate paired contrasts; provenance admission is deliberately external."""
    if type(replicates) is not int or replicates < 100 or type(seed) is not int or seed < 0:
        raise ValueError("Invalid bootstrap controls")
    families, values = validated_values(counts, membership)
    rng = np.random.Generator(np.random.PCG64(seed))
    bins, point, draws = {}, {}, {}
    for name in BINS:
        indices = [i for i, f in enumerate(families) if membership[f] == name]
        n = len(indices)
        selected = values[:, indices, :]
        point[name] = aggregate(selected.mean(axis=1)) if n else None
        draws[name] = None
        if name != "missing" and n >= 5:
            weights = rng.multinomial(n, np.full(n, 1 / n), size=replicates)
            draws[name] = np.asarray([aggregate(weights @ method / n) for method in selected])
        comparisons = []
        for candidate, reference in CONTRASTS:
            family_differences = aggregate(selected[candidate]) - aggregate(selected[reference])
            metrics = {}
            for j, metric in enumerate(METRICS):
                difference = None if not n else float(point[name][candidate, j] - point[name][reference, j])
                differences = None if draws[name] is None else draws[name][candidate, :, j] - draws[name][reference, :, j]
                metrics[metric] = dict(
                    difference=difference, **_intervals(differences),
                    family_wins=int(np.sum(family_differences[:, j] > 1e-10)),
                    family_ties=int(np.sum(np.abs(family_differences[:, j]) <= 1e-10)),
                    family_losses=int(np.sum(family_differences[:, j] < -1e-10)))
            comparisons.append(dict(candidate=METHODS[candidate], reference=METHODS[reference], metrics=metrics))
        bins[name] = dict(families=[families[i] for i in indices], interval_eligible=draws[name] is not None,
                          point_estimates={method: None if not n else dict(zip(METRICS, point[name][i].tolist()))
                                           for i, method in enumerate(METHODS)}, comparisons=comparisons)
    interactions = []
    for candidate, reference in CONTRASTS:
        metrics = {}
        for j, metric in enumerate(METRICS):
            difference, differences = None, None
            if point["lower"] is not None and point["higher"] is not None:
                difference = float((point["higher"][candidate, j] - point["higher"][reference, j])
                                   - (point["lower"][candidate, j] - point["lower"][reference, j]))
            if draws["lower"] is not None and draws["higher"] is not None:
                differences = ((draws["higher"][candidate, :, j] - draws["higher"][reference, :, j])
                               - (draws["lower"][candidate, :, j] - draws["lower"][reference, :, j]))
            metrics[metric] = dict(difference=difference, **_intervals(differences))
        interactions.append(dict(candidate=METHODS[candidate], reference=METHODS[reference], metrics=metrics))
    return dict(status="numerical_strata_result_pending_provenance_admission", publication_ready=False,
                scientific_inputs_admitted=False, replicates=replicates, seed=seed,
                multiplicity_endpoints=ENDPOINTS, alpha=.05, numpy_version=np.__version__,
                quantile_method="linear", rng="PCG64 multinomial; lower then higher; paired across methods",
                units="raw 0-to-1 differences", bins=bins, interactions=interactions,
                interaction_direction="higher contrast minus lower contrast",
                limitations=["No prediction or input provenance is admitted by this numerical kernel.",
                             "Development-exposed, conditional approximate family bootstrap; not independent validation.",
                             "Configuration contrasts are not pure phylogeny ablations; strata are not causal explanations.",
                             "Adjustment covers the 27 primary endpoints, not prior development or other QfO metrics."])
