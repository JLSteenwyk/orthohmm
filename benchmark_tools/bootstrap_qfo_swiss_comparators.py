"""Frozen paired family bootstrap for the eight historical SwissTrees comparators."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_qfo_swiss_comparators import METHODS
from benchmark_tools.audit_qfo_swiss_counts import LABELS, statistics
from benchmark_tools.bootstrap_qfo_swiss_stages import aggregate, METRICS
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check

COUNTS_SHA = "2995868b0407ceda2e7422db6c3fc3b99523716e4bec8b1b43e15f296e40769b"
PROTOCOL_SHA = "a1ba79268b732824e0df9c56dc1b14dce7b9f33c955523e8b87133c5440161b8"
CONTRASTS = ((0, 2), (1, 2), (3, 2), (4, 2), (5, 2), (6, 2), (7, 2), (1, 0))


def validated_values(report):
    if report["status"] != "eight_comparator_swiss_family_counts_verified" or report["shared_represented_genes"]:
        raise ValueError("Require verified disjoint family counts")
    families = report["families"]
    if len(families) != 18 or len(set(families)) != 18:
        raise ValueError("Require 18 distinct families")
    if tuple(r["method"] for r in report["methods"]) != METHODS:
        raise ValueError("Changed method inventory")
    values, reference = [], None
    for method in report["methods"]:
        if [r["family"] for r in method["families"]] != families:
            raise ValueError("Changed family inventory")
        current, universe, seen = [], [], set()
        for row in method["families"]:
            counts, genes = row["counts_without_prior"], row["represented_genes"]
            if set(counts) != set(LABELS) or any(type(v) is not int or v < 0 for v in counts.values()):
                raise ValueError("Invalid confusion counts")
            if not sum(counts.values()) or len(genes) <= 5 or len(set(genes)) != len(genes) or seen.intersection(genes):
                raise ValueError("Empty or overlapping family")
            seen.update(genes)
            universe.append((tuple(sorted(genes)), counts["TP"] + counts["FN"], counts["FP"] + counts["TN"]))
            scores = statistics(counts)
            stored = row["statistics_with_prior"]
            if set(stored) != set(METRICS) or any(not np.isfinite(stored[m]) or abs(scores[m] - stored[m]) > 1e-12 for m in METRICS):
                raise ValueError("Stored family statistics disagree")
            current.append([scores["PPV"], scores["TPR"]])
        if reference is None:
            reference = universe
        if universe != reference or sum(r[1] + r[2] for r in universe) != report["reference_relation_count"]:
            raise ValueError("Reference members or truth totals differ")
        current = np.asarray(current)
        point = aggregate(current.mean(axis=0))
        if any(not np.isfinite(method["aggregate"][m]) or abs(point[j] - method["aggregate"][m]) > 1e-12
               for j, m in enumerate(METRICS)):
            raise ValueError("Stored aggregate disagrees")
        values.append(current)
    return np.asarray(values)


def bootstrap(report, replicates=100000, seed=20260920):
    if type(replicates) is not int or replicates < 100 or type(seed) is not int or seed < 0:
        raise ValueError("Invalid bootstrap controls")
    values = validated_values(report)
    n = len(report["families"])
    weights = np.random.Generator(np.random.PCG64(seed)).multinomial(n, np.full(n, 1 / n), size=replicates)
    draws = np.asarray([aggregate(weights @ method / n) for method in values])
    observed, family_scores = aggregate(values.mean(axis=1)), aggregate(values)
    comparisons = []
    endpoints = len(CONTRASTS) * len(METRICS)
    for candidate, reference in CONTRASTS:
        differences = draws[candidate] - draws[reference]
        family_differences = family_scores[candidate] - family_scores[reference]
        metrics = {}
        for j, metric in enumerate(METRICS):
            metrics[metric] = {
                "difference": float(observed[candidate, j] - observed[reference, j]),
                "paired_percentile_ci": np.quantile(differences[:, j], [0.025, 0.975], method="linear").tolist(),
                "bonferroni_percentile_ci": np.quantile(differences[:, j], [0.05 / (2 * endpoints), 1 - 0.05 / (2 * endpoints)], method="linear").tolist(),
                "family_wins": int(np.sum(family_differences[:, j] > 1e-10)),
                "family_ties": int(np.sum(np.abs(family_differences[:, j]) <= 1e-10)),
                "family_losses": int(np.sum(family_differences[:, j] < -1e-10)),
            }
        comparisons.append({"candidate": METHODS[candidate], "reference": METHODS[reference], "metrics": metrics,
                            "family_differences": [dict(family=name, **dict(zip(METRICS, row.tolist())))
                                                   for name, row in zip(report["families"], family_differences)]})
    return {"status": "paired_swiss_comparator_intervals", "publication_ready": False,
            "replicates": replicates, "seed": seed, "alpha": 0.05, "multiplicity_endpoints": endpoints,
            "rng": "numpy.PCG64 multinomial; shared family multiplicities across all methods",
            "numpy_version": np.__version__, "quantile_method": "linear", "units": "raw 0-to-1 differences",
            "families": report["families"], "point_estimates": {name: dict(zip(METRICS, row.tolist())) for name, row in zip(METHODS, observed)},
            "comparisons": comparisons,
            "limitations": ["Retrospective, development-exposed; not independent or model-selection-adjusted confirmation.",
                            "Only 18 curated families; shared history and merged predictions can violate exchangeability.",
                            "Approximate conditional family-bootstrap sensitivity analysis; adjustment covers these 24 endpoints only.",
                            "Family wins/ties/losses are descriptive; inclusion of zero does not establish equivalence.",
                            "Sequence-only OrthoFinder is an MCL-checkpoint diagnostic; FastOMA used a supplied OrthoFinder tree.",
                            "Phylogenetic versus sensitive OrthoHMM is not a pure reconciliation ablation.",
                            "No interval transfer to recovered ablations, other QfO challenges or the secondary six-metric mean."]}


def markdown(report):
    lines = ["# SwissTrees Comparator Paired Intervals", "",
             "18 families; 100000 shared draws; seed20260920. Candidate minus reference in raw0-to1 units.",
             "Adjustment covers24endpoints. Historical development-exposed comparison, not independent confirmation.", "",
             "| Candidate | Reference | Metric | Difference | Nominal95% CI | Adjusted CI | Wins/ties/losses |",
             "| --- | --- | --- | ---: | --- | --- | --- |"]
    for row in report["comparisons"]:
        for metric in METRICS:
            m = row["metrics"][metric]
            intervals = ["[%.6f, %.6f]" % tuple(m[k]) for k in ("paired_percentile_ci", "bonferroni_percentile_ci")]
            lines.append(f"| {row['candidate']} | {row['reference']} | {metric} | {m['difference']:+.6f} | " +
                         " | ".join(intervals) + f" | {m['family_wins']}/{m['family_ties']}/{m['family_losses']} |")
    return "\n".join([*lines, "", *["- " + text for text in report["limitations"]], ""])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("counts", "protocol", "output", "markdown"):
        parser.add_argument("--" + flag, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.markdown.exists():
        raise FileExistsError("Require fresh output paths")
    counts_identity, protocol_identity = record(args.counts), record(args.protocol)
    if counts_identity["sha256"] != COUNTS_SHA or protocol_identity["sha256"] != PROTOCOL_SHA:
        raise ValueError("Changed frozen counts or protocol")
    result = bootstrap(json.loads(args.counts.read_text()))
    check(counts_identity)
    check(protocol_identity)
    result.update(counts=counts_identity, protocol=protocol_identity, source=record(__file__),
                  helpers=[record(Path(__file__).with_name(name)) for name in
                           ("audit_qfo_swiss_counts.py", "audit_qfo_swiss_comparators.py", "bootstrap_qfo_swiss_stages.py")])
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    args.markdown.write_text(markdown(result))
