"""Frozen 42-endpoint paired SwissTrees analysis for the eight-cell QfO factorial."""

import argparse
import itertools
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_qfo_swiss_counts import LABELS, REFERENCE_SHA, statistics
from benchmark_tools.bootstrap_qfo_swiss_stages import aggregate, METRICS
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check

CELLS = tuple(f"p{p}_c{c}_r{r}" for p, c, r in itertools.product((0, 1), repeat=3))
PROTOCOL_SHA = "f8946e12cefcf84abbee0fb9492f240c05508e045efe00a3006304d34c1fd115"


def contrasts():
    rows = []
    factors = ("P", "C", "R")
    for axis, factor in enumerate(factors):
        others = [i for i in range(3) if i != axis]
        for levels in itertools.product((0, 1), repeat=2):
            bits = [0, 0, 0]
            for i, level in zip(others, levels):
                bits[i] = level
            baseline = 4 * bits[0] + 2 * bits[1] + bits[2]
            candidate = baseline + (4, 2, 1)[axis]
            weights = [0] * 8
            weights[candidate], weights[baseline] = 1, -1
            rows.append({"name": f"{factor}_at_" + "_".join(f"{factors[i]}{level}" for i, level in zip(others, levels)),
                         "kind": "simple_effect", "weights": weights,
                         "candidate": CELLS[candidate], "reference": CELLS[baseline]})
    for p in (0, 1):
        weights = [0] * 8
        for offset, value in ((0, 1), (1, -1), (2, -1), (3, 1)):
            weights[4 * p + offset] = value
        rows.append({"name": f"C_by_R_at_P{p}", "kind": "interaction", "weights": weights,
                     "definition": "(C1R1-C0R1) - (C1R0-C0R0) at fixed P"})
    return rows


def validated_values(report):
    if report["status"] != "qfo_factorial_swiss_counts_verified" or report["shared_represented_genes"]:
        raise ValueError("Require verified disjoint SwissTrees family counts")
    if report["reference"]["sha256"] != REFERENCE_SHA:
        raise ValueError("Changed SwissTrees reference")
    families = report["families"]
    if len(families) != 18 or len(set(families)) != 18 or [r["cell"] for r in report["cells"]] != list(CELLS):
        raise ValueError("Require all 18 families and eight cells in frozen order")
    values, reference = [], None
    for cell in report["cells"]:
        if [r["family"] for r in cell["families"]] != families:
            raise ValueError("Changed family order")
        seen, universe, current, relations = set(), [], [], 0
        for row in cell["families"]:
            counts, genes = row["counts_without_prior"], row["represented_genes"]
            if set(counts) != set(LABELS) or any(type(v) is not int or v < 0 for v in counts.values()):
                raise ValueError("Invalid raw confusion counts")
            if not sum(counts.values()) or len(genes) <= 5 or len(set(genes)) != len(genes) or seen.intersection(genes):
                raise ValueError("Empty or overlapping family")
            seen.update(genes)
            universe.append((tuple(sorted(genes)), counts["TP"] + counts["FN"], counts["FP"] + counts["TN"]))
            calculated = statistics(counts)
            if set(row["statistics_with_prior"]) != set(METRICS) or any(
                not np.isfinite(row["statistics_with_prior"][m]) or abs(calculated[m] - row["statistics_with_prior"][m]) > 1e-12 for m in METRICS):
                raise ValueError("Stored family statistic differs from counts")
            current.append([calculated["PPV"], calculated["TPR"]])
            relations += sum(counts.values())
        if reference is None:
            reference = universe
        if universe != reference or relations != report["reference_relation_count"]:
            raise ValueError("Reference members or truth totals differ")
        values.append(current)
        point = aggregate(np.mean(current, axis=0))
        if any(not np.isfinite(cell["aggregate"][m]) or abs(point[j] - cell["aggregate"][m]) > 1e-12 for j, m in enumerate(METRICS)):
            raise ValueError("Stored aggregate differs from actual benchmark statistic")
    return np.asarray(values)


def bootstrap(report, replicates=100000, seed=20260922):
    if type(replicates) is not int or replicates < 100 or type(seed) is not int or seed < 0:
        raise ValueError("Invalid bootstrap controls")
    values = validated_values(report)
    n = len(report["families"])
    multiplicities = np.random.Generator(np.random.PCG64(seed)).multinomial(n, np.full(n, 1 / n), size=replicates)
    draws = np.asarray([aggregate(multiplicities @ cell / n) for cell in values])
    observed, families = aggregate(values.mean(axis=1)), aggregate(values)
    comparisons = []
    for contrast in contrasts():
        weights = np.asarray(contrast["weights"])
        delta = np.tensordot(weights, draws, axes=1)
        point = weights @ observed
        per_family = np.tensordot(weights, families, axes=1)
        metrics = {}
        for j, metric in enumerate(METRICS):
            metrics[metric] = {"difference": float(point[j]),
                "paired_percentile_ci": np.quantile(delta[:, j], [.025, .975], method="linear").tolist(),
                "bonferroni_percentile_ci": np.quantile(delta[:, j], [.05 / 84, 1 - .05 / 84], method="linear").tolist(),
                "family_wins": int(np.sum(per_family[:, j] > 1e-10)),
                "family_ties": int(np.sum(np.abs(per_family[:, j]) <= 1e-10)),
                "family_losses": int(np.sum(per_family[:, j] < -1e-10))}
        comparisons.append({**contrast, "metrics": metrics, "family_differences":
            [dict(family=name, **dict(zip(METRICS, row.tolist()))) for name, row in zip(report["families"], per_family)]})
    return {"status": "paired_qfo_factorial_swiss_intervals", "publication_ready": False,
            "replicates": replicates, "seed": seed, "rng": "numpy.PCG64 multinomial; shared family draws across all eight cells",
            "numpy_version": np.__version__, "alpha": .05, "multiplicity_endpoints": 42, "quantile_method": "linear",
            "units": "raw 0-to-1 metric units", "families": report["families"],
            "point_estimates": {cell: dict(zip(METRICS, row.tolist())) for cell, row in zip(CELLS, observed)},
            "comparisons": comparisons, "limitations": [
                "Development-exposed conditional analysis, not independent or selection-adjusted validation.",
                "Only 18 families; shared history and merged predictions may violate family exchangeability.",
                "Adjustment covers these 42 endpoints only; no intervals for other QfO challenges or the secondary mean.",
                "Profile-off retains initial HMM search; P includes downstream sequence refinement.",
                "R compares complete native inferred-pair strategy against group-derived pairs, not only group splitting.",
                "For interactions, family wins/losses denote positive/negative interaction, not a method ranking."]}


def markdown(report):
    lines = ["# QfO Factorial SwissTrees Intervals", "", f"{report['replicates']} shared draws; seed {report['seed']}; 42 adjusted endpoints.", "",
             "| Contrast | Metric | Difference | Nominal 95% CI | Adjusted CI | Positive/tie/negative families |",
             "|---|---|---:|---|---|---|"]
    for row in report["comparisons"]:
        for name, metric in row["metrics"].items():
            intervals = ["[%.6f, %.6f]" % tuple(metric[key]) for key in ("paired_percentile_ci", "bonferroni_percentile_ci")]
            lines.append(f"| {row['name']} | {name} | {metric['difference']:+.6f} | " + " | ".join(intervals) +
                         f" | {metric['family_wins']}/{metric['family_ties']}/{metric['family_losses']} |")
    return "\n".join([*lines, "", *["- " + s for s in report["limitations"]], ""])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("counts", "protocol", "output", "markdown"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--counts-sha256", required=True)
    args = parser.parse_args()
    if args.output.exists() or args.markdown.exists():
        raise FileExistsError("Require fresh result paths")
    counts, protocol = record(args.counts), record(args.protocol)
    if counts["sha256"] != args.counts_sha256 or protocol["sha256"] != PROTOCOL_SHA:
        raise ValueError("Changed counts or frozen protocol")
    data = json.loads(args.counts.read_text())
    for item in data["checked_inputs"]:
        check(item)
    result = bootstrap(data)
    for item in [counts, protocol, *data["checked_inputs"]]:
        check(item)
    result.update(counts=counts, protocol=protocol, source=record(__file__), helpers=[record(Path(__file__).with_name(n))
        for n in ("audit_qfo_swiss_counts.py", "bootstrap_qfo_swiss_stages.py")])
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    with args.markdown.open("x") as stream:
        stream.write(markdown(result))
