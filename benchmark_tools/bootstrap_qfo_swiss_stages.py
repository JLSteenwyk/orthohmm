"""Paired family resampling for the frozen four-stage SwissTrees analysis."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_qfo_swiss_counts import LABELS, statistics
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.report_qfo_recovered_stages import STAGES, CONTRASTS

COUNTS_SHA = "546bb5bd6957c8ea990324b79fa31f22b0ed721bc7d6b94b609ab15258f97183"
PROTOCOL_SHA = "427b3079ef69d2f8b16ceb078f815c8e62b4a182f9402916fb5cba1471cf5ebe"
METRICS = ("F1", "PPV", "TPR")


def aggregate(means):
    precision, recall = means[..., 0], means[..., 1]
    return np.stack((2 * precision * recall / (precision + recall), precision, recall), axis=-1)


def validated_values(report):
    if report["status"] != "raw_swiss_family_counts_verified" or report["shared_represented_genes"]:
        raise ValueError("Require verified disjoint family counts")
    families = report["families"]
    if len(families) != 18 or len(set(families)) != 18:
        raise ValueError("Require all 18 unique families")
    if [r["stage"] for r in report["stages"]] != list(STAGES):
        raise ValueError("Changed stage inventory")
    values, reference = [], None
    for stage in report["stages"]:
        if [r["family"] for r in stage["families"]] != families:
            raise ValueError("Changed family inventory/order")
        current, universe, seen, stats = [], [], set(), []
        for row in stage["families"]:
            counts, genes = row["counts_without_prior"], row["represented_genes"]
            if set(counts) != set(LABELS) or any(type(v) is not int or v < 0 for v in counts.values()):
                raise ValueError("Invalid raw confusion counts")
            if not sum(counts.values()) or len(genes) <= 5 or len(set(genes)) != len(genes) or seen.intersection(genes):
                raise ValueError("Empty family, duplicate genes or overlapping families")
            seen.update(genes)
            truth_totals = (counts["TP"] + counts["FN"], counts["FP"] + counts["TN"])
            universe.append((tuple(sorted(genes)), truth_totals))
            calculated = statistics(counts)
            stored = row["statistics_with_prior"]
            if set(stored) != set(METRICS) or any(
                    not np.isfinite(stored[m]) or abs(calculated[m] - stored[m]) > 1e-12 for m in METRICS):
                raise ValueError("Stored family statistics disagree with counts")
            current.append([calculated["PPV"], calculated["TPR"]])
            stats.append(sum(counts.values()))
        if reference is None:
            reference = universe
        if universe != reference or sum(stats) != report["reference_relation_count"]:
            raise ValueError("Reference members or truth totals differ")
        array = np.asarray(current)
        point = aggregate(array.mean(axis=0))
        if any(not np.isfinite(stage["aggregate"][m]) or abs(point[j] - stage["aggregate"][m]) > 1e-12
               for j, m in enumerate(METRICS)):
            raise ValueError("Stored aggregate disagrees with native macro statistic")
        values.append(array)
    return np.asarray(values)


def bootstrap(report, replicates=100000, seed=20260919):
    if type(replicates) is not int or replicates < 100 or type(seed) is not int or seed < 0:
        raise ValueError("Invalid bootstrap controls")
    values = validated_values(report)
    n = values.shape[1]
    multiplicities = np.random.Generator(np.random.PCG64(seed)).multinomial(
        n, np.full(n, 1 / n), size=replicates)
    draws = np.asarray([aggregate(multiplicities @ stage / n) for stage in values])
    observed = aggregate(values.mean(axis=1))
    family_statistics = aggregate(values)
    comparisons = []
    for candidate, baseline in CONTRASTS:
        differences = draws[candidate] - draws[baseline]
        family_differences = family_statistics[candidate] - family_statistics[baseline]
        metrics = {}
        for j, metric in enumerate(METRICS):
            metrics[metric] = {
                "difference": float(observed[candidate, j] - observed[baseline, j]),
                "paired_percentile_ci": np.quantile(differences[:, j], [0.025, 0.975], method="linear").tolist(),
                "bonferroni_percentile_ci": np.quantile(differences[:, j], [0.05 / 24, 1 - 0.05 / 24], method="linear").tolist(),
                "family_wins": int(np.sum(family_differences[:, j] > 1e-10)),
                "family_ties": int(np.sum(np.abs(family_differences[:, j]) <= 1e-10)),
                "family_losses": int(np.sum(family_differences[:, j] < -1e-10)),
            }
        comparisons.append({"candidate": STAGES[candidate], "reference": STAGES[baseline], "metrics": metrics,
                            "family_differences": [dict(family=name, **dict(zip(METRICS, row.tolist())))
                                                   for name, row in zip(report["families"], family_differences)]})
    return {"status": "paired_swiss_stage_intervals", "publication_ready": False,
            "replicates": replicates, "seed": seed, "alpha": 0.05, "multiplicity_endpoints": 12,
            "rng": "numpy.PCG64 multinomial; shared family draws across all stages",
            "numpy_version": np.__version__, "quantile_method": "linear", "units": "raw 0-to-1 metric units",
            "families": report["families"], "point_estimates": {stage: dict(zip(METRICS, row.tolist()))
                                                                 for stage, row in zip(STAGES, observed)},
            "comparisons": comparisons,
            "limitations": ["Development-exposed follow-up, not independent or selection-adjusted confirmation.",
                            "Only 18 curated families; shared history and merged predictions can violate exchangeability.",
                            "Percentile intervals are approximate; adjustment covers these 12 endpoints only.",
                            "Family wins/ties/losses and family differences are descriptive.",
                            "Profiles include downstream graph/assignment responses; refined is not phylogeny.",
                            "These intervals do not cover other QfO challenges or the six-metric secondary mean."]}


def markdown(report):
    lines = ["# Recovered SwissTrees Paired Intervals", "",
             "18 families; 100000 shared draws; seed20260919. Raw0-to1 differences, candidate minus reference.",
             "Development-exposed follow-up. Adjustment covers all12 planned endpoints.", "",
             "| Contrast | Metric | Difference | Nominal95% CI | Adjusted CI | Wins/ties/losses |",
             "| --- | --- | ---: | --- | --- | --- |"]
    for row in report["comparisons"]:
        for metric in METRICS:
            m = row["metrics"][metric]
            intervals = ["[%.6f, %.6f]" % tuple(m[k]) for k in ("paired_percentile_ci", "bonferroni_percentile_ci")]
            lines.append(f"| {row['candidate']} - {row['reference']} | {metric} | {m['difference']:+.6f} | " +
                         " | ".join(intervals) + f" | {m['family_wins']}/{m['family_ties']}/{m['family_losses']} |")
    return "\n".join([*lines, "", *["- " + text for text in report["limitations"]], ""])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("counts", "protocol", "output", "markdown"):
        parser.add_argument("--" + flag, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.markdown.exists():
        raise FileExistsError("Require fresh result paths")
    counts_record, protocol_record = record(args.counts), record(args.protocol)
    if counts_record["sha256"] != COUNTS_SHA or protocol_record["sha256"] != PROTOCOL_SHA:
        raise ValueError("Changed frozen counts/protocol")
    result = bootstrap(json.loads(args.counts.read_text()))
    check(counts_record)
    check(protocol_record)
    result.update(counts=counts_record, protocol=protocol_record, source=record(Path(__file__)),
                  helpers=[record(Path(__file__).with_name(name)) for name in
                           ("audit_qfo_swiss_counts.py", "report_qfo_recovered_stages.py")])
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    with args.markdown.open("x") as handle:
        handle.write(markdown(result))
