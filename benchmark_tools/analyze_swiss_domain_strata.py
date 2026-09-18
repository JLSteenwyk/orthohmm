"""Execute the frozen annotation-defined SwissTrees error-strata analysis."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.bootstrap_qfo_swiss_comparators import validated_values, aggregate, METHODS, METRICS, COUNTS_SHA
from benchmark_tools.inventory_swiss_annotations import summarize
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check

ANNOTATION_SHA = "d5e269158c2fb603acc0805a1c140a88758342d7b8e32b3094750ade22fbc06c"
PROTOCOL_SHA = "8316866d11f988cef6a3e88f07e2802d1557a7ae581b100b85884f704b0afee8"
CONTRASTS = ((0, 2), (1, 2), (1, 0))
PRIMARY = ("median_pfam_types_below_two", "median_pfam_types_at_least_two")
SECONDARY = ("repeated_type_fraction_below_quarter", "repeated_type_fraction_at_least_quarter")


def strata(counts, annotations):
    if annotations["status"] != "prediction_independent_swiss_annotation_inventory":
        raise ValueError("Unverified annotation inventory")
    rows = counts["methods"][0]["families"]
    universe = {g for row in rows for g in row["represented_genes"]}
    if set(annotations["genes"]) != universe or set(annotations["families"]) != set(counts["families"]):
        raise ValueError("Changed annotation coverage or families")
    bins = {key: [] for key in (*PRIMARY, *SECONDARY)}
    for row in rows:
        family, genes = row["family"], row["represented_genes"]
        summary = summarize(genes, annotations["genes"])
        if summary != annotations["families"][family] or summary["missing_annotation_genes"]:
            raise ValueError("Annotation summary differs from selected genes")
        bins[PRIMARY[int(summary["median_pfam_types_among_annotated"] >= 2)]].append(family)
        bins[SECONDARY[int(summary["annotated_repeated_pfam_type"] / len(genes) >= .25)]].append(family)
    if [len(bins[k]) for k in (*PRIMARY, *SECONDARY)] != [12, 6, 15, 3]:
        raise ValueError("Frozen stratum sizes changed")
    if set(bins[PRIMARY[1]]) != {"APP", "BAR", "HOX", "NOX", "TRFE", "VATB"} or set(bins[SECONDARY[1]]) != {"MAPT", "PSEN", "TRFE"}:
        raise ValueError("Frozen stratum membership changed")
    return bins


def interval_metrics(point, draws):
    return {metric: {"difference": float(point[j]),
                     "nominal95": np.quantile(draws[:, j], [.025, .975], method="linear").tolist(),
                     "bonferroni27": np.quantile(draws[:, j], [.05/54, 1-.05/54], method="linear").tolist()}
            for j, metric in enumerate(METRICS)}


def analyze(counts, annotations, replicates=100000, seed=20260921):
    if type(replicates) is not int or replicates < 100 or type(seed) is not int or seed < 0:
        raise ValueError("Invalid resampling controls")
    values = validated_values(counts)
    bins = strata(counts, annotations)
    rng = np.random.Generator(np.random.PCG64(seed))
    point, draws, selected = {}, {}, {}
    for name, families in bins.items():
        indices = [counts["families"].index(f) for f in families]
        selected[name] = values[:, indices, :]
        point[name] = aggregate(selected[name].mean(axis=1))
        if name in PRIMARY:
            n = len(families)
            weights = rng.multinomial(n, np.full(n, 1/n), size=replicates)
            draws[name] = np.asarray([aggregate(weights @ method / n) for method in selected[name]])
    contrasts = []
    for candidate, reference in CONTRASTS:
        effects, observed = {}, {}
        records = []
        for name in PRIMARY:
            observed[name] = point[name][candidate] - point[name][reference]
            effects[name] = draws[name][candidate] - draws[name][reference]
            family_delta = aggregate(selected[name][candidate]) - aggregate(selected[name][reference])
            records.append({"stratum": name, "metrics": interval_metrics(observed[name], effects[name]),
                            "family_wins_ties_losses": {metric: [int(np.sum(family_delta[:, j] > 1e-10)),
                                int(np.sum(np.abs(family_delta[:, j]) <= 1e-10)), int(np.sum(family_delta[:, j] < -1e-10))]
                                for j, metric in enumerate(METRICS)}})
        interaction = interval_metrics(observed[PRIMARY[1]] - observed[PRIMARY[0]], effects[PRIMARY[1]] - effects[PRIMARY[0]])
        family_delta = aggregate(values[candidate]) - aggregate(values[reference])
        contrasts.append({"candidate": METHODS[candidate], "reference": METHODS[reference], "primary_strata": records,
                          "interaction_higher_minus_lower": interaction,
                          "secondary_descriptive_differences": {name: dict(zip(METRICS, (point[name][candidate] - point[name][reference]).tolist())) for name in SECONDARY},
                          "family_differences": [dict(family=f, **dict(zip(METRICS, row.tolist()))) for f, row in zip(counts["families"], family_delta)]})
    return {"status": "frozen_swiss_domain_strata_evaluated", "replicates": replicates, "seed": seed,
            "multiplicity_endpoints": 27, "rng": "NumPy PCG64, independent primary-bin multinomial draws; paired across methods",
            "numpy_version": np.__version__, "units": "raw 0-to-1 differences", "strata": bins,
            "point_estimates": {name: {method: dict(zip(METRICS, row.tolist())) for method, row in zip(METHODS, matrix)} for name, matrix in point.items()},
            "contrasts": contrasts, "publication_ready": False,
            "limitations": ["Retrospective development-exposed explanatory analysis; no method retuning.",
                            "Six/twelve curated-family primary bins; approximate conditional bootstrap with possible family dependence.",
                            "Domain annotations do not establish causality, fragments, domain loss or ancestral duplication.",
                            "The repeated-type split is descriptive only; its higher bin contains three families.",
                            "Interactions are differences of method contrasts, not causal domain effects.",
                            "Adjustment covers these 27 endpoints, not prior tuning or all publication analyses.",
                            "MCL checkpoint and supplied-tree FastOMA remain diagnostics; OrthoHMM contrast is not a pure reconciliation ablation.",
                            "FAS annotations are not independent FAS validation; intervals do not cover other QfO metrics or the secondary mean."]}


def markdown(report):
    lines = ["# SwissTrees Domain-Stratified Results", "", "Raw0-to1 differences. Nominal95% and27-endpoint adjusted intervals.",
             "Primary bins:12families with median Pfam types<2,6with median>=2. Repeat bins:15/3families, descriptive only.", "",
             "| Candidate | Reference | Stratum | Metric | Difference | Nominal CI | Adjusted CI |", "| --- | --- | --- | --- | ---: | --- | --- |"]
    for contrast in report["contrasts"]:
        rows = [(r["stratum"], r["metrics"]) for r in contrast["primary_strata"]]
        rows.append(("interaction: higher minus lower", contrast["interaction_higher_minus_lower"]))
        for name, metrics in rows:
            for metric, values in metrics.items():
                ci = ["[%.6f, %.6f]" % tuple(values[k]) for k in ("nominal95", "bonferroni27")]
                lines.append(f"| {contrast['candidate']} | {contrast['reference']} | {name} | {metric} | {values['difference']:+.6f} | " + " | ".join(ci) + " |")
    lines += ["", "## All Method Point Estimates", "", "| Stratum | Method | F1 | Precision | Recall |", "| --- | --- | ---: | ---: | ---: |"]
    for name, methods in report["point_estimates"].items():
        for method, values in methods.items():
            lines.append(f"| {name} | {method} | {values['F1']:.6f} | {values['PPV']:.6f} | {values['TPR']:.6f} |")
    return "\n".join([*lines, "", *["- " + text for text in report["limitations"]], ""])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("counts", "annotations", "protocol", "output", "markdown"):
        parser.add_argument("--" + flag, required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists() or args.markdown.exists():
        raise FileExistsError("Require new result paths")
    identities = [record(p) for p in (args.counts, args.annotations, args.protocol)]
    if [r["sha256"] for r in identities] != [COUNTS_SHA, ANNOTATION_SHA, PROTOCOL_SHA]:
        raise ValueError("Frozen input identities differ")
    report = analyze(json.loads(args.counts.read_text()), json.loads(args.annotations.read_text()))
    for item in identities:
        check(item)
    report.update(inputs=identities, source=record(__file__), helpers=[record(Path(__file__).with_name(p)) for p in
                  ("bootstrap_qfo_swiss_comparators.py", "bootstrap_qfo_swiss_stages.py", "inventory_swiss_annotations.py")])
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    args.markdown.write_text(markdown(report))
