"""Display retained GO/EC arithmetic without rescoring or changing endpoints."""

import argparse
import csv
from fractions import Fraction
from itertools import combinations
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

BASELINE = "orthofinder_3_1_5_full"
METHODS = {
    "orthohmm_high_sensitivity": "OrthoHMM sensitive",
    "orthohmm_phylogeny_satellite_v2": "OrthoHMM phylogenetic",
    "orthofinder_3_1_5_sequence_only": "OrthoFinder sequence-only",
    "sonicparanoid_2_0_9": "SonicParanoid",
    "proteinortho_6_3_6": "ProteinOrtho",
    "fastoma_0_3_5": "FastOMA",
    "orthomcl_1_4": "OrthoMCL",
}
METRICS = ("GO", "EC")
COMPONENTS = ("shared_with_original_denominators", "left_only", "negative_right_only")


def close(actual, expected):
    if type(actual) not in (int, float) or not math.isfinite(actual) or abs(actual - float(expected)) > 1e-12:
        raise ValueError("Retained arithmetic differs from integer score sums")


def counts(result):
    names = ("left_pairs", "right_pairs", "shared_pairs", "left_only_pairs", "right_only_pairs",
             "left_score_sum_millionths", "right_score_sum_millionths",
             "shared_left_sum_millionths", "shared_right_sum_millionths")
    if any(type(result[n]) is not int or result[n] < 0 for n in names):
        raise ValueError("Require nonnegative integer counts and sums")
    nl, nr, ns = (result[n] for n in names[:3])
    sl, sr, ssl, ssr = (result[n] for n in names[5:])
    if (not nl or not nr or ns > min(nl, nr)
            or result["left_only_pairs"] != nl - ns or result["right_only_pairs"] != nr - ns
            or not 0 <= ssl <= ns * 1000000 or not 0 <= ssr <= ns * 1000000
            or not 0 <= sl - ssl <= (nl - ns) * 1000000
            or not 0 <= sr - ssr <= (nr - ns) * 1000000):
        raise ValueError("Impossible scored-pair partition")
    # This figure is specifically for the retained identical-serialized-score panel.
    if (ssl != ssr or type(result["shared_pairs_with_different_serialized_scores"]) is not int
            or result["shared_pairs_with_different_serialized_scores"] != 0
            or result["maximum_shared_absolute_difference_millionths"] != (0 if ns else None)):
        raise ValueError("Shared serialized scores are not identical")
    ml, mr = Fraction(sl, nl * 1000000), Fraction(sr, nr * 1000000)
    components = (Fraction(ssl, nl * 1000000) - Fraction(ssr, nr * 1000000),
                  Fraction(sl - ssl, nl * 1000000), -Fraction(sr - ssr, nr * 1000000))
    if sum(components) != ml - mr:
        raise ValueError("Exact decomposition fails")
    for key, value in (("left_mean", ml), ("right_mean", mr), ("original_mean_difference", ml - mr)):
        close(result[key], value)
    if set(result["original_mean_difference_components"]) != set(COMPONENTS):
        raise ValueError("Wrong component inventory")
    for key, value in zip(COMPONENTS, components):
        close(result["original_mean_difference_components"][key], value)
    if ns:
        close(result["shared_conditional_mean_difference"], 0)
    elif result["shared_conditional_mean_difference"] is not None:
        raise ValueError("Empty intersection must have no conditional mean")
    return nl, nr, ns, sl, sr, ssl, ssr


def endpoints(panel):
    if (panel["status"] != "corrected_all_method_scored_pair_panel"
            or panel["uncertainty_admitted"] is not False or panel["publication_ready"] is not False):
        raise ValueError("Require the retained descriptive panel")
    keys = [m["key"] for m in panel["methods"]]
    if len(keys) != 8 or set(keys) != {BASELINE, *METHODS}:
        raise ValueError("Require all eight distinct methods")
    summaries = {}
    for endpoint in panel["endpoints"]:
        key = endpoint["metric"], endpoint["method"]
        if key in summaries or key[0] not in METRICS or key[1] not in keys:
            raise ValueError("Duplicate or unknown endpoint")
        if endpoint["historical_execution_output_bound"] is not True:
            raise ValueError("Missing historical raw-output binding")
        summaries[key] = endpoint
    if len(summaries) != 16:
        raise ValueError("Require all 16 endpoint summaries")
    comparisons, selected = {}, {}
    for row in panel["comparisons"]:
        metric, left, right = row["metric"], row["left"], row["right"]
        key = metric, frozenset((left, right))
        if metric not in METRICS or left not in keys or right not in keys or left == right or key in comparisons:
            raise ValueError("Duplicate or invalid comparison")
        nl, nr, ns, sl, sr, ssl, ssr = counts(row["result"])
        for method, n, score_sum in ((left, nl, sl), (right, nr, sr)):
            endpoint = summaries[metric, method]
            if type(endpoint["assessed_pairs"]) is not int or endpoint["assessed_pairs"] != n:
                raise ValueError("Count differs from endpoint summary")
            close(endpoint["rounded_mean"], Fraction(score_sum, n * 1000000))
            admitted = endpoint["admitted_mean"]
            if (type(admitted) not in (int, float) or not math.isfinite(admitted)
                    or not 0 <= admitted <= 1 or abs(admitted - endpoint["rounded_mean"]) > 5.051e-7):
                raise ValueError("Serialized mean differs from admitted score")
        comparisons[key] = row
        if BASELINE not in (left, right):
            continue
        if left == BASELINE:
            method, nl, nr, sl, sr, ssl, ssr = right, nr, nl, sr, sl, ssr, ssl
        else:
            method = left
        ml, mr = Fraction(sl, nl * 1000000), Fraction(sr, nr * 1000000)
        selected[metric, method] = dict(metric=metric, method=method, reference=BASELINE,
            method_pairs=nl, reference_pairs=nr, shared_pairs=ns,
            method_shared_fraction=ns / nl, reference_shared_fraction=ns / nr,
            method_mean=float(ml), reference_mean=float(mr), difference_points=float(100 * (ml - mr)),
            shared_denominator_points=float(100 * (Fraction(ssl, nl * 1000000) - Fraction(ssr, nr * 1000000))),
            method_only_points=float(Fraction(100 * (sl - ssl), nl * 1000000)),
            negative_reference_only_points=float(-Fraction(100 * (sr - ssr), nr * 1000000)),
            shared_mean=ssl / (ns * 1000000) if ns else None,
            method_only_mean=(sl - ssl) / ((nl - ns) * 1000000) if nl > ns else None,
            reference_only_mean=(sr - ssr) / ((nr - ns) * 1000000) if nr > ns else None)
    expected = {(metric, frozenset(pair)) for metric in METRICS for pair in combinations(keys, 2)}
    if set(comparisons) != expected:
        raise ValueError("Incomplete 56-comparison panel")
    return [selected[metric, method] for metric in METRICS for method in METHODS]


def run(path, digest, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    source = record(path)
    if source["sha256"] != digest:
        raise ValueError("Panel checksum differs")
    panel = json.loads(path.read_text())
    rows = endpoints(panel)
    refs = [source, record(__file__), record(Path(__file__).with_name("prepare_ob_candidate_neighborhood.py")),
            *panel["checked_records"]]
    for pin in refs:
        check(pin)
    output.mkdir(parents=True, exist_ok=False)
    with (output / "endpoints.tsv").open("x") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    table = ["# QfO Native Scored-Pair Decomposition", "",
             "All contrasts are method minus full OrthoFinder 3.1.5. Values use serialized six-decimal pair scores.",
             "Terms sum to the original rounded-mean difference; shared-pair fractions use each method's own eligible-pair count.",
             "No new inference, scoring, confidence interval or causal attribution. GO/EC similarities are not F1.", "",
             "| Metric | Method | Scored Pairs | Shared (% Method / Reference) | Difference (Points) | Shared Denominators | Method-Only | Negative Reference-Only |",
             "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for row in rows:
        table.append(f"| {row['metric']} | {METHODS[row['method']]} | {row['method_pairs']:,} | "
                     f"{100 * row['method_shared_fraction']:.2f} / {100 * row['reference_shared_fraction']:.2f} | "
                     + " | ".join(f"{row[k]:+.6f}" for k in ("difference_points", "shared_denominator_points",
                                                             "method_only_points", "negative_reference_only_points")) + " |")
    (output / "scores.md").write_text("\n".join(table) + "\n")
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 10, "svg.hashsalt": "qfo-pair-decomposition"}):
        fig, axes = plt.subplots(1, 2, figsize=(16, 7), sharey=True)
        for axis, metric in zip(axes, METRICS):
            axis.axvline(0, color="#777777", linewidth=.8)
            for y, row in enumerate(r for r in rows if r["metric"] == metric):
                positive, negative = 0., 0.
                for key, color in zip(("shared_denominator_points", "method_only_points", "negative_reference_only_points"),
                                      ("#217f87", "#cb9624", "#bb4a57")):
                    value = row[key]
                    start = positive if value >= 0 else negative
                    axis.barh(y, value, left=start, height=.58, color=color)
                    if value >= 0:
                        positive += value
                    else:
                        negative += value
                axis.plot(row["difference_points"], y, "D", color="#222222", markersize=5, zorder=5)
            axis.set(title=metric, xlim=(-105, 105), xlabel="Additive score terms (100 x similarity difference)")
            axis.set_xticks([-100, -50, 0, 50, 100])
            axis.spines[["top", "right", "left"]].set_visible(False)
            axis.tick_params(axis="y", length=0)
        axes[0].set_yticks(range(7), METHODS.values())
        axes[0].set_ylim(6.6, -.6)
        handles = [plt.Rectangle((0, 0), 1, 1, color=color) for color in ("#217f87", "#cb9624", "#bb4a57")]
        handles.append(plt.Line2D([0], [0], marker="D", color="#222222", linestyle="none", markersize=5))
        fig.legend(handles, ["Shared-pair denominator term", "Method-only term", "Negative reference-only term", "Net difference"],
                   loc="upper center", bbox_to_anchor=(.6, .93), ncol=2, frameon=False)
        fig.suptitle("QfO: eligible-pair composition of GO/EC score differences", x=.025, ha="left", fontsize=15)
        fig.text(.025, .055, "Seven methods minus full OrthoFinder 3.1.5; shared pairs have identical serialized scores. Not confidence intervals or causal effects.", fontsize=9)
        fig.text(.025, .02, "Different prediction semantics and eligibility remain. Intersection conditioning changes the endpoint; native scores are not replaced.", fontsize=9)
        fig.subplots_adjust(left=.22, right=.985, top=.77, bottom=.16, wspace=.13)
        for suffix in ("png", "pdf", "svg"):
            fig.savefig(output / f"qfo_pair_decomposition.{suffix}", dpi=180)
        plt.close(fig)
    for pin in refs:
        check(pin)
    manifest = dict(status="retained_qfo_scored_pair_decomposition_plotted", inputs=refs,
        contrasts=len(rows), checked_comparisons=56, outputs=[record(p) for p in sorted(output.iterdir())],
        matplotlib_version=matplotlib.__version__, new_scoring=False, new_uncertainty=False,
        native_endpoints_replaced=False, publication_ready=False, visual_review_complete=False)
    with (output / "manifest.json").open("x") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.results, args.sha256, args.output)
