"""Plot the complete frozen OrthoBench stratified precision-recall panel."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from benchmark_tools.assemble_ob_stratified_errors import BASELINE, COMPARATORS, METHODS, validate_strata
from benchmark_tools.bootstrap_orthobench import METRICS, statistics, weighted_records
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.prepare_ob_error_strata import CATEGORIES
from benchmark_tools.run_simulation_methods import read_frozen

LABELS = ("Size: 2-20 genes", "Size: 21-50 genes", "Size: >50 genes",
          "Copy number: single", "Copy number: multiple",
          "Identity: lower", "Identity: higher", "Identity: missing",
          "Relative length: short", "Relative length: not short", "Relative length: missing",
          "Composition: concentrated", "Composition: not concentrated", "Composition: missing")
ORDER = tuple(d + ":" + c for d, cats in CATEGORIES.items() for c in cats)


def validate(report):
    if (report["status"] != "stratified_analysis_complete" or report["baseline"] != BASELINE
            or report["multiplicity_endpoints"] != 84 or report["replicates"] != 20000
            or report["seed"] != 20260918 or report["alpha"] != .05
            or set(report["strata"]) != set(ORDER) or set(report["scores"]) != set(METHODS)):
        raise ValueError("Incomplete panel or changed uncertainty specification")
    validate_strata(report["feature_admission"])
    counts = {"planned_endpoints": 84, "finite_points": 0, "interval_endpoints": 0, "nonestimable_endpoints": 0}
    for label in ORDER:
        row = report["strata"][label]
        names = report["feature_admission"]["strata"][label]
        n = len(names)
        status = "paired_bootstrap" if n >= 5 else "descriptive_only_lt5" if n else "empty_nonestimable"
        if (row["families"] != names or row["family_count"] != n or row["status"] != status
                or row["multiplicity_endpoints"] != 84 or set(row["comparisons"]) != set(COMPARATORS)
                or set(row["point_estimates_percent"]) != set(METHODS)):
            raise ValueError("Changed stratum membership or status")
        for method in METHODS:
            point = row["point_estimates_percent"][method]
            if not n:
                if point is not None:
                    raise ValueError("Empty stratum has a score")
                continue
            records = [record for record in report["scores"][method]["refog_records"] if record["refog"] in names]
            actual_names, _, weighted = weighted_records(records)
            expected = statistics(weighted.sum(axis=0))
            if (actual_names != sorted(names) or set(point) != set(METRICS)
                    or not np.allclose([point[m] for m in METRICS], expected, rtol=0, atol=1e-10)):
                raise ValueError("Point scores differ from full-reference sufficient statistics")
        if n >= 5:
            bootstrap = row["bootstrap"]
            if (any(bootstrap[key] != report[key] for key in ("baseline", "replicates", "seed", "alpha"))
                    or bootstrap["multiplicity"] != "Bonferroni tail adjustment over 84 planned endpoints"
                    or bootstrap["families"] != names or bootstrap["comparisons"] != row["comparisons"]
                    or bootstrap["point_estimates_percent"] != row["point_estimates_percent"]):
                raise ValueError("Changed bootstrap evidence")
        for method in COMPARATORS:
            contrast = row["comparisons"][method]
            if (contrast["versus"] != BASELINE or set(contrast["metrics"]) != set(METRICS)
                    or sum(contrast[key] for key in ("family_f1_wins", "family_f1_ties", "family_f1_losses")) != n):
                raise ValueError("Changed contrast inventory")
            for metric in METRICS:
                value = contrast["metrics"][metric]
                point, nominal, adjusted = (value[key] for key in
                    ("difference_percentage_points", "paired_percentile_ci", "bonferroni_percentile_ci"))
                if not n:
                    if any(item is not None for item in (point, nominal, adjusted)):
                        raise ValueError("Empty stratum has an effect or interval")
                    counts["nonestimable_endpoints"] += 1
                    continue
                expected = row["point_estimates_percent"][method][metric] - row["point_estimates_percent"][BASELINE][metric]
                if not np.isfinite(point) or not np.isclose(point, expected, atol=1e-10, rtol=0):
                    raise ValueError("Contrast differs from observed scores")
                counts["finite_points"] += 1
                if n < 5:
                    if nominal is not None or adjusted is not None:
                        raise ValueError("Small stratum must not have an interval")
                else:
                    if (len(nominal) != 2 or len(adjusted) != 2 or not np.isfinite([*nominal, *adjusted]).all()
                            or not -100 <= adjusted[0] <= nominal[0] <= nominal[1] <= adjusted[1] <= 100):
                        raise ValueError("Invalid or nonnested intervals")
                    counts["interval_endpoints"] += 1
    return counts


def plot(report):
    validate(report)
    fig, axes = plt.subplots(1, 3, figsize=(16, 11.5), sharey=True)
    fig.subplots_adjust(left=.24, right=.98, bottom=.19, top=.84, wspace=.10)
    fig.suptitle("OrthoBench: stratified accuracy tradeoffs", x=.035, ha="left", y=.97, fontsize=19)
    fig.text(.035, .925, "OrthoHMM minus full OrthoFinder 3.1.5 | all 14 prespecified strata | 70 reference families", fontsize=11)
    colors, markers = ("#007d83", "#b04f39"), ("o", "s")
    limits = [0.]
    for row in report["strata"].values():
        for comparison in row["comparisons"].values():
            for value in comparison["metrics"].values():
                if value["difference_percentage_points"] is not None:
                    limits.append(value["difference_percentage_points"])
                limits.extend(value["bonferroni_percentile_ci"] or [])
    bounds = (5 * np.floor(min(limits) / 5) - 5, 5 * np.ceil(max(limits) / 5) + 5)
    for k, (ax, metric) in enumerate(zip(axes, METRICS)):
        ax.axvline(0, color="#888888", linewidth=1)
        for i, label in enumerate(ORDER):
            row = report["strata"][label]
            if not row["family_count"]:
                ax.text(.5, i, "No families", transform=ax.get_yaxis_transform(), ha="center", va="center", color="#666666", fontsize=9)
                continue
            for j, method in enumerate(COMPARATORS):
                value = row["comparisons"][method]["metrics"][metric]
                y = i + (-.16 if j == 0 else .16)
                if row["status"] == "paired_bootstrap":
                    ax.plot(value["bonferroni_percentile_ci"], [y, y], color=colors[j], linewidth=1)
                    ax.plot(value["paired_percentile_ci"], [y, y], color=colors[j], linewidth=3.5, solid_capstyle="butt")
                ax.plot(value["difference_percentage_points"], y, marker=markers[j], linestyle="none", color=colors[j],
                        markerfacecolor="white" if row["family_count"] < 5 else colors[j], markersize=5)
        for boundary in (2.5, 4.5, 7.5, 10.5):
            ax.axhline(boundary, color="#dddddd", linewidth=.8)
        ax.set_xlim(*bounds)
        ax.set_ylim(13.6, -.6)
        ax.set_title(f"{'ABC'[k]}  {('F1', 'Precision', 'Recall')[k]}", loc="left", pad=12, fontsize=13)
        ax.set_xlabel("Difference (percentage points)", fontsize=10)
        ax.grid(axis="x", alpha=.15)
        ax.tick_params(length=0, labelsize=9)
        for spine in ax.spines.values():
            spine.set_visible(False)
    axes[0].set_yticks(range(14), [f"{label}  (n={report['strata'][key]['family_count']})" for key, label in zip(ORDER, LABELS)])
    handles = [Line2D([], [], color=color, marker=marker, linestyle="none", label=label)
               for color, marker, label in zip(colors, markers, ("OrthoHMM sensitive", "OrthoHMM phylogenetic"))]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(.03, .9), ncol=2, frameon=False, fontsize=11)
    f1_intervals = [comparison["metrics"]["f_score"]["bonferroni_percentile_ci"]
                    for row in report["strata"].values() for comparison in row["comparisons"].values()
                    if row["status"] == "paired_bootstrap"]
    excluding_zero = sum(lo > 0 or hi < 0 for lo, hi in f1_intervals)
    notes = ["Thick: nominal 95% CI. Thin: Bonferroni CI over all 84 planned endpoints. Open points: n<5, descriptive only.",
             "20,000 paired RefOG draws; seed 20260918. Adjusted tails contain about six draws each; intervals are approximate.",
             "Exploratory, development-exposed analysis. Strata overlap; descriptors do not establish duplication, fragment or domain mechanisms.",
             f"Full-reference scoring retained. Empty strata are not zero. Adjusted F1 intervals excluding zero: {excluding_zero}/{len(f1_intervals)}."]
    for y, note in zip((.13, .098, .066, .034), notes):
        fig.text(.035, y, note, fontsize=10)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    report = read_frozen(args.results, args.sha256)
    counts = validate(report)
    if args.output.exists():
        raise FileExistsError(args.output)
    fig = plot(report)
    args.output.mkdir(parents=True)
    outputs = []
    for extension in ("png", "pdf", "svg"):
        path = args.output / ("stratified_errors." + extension)
        fig.savefig(path, dpi=180)
        outputs.append(file_provenance(path))
    plt.close(fig)
    (args.output / "manifest.json").write_text(json.dumps({"source_results": file_provenance(args.results),
        "plotter": file_provenance(Path(__file__)), "outputs": outputs, "matplotlib": matplotlib.__version__,
        "endpoint_inventory": counts, "stratum_order": ORDER}, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
