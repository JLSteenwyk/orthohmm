"""Render all prespecified simulation tree contrasts from reproducible summaries."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.summarize_simulation_tree_panel import CONDITIONS, METHODS, CONTRASTS, METRICS, summarize

COLORS = ("#007d83", "#bd641c")
LABELS = ("Baseline", "Divergent", "Turnover", "Divergent + turnover", "Missing 20%",
          "Uneven taxa", "Taxon-count control")


def validate(report):
    expected = summarize(report["records"])
    for key in ("status", "contrasts", "bootstrap", "statistic", "inference_role"):
        if report[key] != expected[key]:
            raise ValueError("Figure evidence differs from the fixed complete-panel analysis")
    return {(r["condition"], r["method"], r["target"], r["reference"]): r for r in report["contrasts"]}


def plot(report):
    indexed = validate(report)
    fig, axes = plt.subplots(3, 3, figsize=(16, 13))
    fig.subplots_adjust(left=.16, right=.98, top=.84, bottom=.17, hspace=.37, wspace=.25)
    fig.suptitle("Simulation sensitivity to supplied species trees", x=.025, ha="left", y=.98, fontsize=18)
    failed = sum(r["status"] == "failed" for r in report["records"])
    fig.text(.025, .945, f"70 datasets | 560 planned method/tree outcomes | {failed} failures retained | 126 exploratory endpoints", fontsize=11)
    fig.legend([Line2D([], [], color=c, marker=m, linestyle="") for c, m in zip(COLORS, ("o", "s"))],
               ["OrthoHMM satellite_v2", "OrthoFinder full 3.1.5"], loc="upper left",
               bbox_to_anchor=(.02, .928), ncol=2, frameon=False, fontsize=11)
    titles = ("Generating minus inferred", "NNI1 minus generating", "NNI2 minus generating")
    for col, (target, reference) in enumerate(CONTRASTS):
        bounds = [v for row in report["contrasts"] if (row["target"], row["reference"]) == (target, reference)
                  for value in row["metrics"].values() for v in (value["bonferroni_126_ci"] or [])]
        bounds.extend(value["difference_percentage_points"] for row in report["contrasts"]
                      if (row["target"], row["reference"]) == (target, reference)
                      for value in row["metrics"].values() if value["difference_percentage_points"] is not None)
        lo, hi = min([0., *bounds]), max([0., *bounds])
        pad = max((hi - lo) * .12, .1)
        for k, metric in enumerate(METRICS):
            ax = axes[k, col]
            ax.axvline(0, color="#999999", linewidth=.8)
            for i, condition in enumerate(CONDITIONS):
                for j, method in enumerate(METHODS):
                    value = indexed[condition, method, target, reference]["metrics"][metric]
                    y = i + (-.14 if j == 0 else .14)
                    point = value["difference_percentage_points"]
                    if point is None:
                        ax.text(.98, y, "unavailable", transform=ax.get_yaxis_transform(),
                                ha="right", va="center", fontsize=6, color=COLORS[j])
                        continue
                    if value["bonferroni_126_ci"] is not None:
                        ax.plot(value["bonferroni_126_ci"], [y, y], color=COLORS[j], linewidth=1)
                        ax.plot(value["paired_95_percent_ci"], [y, y], color=COLORS[j], linewidth=3)
                    ax.plot(point, y, "o" if j == 0 else "s", color=COLORS[j], markersize=4)
            ax.set_xlim(lo - pad, hi + pad)
            ax.set_ylim(6.6, -.6)
            ax.set_yticks(range(7), LABELS if col == 0 else [""] * 7, fontsize=9)
            ax.set_title(f"{'ABCDEFGHI'[k * 3 + col]}  {titles[col]}", loc="left", fontsize=10, pad=10)
            ax.set_xlabel(("F1", "Precision", "Recall")[k] + " difference (percentage points)", fontsize=9)
            ax.grid(axis="x", alpha=.15)
            ax.tick_params(length=0, labelsize=9)
            for spine in ax.spines.values():
                spine.set_visible(False)
    counts = sorted({c["paired_seed_count"] for c in report["contrasts"]})
    fig.text(.025, .115, "Dots/squares: paired-seed means. Thick: nominal 95% CI. Thin: Bonferroni CI across all 126 endpoints.", fontsize=10)
    fig.text(.025, .085, "20,000 paired-seed draws; seed 20260918. Available paired counts: " + ", ".join(map(str, counts)) + ". Failures are not assigned zero accuracy.", fontsize=10)
    fig.text(.025, .055, "Generating tree is an oracle diagnostic. Retained upstream differences and per-contrast exclusions are reported in the accompanying audit.", fontsize=9)
    fig.text(.025, .027, "Development-exposed simplified simulations; small seed counts and approximate bootstrap tails limit generalization and simultaneous coverage.", fontsize=9)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = read_frozen(args.results, args.sha256)
    fig = plot(report)
    args.output.mkdir(parents=True)
    outputs = []
    for extension in ("png", "pdf", "svg"):
        path = args.output / ("simulation_tree_robustness." + extension)
        fig.savefig(path, dpi=180)
        outputs.append(record(path))
    plt.close(fig)
    manifest = {"source_results": record(args.results), "plotter": record(__file__),
                "summary_source": record(Path(__file__).with_name("summarize_simulation_tree_panel.py")),
                "outputs": outputs, "matplotlib": matplotlib.__version__, "endpoints_shown": 126}
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
