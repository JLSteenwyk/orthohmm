"""Plot every prespecified parameter contrast with verified paired uncertainty."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from benchmark_tools.bootstrap_orthobench import METRICS
from benchmark_tools.parameter_neighborhood_statistics import BASELINE, VARIANTS, summarize
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.run_simulation_methods import read_frozen

LABELS = {"control": "Frozen control", "cpm_low": "CPM = 0.08", "cpm_high": "CPM = 0.12",
          "norm_low": "Min. norm = 0.024", "norm_high": "Min. norm = 0.036",
          "margin_low": "Min. margin = 1.2", "margin_high": "Min. margin = 1.8"}


def validate(report):
    if (set(report["scores"]) != {BASELINE, *VARIANTS} or report["failed_variants"] != {}
            or len(report["families"]) != 70 or report["publication_ready"] is not False):
        raise ValueError("Require complete development-exposed parameter panel")
    fresh = summarize({label: score["refog_records"] for label, score in report["scores"].items()}, {})
    for key in ("baseline", "families", "replicates", "seed", "alpha", "point_estimates_percent",
                "comparisons", "multiplicity", "planned_variants", "planned_endpoints", "missing_intervals"):
        if report[key] != fresh[key]:
            raise ValueError("Displayed statistics differ from fresh paired bootstrap: " + key)
    for label in (BASELINE, *VARIANTS):
        coverage = report["coverage_resources"][label]
        if coverage["input_genes"] != 251378 or coverage["assigned_genes"] != 251378:
            raise ValueError("Incomplete gene coverage")


def plot(report):
    validate(report)
    fig = plt.figure(figsize=(12, 9))
    grid = fig.add_gridspec(2, 3, left=.21, right=.97, top=.83, bottom=.19, hspace=.5, wspace=.32,
                          height_ratios=[1, 1.2])
    table_ax = fig.add_subplot(grid[0, :])
    table_ax.axis("off")
    rows = [[LABELS[label], *[f"{report['point_estimates_percent'][label][metric]:.3f}" for metric in METRICS]]
            for label in (BASELINE, *VARIANTS)]
    table = table_ax.table(cellText=rows, colLabels=["Configuration", "F1 (%)", "Precision (%)", "Recall (%)"],
                          colWidths=[.43, .19, .19, .19], cellLoc="center", bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#dddddd")
        cell.set_linewidth(.5)
        cell.set_facecolor("#eeeeee" if row == 0 else "white")
        if row == 0:
            cell.set_text_props(weight="bold")
    extent = max(abs(value) for label in VARIANTS for metric in METRICS
                 for value in report["comparisons"][label]["metrics"][metric]["bonferroni_percentile_ci"])
    limit = max(1., float(np.ceil(extent * 1.1)))
    for i, metric in enumerate(METRICS):
        ax = fig.add_subplot(grid[1, i])
        ax.axvline(0, color="#888888", linewidth=1)
        for y, label in enumerate(VARIANTS):
            color = ("#007d83", "#b64b42", "#6b5b95")[y // 2]
            row = report["comparisons"][label]["metrics"][metric]
            ax.plot(row["bonferroni_percentile_ci"], [y, y], color=color, linewidth=1.3)
            ax.plot(row["paired_percentile_ci"], [y, y], color=color, linewidth=4, solid_capstyle="butt")
            ax.plot(row["difference_percentage_points"], y, "o", color=color, markersize=4)
        ax.set_xlim(-limit, limit)
        ax.set_ylim(5.6, -.6)
        ax.set_yticks(range(6), [LABELS[label] for label in VARIANTS] if i == 0 else [""] * 6)
        ax.set_title(("F1", "Precision", "Recall")[i], loc="left", fontsize=12)
        ax.set_xlabel("Variant minus control (pp)", fontsize=10)
        ax.grid(axis="x", alpha=.2)
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
    fig.suptitle("OrthoBench: limited parameter sensitivity", x=.04, ha="left", y=.97, fontsize=17)
    fig.text(.04, .92, "Six one-at-a-time changes of 20% | full inferred-tree pipeline | all 251,378 proteins retained", fontsize=11)
    fig.text(.04, .88, "Frozen HMM-centered method; 70 reference families. No best-variant selection or default change.", fontsize=10)
    fig.text(.04, .12, "Thick: paired 95% CI. Thin: Bonferroni CI over 18 planned endpoints. 20,000 RefOG draws.", fontsize=10)
    fig.text(.04, .07, "Development-exposed sensitivity analysis, not independent confirmation or an OrthoFinder superiority test.", fontsize=10)
    fig.text(.04, .03, "Intervals including zero do not establish equivalence or general robustness.", fontsize=10)
    return fig


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = read_frozen(args.results, args.sha256)
    figure = plot(report)
    args.output.mkdir(parents=True)
    outputs = []
    for extension in ("png", "pdf", "svg"):
        path = args.output / ("parameter_neighborhood." + extension)
        figure.savefig(path, dpi=180)
        outputs.append(record(path))
    plt.close(figure)
    (args.output / "manifest.json").write_text(json.dumps({"results": record(args.results), "source": record(__file__),
        "outputs": outputs, "endpoints": 18, "matplotlib": matplotlib.__version__}, indent=2, sort_keys=True) + "\n")
