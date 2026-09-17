"""Plot the admitted initial-search replacement control and all paired effects."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from benchmark_tools.bootstrap_orthobench import METRICS, statistics, weighted_records
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.run_simulation_methods import read_frozen

BASELINE = "p0_c0_r0"
VARIANTS = ("all_hits", "top100")
LABELS = {BASELINE: "HMM search", "all_hits": "DIAMOND: all hits", "top100": "DIAMOND: top 100"}
RESULT_SHA = "b374010483f56d1ceacec956fbd4ad67823e53d81f4219f81704df8409ffd829"


def validate(report):
    if (report["baseline"] != BASELINE or report["replicates"] != 20000 or report["seed"] != 20260918
            or report["alpha"] != .05 or len(report["families"]) != 70 or len(set(report["families"])) != 70
            or report["multiplicity"] != "Bonferroni tail adjustment over 6 reported contrasts/metrics"
            or set(report["point_estimates_percent"]) != {BASELINE, *VARIANTS}
            or set(report["comparisons"]) != set(VARIANTS)):
        raise ValueError("Changed search-control panel or uncertainty specification")
    for label in (BASELINE, *VARIANTS):
        names, _, counts = weighted_records(report["scores"][label]["refog_records"])
        point = np.array([report["point_estimates_percent"][label][metric] for metric in METRICS])
        if names != sorted(report["families"]) or not np.allclose(point, statistics(counts.sum(axis=0)), atol=1e-10, rtol=0):
            raise ValueError("Displayed scores differ from sufficient statistics")
        coverage = report["gene_coverage"][label]
        if coverage["input_genes"] != 251378 or coverage["assigned_genes"] != 251378:
            raise ValueError("Incomplete gene coverage")
    for label in VARIANTS:
        comparison = report["comparisons"][label]
        if comparison["versus"] != BASELINE or sum(comparison[key] for key in
                ("family_f1_wins", "family_f1_ties", "family_f1_losses")) != 70:
            raise ValueError("Changed paired comparison inventory")
        for metric in METRICS:
            row = comparison["metrics"][metric]
            point = row["difference_percentage_points"]
            nominal, adjusted = row["paired_percentile_ci"], row["bonferroni_percentile_ci"]
            if (len(nominal) != 2 or len(adjusted) != 2 or not np.isfinite([point, *nominal, *adjusted]).all()
                    or not adjusted[0] <= nominal[0] <= nominal[1] <= adjusted[1]):
                raise ValueError("Invalid uncertainty interval")
            expected = report["point_estimates_percent"][label][metric] - report["point_estimates_percent"][BASELINE][metric]
            if not np.isclose(point, expected, atol=1e-10, rtol=0):
                raise ValueError("Displayed effect differs from scores")


def plot(report):
    validate(report)
    fig = plt.figure(figsize=(12, 7))
    grid = fig.add_gridspec(2, 3, left=.19, right=.97, top=.78, bottom=.23, hspace=.72, wspace=.32)
    table_ax = fig.add_subplot(grid[0, :])
    table_ax.axis("off")
    rows = [[LABELS[label], *[f"{report['point_estimates_percent'][label][metric]:.2f}" for metric in METRICS]]
            for label in (BASELINE, *VARIANTS)]
    table = table_ax.table(cellText=rows, colLabels=["Initial search", "F1 (%)", "Precision (%)", "Recall (%)"],
                           colWidths=[.43, .19, .19, .19], cellLoc="center", bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#dddddd")
        cell.set_linewidth(.5)
        cell.set_facecolor("#eeeeee" if row == 0 else "white")
        if row == 0:
            cell.set_text_props(weight="bold")
    axes = [fig.add_subplot(grid[1, i]) for i in range(3)]
    for i, ax in enumerate(axes):
        ax.axvline(0, color="#888888", linewidth=1)
        for y, (label, color) in enumerate(zip(VARIANTS, ("#007d83", "#b64b42"))):
            value = report["comparisons"][label]["metrics"][METRICS[i]]
            ax.plot(value["bonferroni_percentile_ci"], [y, y], color=color, linewidth=1.4)
            ax.plot(value["paired_percentile_ci"], [y, y], color=color, linewidth=4, solid_capstyle="butt")
            ax.plot(value["difference_percentage_points"], y, "o", color=color, markersize=5)
        ax.set_xlim(-40, 40)
        ax.set_xticks([-40, -20, 0, 20, 40])
        ax.set_ylim(1.65, -.65)
        ax.set_yticks([0, 1], [LABELS[label] for label in VARIANTS] if i == 0 else ["", ""])
        ax.set_title(("F1", "Precision", "Recall")[i], loc="left", fontsize=12)
        ax.set_xlabel("DIAMOND minus HMM (pp)", fontsize=10)
        ax.grid(axis="x", alpha=.2)
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
    fig.suptitle("OrthoBench: initial-search replacement control", x=.04, ha="left", y=.97, fontsize=17)
    fig.text(.04, .90, "70 reference families | all 251,378 proteins retained | fixed downstream grouping and refinement", fontsize=11)
    fig.text(.04, .85, "Profile expansion, candidate expansion and phylogenetic reconciliation are off in every arm.", fontsize=10)
    fig.text(.04, .14, "Thick: paired 95% CI. Thin: Bonferroni CI over six endpoints. 20,000 RefOG draws; seed 20260918.", fontsize=10)
    fig.text(.04, .09, "DIAMOND increases recall and reduces precision. Both adjusted F1 intervals include zero.", fontsize=10)
    fig.text(.04, .04, "Development-exposed control; search sensitivity and score calibration are not matched. No established HMM F1 advantage.", fontsize=10)
    return fig


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = read_frozen(args.results, RESULT_SHA)
    figure = plot(report)
    args.output.mkdir(parents=True)
    outputs = []
    for extension in ("png", "pdf", "svg"):
        path = args.output / ("sequence_search_control." + extension)
        figure.savefig(path, dpi=180)
        outputs.append(record(path))
    plt.close(figure)
    (args.output / "manifest.json").write_text(json.dumps({"results": record(args.results), "source": record(__file__),
        "outputs": outputs, "endpoints": 6, "matplotlib": matplotlib.__version__}, indent=2, sort_keys=True) + "\n")
