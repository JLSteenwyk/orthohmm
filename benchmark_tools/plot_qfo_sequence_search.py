"""Render every prespecified corrected SwissTrees sequence-control contrast."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen

RESULT_SHA = "e4e4696e6f090d7d83b99e5190bd21606a980fac294c62cdf260582d822d89f3"
VARIANTS = ("p0_c0_r0", "all_hits", "top100")
LABELS = ("Initial HMM", "DIAMOND: all hits", "DIAMOND: top 100")
METRICS = ("F1", "PPV", "TPR")


def validate(report):
    if (report.get("status") != "paired_corrected_sequence_swiss_intervals"
            or report.get("uncertainty_admitted") is not True
            or report.get("publication_ready") is not False
            or report["replicates"] != 100000 or report["seed"] != 20260923
            or report["multiplicity_endpoints"] != 6 or report["alpha"] != .05
            or report["quantile_method"] != "linear"
            or report["units"] != "raw 0-to-1 metric units"
            or len(report["families"]) != 18 or len(set(report["families"])) != 18
            or set(report["point_estimates"]) != set(VARIANTS)
            or [c["candidate"] for c in report["comparisons"]] != list(VARIANTS[1:])):
        raise ValueError("Changed admitted panel or uncertainty protocol")
    points = report["point_estimates"]
    for variant in VARIANTS:
        row = points[variant]
        if set(row) != set(METRICS) or any(type(v) not in (int, float)
                or not np.isfinite(v) or not 0 <= v <= 1 for v in row.values()):
            raise ValueError("Invalid score")
        p, r = row["PPV"], row["TPR"]
        if not np.isclose(row["F1"], 2 * p * r / (p + r) if p + r else 0, atol=1e-12, rtol=0):
            raise ValueError("F1 differs from macro precision and recall")
    for contrast in report["comparisons"]:
        if contrast["reference"] != VARIANTS[0] or set(contrast["metrics"]) != set(METRICS):
            raise ValueError("Changed contrast orientation or endpoints")
        for metric, row in contrast["metrics"].items():
            nominal, adjusted = row["paired_percentile_ci"], row["bonferroni_percentile_ci"]
            if (len(nominal) != 2 or len(adjusted) != 2
                    or not np.isfinite([row["difference"], *nominal, *adjusted]).all()
                    or not -1 <= adjusted[0] <= nominal[0] <= nominal[1] <= adjusted[1] <= 1):
                raise ValueError("Invalid interval")
            expected = points[contrast["candidate"]][metric] - points[VARIANTS[0]][metric]
            if not np.isclose(row["difference"], expected, atol=1e-12, rtol=0):
                raise ValueError("Effect disagrees with point estimates")
            counts = [row[k] for k in ("family_wins", "family_ties", "family_losses")]
            if any(type(v) is not int or v < 0 for v in counts) or sum(counts) != 18:
                raise ValueError("Invalid family inventory")
            if any(abs(v) > .30 for v in (*nominal, *adjusted, row["difference"])):
                raise ValueError("Fixed plotting range would clip evidence")


def plot(report):
    validate(report)
    fig = plt.figure(figsize=(12, 7))
    grid = fig.add_gridspec(2, 3, left=.19, right=.97, top=.78, bottom=.23, hspace=.72, wspace=.32)
    table_ax = fig.add_subplot(grid[0, :])
    table_ax.axis("off")
    rows = [[label, *[f"{100 * report['point_estimates'][variant][m]:.2f}" for m in METRICS]]
            for variant, label in zip(VARIANTS, LABELS)]
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
    for i, metric in enumerate(METRICS):
        ax = fig.add_subplot(grid[1, i])
        ax.axvline(0, color="#888888", linewidth=1)
        for y, (contrast, color) in enumerate(zip(report["comparisons"], ("#007d83", "#b64b42"))):
            row = contrast["metrics"][metric]
            ax.plot(100 * np.array(row["bonferroni_percentile_ci"]), [y, y], color=color, linewidth=1.4)
            ax.plot(100 * np.array(row["paired_percentile_ci"]), [y, y], color=color, linewidth=4,
                    solid_capstyle="butt")
            ax.plot(100 * row["difference"], y, "o", color=color, markersize=5)
        ax.set_xlim(-30, 30)
        ax.set_xticks([-30, -15, 0, 15, 30])
        ax.set_ylim(1.65, -.65)
        ax.set_yticks([0, 1], LABELS[1:] if i == 0 else ["", ""])
        ax.set_title(("F1", "Precision", "Recall")[i], loc="left", fontsize=12)
        ax.set_xlabel("DIAMOND minus HMM (pp)", fontsize=10)
        ax.grid(axis="x", alpha=.2)
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
    fig.suptitle("Corrected QfO SwissTrees: initial-search replacement", x=.04, ha="left", y=.97, fontsize=17)
    fig.text(.04, .90, "18 reference families | fixed downstream grouping and refinement | native macro statistic", fontsize=11)
    fig.text(.04, .85, "Profile expansion, candidate expansion and phylogenetic reconciliation are off in every arm.", fontsize=10)
    fig.text(.04, .14, "Thick: paired 95% CI. Thin: Bonferroni CI over six endpoints. 100,000 family draws; seed 20260923.", fontsize=10)
    fig.text(.04, .09, "Negative differences favor HMM. Adjusted F1 and recall intervals include zero; precision intervals do not.", fontsize=10)
    fig.text(.04, .04, "Development-exposed control; equal E-values do not establish matched sensitivity or cost. No established HMM F1 advantage.", fontsize=10)
    return fig


def render(results, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    report = read_frozen(results, RESULT_SHA)
    source = record(results)
    figure = plot(report)
    check(source)
    output.mkdir(parents=True)
    outputs = []
    for extension in ("png", "pdf", "svg"):
        path = output / ("sequence_search_control." + extension)
        figure.savefig(path, dpi=180)
        outputs.append(record(path))
    plt.close(figure)
    check(source)
    result = {"results": source, "source": record(__file__), "outputs": outputs,
        "endpoints": 6, "axis_units": "percentage points", "score_units": "percent",
        "matplotlib": matplotlib.__version__, "publication_ready": False,
        "limitations": "Renders the pinned admitted analysis; does not rerun its upstream source audit."}
    (output / "manifest.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    render(args.results.resolve(), args.output.absolute())
