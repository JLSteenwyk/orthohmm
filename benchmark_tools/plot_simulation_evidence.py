"""Plot corrected simulation admission counts and complete-case F1 contrasts."""

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
from benchmark_tools.summarize_simulation_panel import CONDITIONS, METHODS, PANELS
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.orthobench_stage_diagnostics import file_provenance

NAMES = ("Baseline", "Divergence", "Turnover", "Divergence + turnover", "20% missing", "Uneven taxa", "Taxon-count control")
METHOD_NAMES = ("OrthoHMM\nhigh sensitivity", "OrthoHMM\nsatellite_v2", "OrthoFinder\nfull", "OrthoFinder\ncheckpoint")
COLORS = ("#007f86", "#b64062")


def plot(report):
    panel = report["panel"]
    if panel not in PANELS:
        raise ValueError("Unknown panel")
    counts = np.array([[len(report["conditions"][c]["methods"][m]["complete_seeds"]) for m in METHODS] for c in CONDITIONS])
    fig, (left, right) = plt.subplots(1, 2, figsize=(14, 7), gridspec_kw={"width_ratios": [1, 1.5]})
    fig.subplots_adjust(left=.15, right=.98, bottom=.19, top=.82, wspace=.38)
    left.imshow(counts, cmap="Greys", vmin=0, vmax=12, aspect="auto")
    left.set_xticks(range(4), METHOD_NAMES, fontsize=9)
    left.set_yticks(range(7), NAMES, fontsize=10)
    left.set_title("A  Admitted datasets (out of 10)", loc="left", fontsize=12, pad=14)
    for i in range(7):
        for j in range(4):
            left.text(j, i, str(counts[i, j]), ha="center", va="center", color="white" if counts[i, j] >= 7 else "black", fontsize=12)
    right.set_title("B  Paired F1 difference versus full OrthoFinder", loc="left", fontsize=12, pad=14)
    right.axvline(0, color="#777777", linewidth=1)
    labels, positions = [], []
    estimated = 0
    for i, condition in enumerate(CONDITIONS):
        for j, method in enumerate(METHODS[:2]):
            y = 2 * i + j
            row = report["conditions"][condition]["contrasts"][method]
            n = len(row["included_seeds"])
            labels.append(f"n={n}")
            positions.append(y)
            metric = row.get("metrics", {}).get("f1")
            if metric:
                estimated += 1
                point = metric["difference_percentage_points"]
                adjusted = metric["bonferroni_14_ci"]
                nominal = metric["paired_95_percent_ci"]
                right.plot(adjusted, [y, y], color=COLORS[j], linewidth=1)
                right.plot(nominal, [y, y], color=COLORS[j], linewidth=4, solid_capstyle="butt")
                right.plot(point, y, "o", color=COLORS[j], markersize=5)
    right.set_yticks(positions, labels, fontsize=9)
    right.set_ylim(13.7, -.7)
    right.set_xlabel("OrthoHMM minus OrthoFinder (percentage points)", fontsize=10)
    right.grid(axis="x", alpha=.2)
    for ax in (left, right):
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.tick_params(length=0)
    if estimated:
        right.set_xlim(min(-1, right.get_xlim()[0]), max(1, right.get_xlim()[1]))
    else:
        right.set_xlim(-1, 1)
        right.text(.5, .5, "No admitted paired comparisons\nAll comparator runs failed native admission",
                   transform=right.transAxes, ha="center", va="center", bbox={"facecolor": "white", "edgecolor": "none"})
        right.set_xticks([])
    title = "Variable family lengths" if panel == "variable_length_v2" else "Fixed-length stress panel"
    fig.suptitle(f"Corrected simulations: {title}", fontsize=17, x=.15, ha="left", y=.95)
    fig.legend([Line2D([], [], color=c, marker="o") for c in COLORS],
               ["OrthoHMM high sensitivity", "OrthoHMM satellite_v2"], loc="upper left", bbox_to_anchor=(.145, .92), ncol=2, frameon=False)
    fig.text(.15, .09, "Thick lines: nominal 95% CI. Thin lines: Bonferroni-14 CI. Bootstrap unit: independent simulation seed.", fontsize=10)
    fig.text(.15, .055, "Comparisons condition on joint success; failed runs receive no imputed score. Panels are not pooled.", fontsize=10)
    fig.text(.15, .02, "Ten seeds limit interval precision. The sequence checkpoint requires a valid full parent run.", fontsize=9, color="#555555")
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = read_frozen(args.results, args.sha256)
    if args.output.exists():
        raise FileExistsError(args.output)
    fig = plot(report)
    args.output.mkdir(parents=True)
    outputs = []
    for extension in ("png", "pdf", "svg"):
        path = args.output / f"simulation_evidence.{extension}"
        fig.savefig(path, dpi=180)
        outputs.append(file_provenance(path))
    plt.close(fig)
    manifest = {"source_results": file_provenance(args.results), "plotter": file_provenance(Path(__file__)),
                "matplotlib": matplotlib.__version__, "outputs": outputs, "panel": report["panel"]}
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
