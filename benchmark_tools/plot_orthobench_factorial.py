"""Publication figure for the frozen eight-cell OrthoBench factorial."""

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
from benchmark_tools.bootstrap_orthobench_factorial import CELLS, conditional_contrasts
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.run_simulation_methods import read_frozen

COLORS = ("#007d83", "#a66b0b", "#755297")


def validate(report):
    if (report.get("replicates") != 20000 or report.get("multiplicity_endpoints") != 36
            or report.get("seed") != 20260918 or report.get("failed_cells") != {}
            or report.get("alpha") != .05 or len(report.get("families", [])) != 70
            or set(report["point_estimates_percent"]) != set(CELLS)):
        raise ValueError("Expected complete frozen factorial and prespecified uncertainty")
    contrasts = conditional_contrasts()
    if len(report["comparisons"]) != len(contrasts):
        raise ValueError("Incomplete contrast inventory")
    for row, expected in zip(report["comparisons"], contrasts):
        if row.get("status") != "complete" or any(row.get(key) != value for key, value in expected.items()):
            raise ValueError("Changed contrast identity or order")
        for metric in METRICS:
            values = row["metrics"][metric]
            point = values["difference_percentage_points"]
            nominal, adjusted = (values[key] for key in ("paired_percentile_ci", "bonferroni_percentile_ci"))
            if (len(nominal) != 2 or len(adjusted) != 2
                    or not np.isfinite([point, *nominal, *adjusted]).all()
                    or not adjusted[0] <= nominal[0] <= nominal[1] <= adjusted[1]):
                raise ValueError("Invalid or nonnested intervals")
            difference = (report["point_estimates_percent"][row["on"]][metric]
                          - report["point_estimates_percent"][row["off"]][metric])
            if not np.isclose(point, difference, atol=1e-10, rtol=0):
                raise ValueError("Contrast disagrees with cell estimates")
    scores = np.array([[report["point_estimates_percent"][cell][key] for key in METRICS] for cell in CELLS])
    if not np.isfinite(scores).all() or np.any((scores < 0) | (scores > 100)):
        raise ValueError("Invalid percentage score")
    return scores


def plot(report):
    scores = validate(report)
    fig, axes = plt.subplots(1, 4, figsize=(19, 8), gridspec_kw={"width_ratios": [1.1, 1, 1, 1]})
    fig.subplots_adjust(left=.07, right=.985, top=.79, bottom=.25, wspace=.7)
    fig.suptitle("OrthoBench: HMM refinement, candidate expansion and reconciliation", x=.025, ha="left", y=.97, fontsize=17)
    fig.text(.025, .91, "70 reference families | 251,378 input genes | development-exposed component analysis", fontsize=11)
    fig.text(.025, .86, "P = profile refinement; C = candidate expansion; R = reconciliation. Cell labels show P C R (0: off, 1: on).", fontsize=10)
    left = axes[0]
    left.imshow(scores, vmin=0, vmax=100, cmap="Greys", aspect="auto")
    left.set_xticks(range(3), ["F1", "Precision", "Recall"])
    left.set_yticks(range(8), [cell.replace("p", "").replace("c", "").replace("r", "").replace("_", " ") for cell in CELLS])
    for i in range(8):
        for j in range(3):
            left.text(j, i, f"{scores[i, j]:.2f}", color="white" if scores[i, j] >= 55 else "black", ha="center", va="center", fontsize=10)
    left.set_title("A  Observed scores (%)", loc="left", pad=12)
    for k, ax in enumerate(axes[1:]):
        ax.axvline(0, color="#888888", linewidth=1)
        labels = []
        for i, row in enumerate(report["comparisons"]):
            value = row["metrics"][METRICS[k]]
            color = COLORS[i // 4]
            ax.plot(value["bonferroni_percentile_ci"], [i, i], color=color, linewidth=1)
            ax.plot(value["paired_percentile_ci"], [i, i], color=color, linewidth=4, solid_capstyle="butt")
            ax.plot(value["difference_percentage_points"], i, "o", color=color, markersize=4)
            fixed = row["fixed"]
            short = {"profile_expansion": "P", "candidate_expansion": "C", "reconciliation": "R"}
            labels.append(short[row["factor"]] + " | " + " ".join(f"{short[f]}{v}" for f, v in fixed.items()))
        ax.set_yticks(range(12), labels if k == 0 else [""] * 12, fontsize=9)
        ax.set_ylim(11.6, -.6)
        for boundary in (3.5, 7.5):
            ax.axhline(boundary, color="#dddddd", linewidth=1)
        ax.set_xlabel("On minus off (percentage points)", fontsize=9)
        ax.set_title(f"{'BCD'[k]}  {('F1', 'Precision', 'Recall')[k]} effects", loc="left", pad=12)
        ax.grid(axis="x", alpha=.2)
    for ax in axes:
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
    fig.text(.025, .16, "Rows: profile effects (teal), candidate effects (ochre), reconciliation effects (purple); remaining factors held fixed.", fontsize=10)
    fig.text(.025, .115, "Thick: nominal 95% CI. Thin: Bonferroni CI across all 36 endpoints. 20,000 paired RefOG draws; seed 20260918.", fontsize=10)
    fig.text(.025, .07, "Profile-off retains initial HMM search. Reconciliation changes candidate co-membership to root HOGs. No independent superiority claim.", fontsize=10)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    report = read_frozen(args.results, args.sha256)
    if args.output.exists():
        raise FileExistsError(args.output)
    fig = plot(report)
    args.output.mkdir(parents=True)
    outputs = []
    for extension in ("png", "pdf", "svg"):
        path = args.output / ("orthobench_factorial." + extension)
        fig.savefig(path, dpi=180)
        outputs.append(file_provenance(path))
    plt.close(fig)
    manifest = {"source_results": file_provenance(args.results), "plotter": file_provenance(Path(__file__)),
                "outputs": outputs, "matplotlib": matplotlib.__version__, "endpoints_shown": 36}
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
