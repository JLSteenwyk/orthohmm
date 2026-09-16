"""Plot every prespecified species-tree perturbation and its paired effects."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from benchmark_tools.assemble_species_tree_robustness import BASELINE, native_methods
from benchmark_tools.bootstrap_orthobench import METRICS
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.run_species_tree_control import PERTURBATIONS


def validate(report):
    native_methods(report["native_validation"])
    if (report["baseline"] != BASELINE or report["replicates"] != 20000 or report["seed"] != 20260918
            or report["alpha"] != .05 or len(set(report["families"])) != 70 or len(report["families"]) != 70
            or report["multiplicity"] != "Bonferroni tail adjustment over 18 reported contrasts/metrics"
            or set(report["point_estimates_percent"]) != {BASELINE, *PERTURBATIONS}
            or set(report["comparisons"]) != set(PERTURBATIONS)):
        raise ValueError("Incomplete panel or changed uncertainty specification")
    for label in PERTURBATIONS:
        row = report["comparisons"][label]
        if row["versus"] != BASELINE or sum(row[key] for key in ("family_f1_wins", "family_f1_ties", "family_f1_losses")) != 70:
            raise ValueError("Changed baseline or family comparison inventory")
        for metric in METRICS:
            value = row["metrics"][metric]
            point = value["difference_percentage_points"]
            nominal, adjusted = (value[key] for key in ("paired_percentile_ci", "bonferroni_percentile_ci"))
            if (len(nominal) != 2 or len(adjusted) != 2 or not np.isfinite([point, *nominal, *adjusted]).all()
                    or not adjusted[0] <= nominal[0] <= nominal[1] <= adjusted[1]):
                raise ValueError("Invalid or nonnested intervals")
            expected = report["point_estimates_percent"][label][metric] - report["point_estimates_percent"][BASELINE][metric]
            if not np.isclose(point, expected, atol=1e-10, rtol=0):
                raise ValueError("Contrast differs from observed scores")
    scores = np.array([[report["point_estimates_percent"][label][key] for key in METRICS]
                       for label in (BASELINE, *PERTURBATIONS)])
    if not np.isfinite(scores).all() or np.any((scores < 0) | (scores > 100)):
        raise ValueError("Invalid percentage score")
    return scores


def plot(report):
    scores = validate(report)
    fig, axes = plt.subplots(1, 4, figsize=(17, 6.8), gridspec_kw={"width_ratios": [1.15, 1, 1, 1]})
    fig.subplots_adjust(left=.08, right=.98, top=.76, bottom=.25, wspace=.65)
    fig.suptitle("OrthoBench: sensitivity to supplied species-tree topology", x=.025, ha="left", y=.97, fontsize=17)
    fig.text(.025, .90, "70 reference families | fixed candidates and reconciliation rules | raw gene trees reused", fontsize=11)
    fig.text(.025, .84, "Unchanged supplied control reproduces the inferred baseline. Six prespecified topology perturbations; no best-tree selection.", fontsize=10)
    labels = ["Control", *["NNI " + name[-1] for name in PERTURBATIONS[:3]],
              *["2-NNI " + name[-1] for name in PERTURBATIONS[3:]]]
    axes[0].imshow(scores, vmin=0, vmax=100, cmap="Greys", aspect="auto")
    axes[0].set_xticks(range(3), ["F1", "Precision", "Recall"])
    axes[0].set_yticks(range(7), labels)
    axes[0].set_title("A  Observed scores (%)", loc="left", pad=12)
    for i in range(7):
        for j in range(3):
            axes[0].text(j, i, f"{scores[i, j]:.2f}", ha="center", va="center", fontsize=10,
                         color="white" if scores[i, j] >= 55 else "black")
    for k, ax in enumerate(axes[1:]):
        ax.axvline(0, color="#888888", linewidth=1)
        for i, label in enumerate(PERTURBATIONS):
            value = report["comparisons"][label]["metrics"][METRICS[k]]
            color = "#007d83" if i < 3 else "#a66b0b"
            ax.plot(value["bonferroni_percentile_ci"], [i, i], color=color, linewidth=1)
            ax.plot(value["paired_percentile_ci"], [i, i], color=color, linewidth=4, solid_capstyle="butt")
            ax.plot(value["difference_percentage_points"], i, "o", color=color, markersize=4)
        ax.set_yticks(range(6), labels[1:] if k == 0 else [""] * 6)
        ax.set_ylim(5.6, -.6)
        ax.axhline(2.5, color="#dddddd", linewidth=1)
        ax.set_title(f"{'BCD'[k]}  {('F1', 'Precision', 'Recall')[k]} effects", loc="left", pad=12)
        ax.set_xlabel("Perturbation minus control (pp)", fontsize=9)
        ax.grid(axis="x", alpha=.2)
    for ax in axes:
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
    fig.text(.025, .16, "Teal: rooted clade distance 2. Ochre: distance 4. Labels retain the frozen variant indices.", fontsize=10)
    fig.text(.025, .11, "Thick: nominal 95% CI. Thin: Bonferroni CI across all 18 endpoints. 20,000 paired RefOG draws; seed 20260918.", fontsize=10)
    fig.text(.025, .06, "Development-exposed topology stress test, not posterior uncertainty, independent validation, or a species-tree selection procedure.", fontsize=10)
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
        path = args.output / ("species_tree_robustness." + extension)
        fig.savefig(path, dpi=180)
        outputs.append(file_provenance(path))
    plt.close(fig)
    (args.output / "manifest.json").write_text(json.dumps({"source_results": file_provenance(args.results),
        "plotter": file_provenance(Path(__file__)), "outputs": outputs,
        "matplotlib": matplotlib.__version__, "endpoints_shown": 18}, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
