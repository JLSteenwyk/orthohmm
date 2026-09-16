"""Plot frozen YGOB group-recovery scores and prespecified paired contrasts."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.run_simulation_methods import read_frozen

METHODS = ("orthohmm_satellite_v2", "orthohmm_high_sensitivity", "orthofinder_full", "orthofinder_sequence_only")
NAMES = ("OrthoHMM satellite_v2", "OrthoHMM high sensitivity", "OrthoFinder full", "OrthoFinder checkpoint\n(diagnostic)")
METRICS = ("f1", "precision", "recall")
COLORS = ("#343434", "#00858a", "#bc7510")


def plot(report):
    if report.get("frozen_evaluation_gates_verified") is not True or set(report["scores"]) != set(METHODS):
        raise ValueError("Expected admitted frozen four-method results")
    if report["uncertainty"]["replicates"] != 20000 or report["uncertainty"]["multiplicity_count"] != 6:
        raise ValueError("Unexpected uncertainty specification")
    scores = np.array([[report["scores"][m]["metrics"][k] * 100 for k in METRICS] for m in METHODS])
    if not np.isfinite(scores).all() or np.any((scores < 0) | (scores > 100)):
        raise ValueError("Invalid score")
    fig, (left, right) = plt.subplots(1, 2, figsize=(15, 7), gridspec_kw={"width_ratios": [1, 1.1]})
    fig.subplots_adjust(left=.17, right=.98, bottom=.23, top=.78, wspace=.48)
    fig.suptitle("Frozen YGOB validation: curated group recovery", x=.035, ha="left", y=.965, fontsize=17)
    fig.text(.035, .915, "16 species | 83,391 reference genes | 10,250 pillars | all methods: 100% reference-gene coverage", fontsize=10)
    for k, metric in enumerate(METRICS):
        positions = np.arange(4) + (k - 1) * .22
        left.barh(positions, scores[:, k], height=.18, color=COLORS[k], label=metric.upper() if k == 0 else metric.title())
        for y, score in zip(positions, scores[:, k]):
            left.text(score - .8, y, f"{score:.2f}", va="center", ha="right", fontsize=8, color="white")
    left.set_yticks(range(4), NAMES, fontsize=10)
    left.set_ylim(3.55, -.55)
    left.set_xlim(0, 100)
    left.set_xlabel("Score (%)")
    left.set_title("A  Observed scores", loc="left", fontsize=12, pad=15)
    left.legend(loc="upper left", bbox_to_anchor=(-.03, 1.24), ncol=3, frameon=False)
    right.axvline(0, color="#888888", linewidth=1)
    labels = []
    for j, method in enumerate(METHODS[:2]):
        for k, metric in enumerate(METRICS):
            y = j * 4 + k
            row = report["uncertainty"]["comparisons"][method][metric]
            point = row["difference_percentage_points"]
            adjusted, nominal = row["bonferroni_ci"], row["paired_95_percent_ci"]
            if not np.isfinite([point, *adjusted, *nominal]).all() or adjusted[0] > adjusted[1] or nominal[0] > nominal[1]:
                plt.close(fig)
                raise ValueError("Invalid interval")
            right.plot(adjusted, [y, y], color=COLORS[k], linewidth=1)
            right.plot(nominal, [y, y], color=COLORS[k], linewidth=4, solid_capstyle="butt")
            right.plot(point, y, "o", color=COLORS[k], markersize=5)
            labels.append((y, ("Satellite_v2" if j == 0 else "High sensitivity") + " / " + metric.upper()))
    right.set_yticks([y for y, _ in labels], [name for _, name in labels], fontsize=9)
    right.set_ylim(6.7, -.7)
    right.set_xlabel("OrthoHMM minus full OrthoFinder (percentage points)")
    right.set_title("B  Paired differences", loc="left", fontsize=12, pad=15)
    right.grid(axis="x", alpha=.2)
    for ax in (left, right):
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
    fig.text(.035, .14, "Thick: nominal 95% CI. Thin: Bonferroni CI across six endpoints. 20,000 paired pillar bootstrap replicates.", fontsize=10)
    fig.text(.035, .095, "Satellite_v2 F1 interval includes zero: neither superiority nor formal equivalence is established.", fontsize=10)
    fig.text(.035, .05, "Novel-taxon transfer, not family-disjoint validation. Co-membership includes within-species pairs; not resolved pairwise orthology.", fontsize=9)
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
        path = args.output / ("ygob_validation." + extension)
        fig.savefig(path, dpi=180)
        outputs.append(file_provenance(path))
    plt.close(fig)
    manifest = {"source_results": file_provenance(args.results), "plotter": file_provenance(Path(__file__)),
                "outputs": outputs, "matplotlib": matplotlib.__version__, "statistic": "frozen curated group co-membership"}
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
