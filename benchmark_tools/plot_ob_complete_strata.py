"""Plot all retained OrthoBench descriptive strata without adding inference."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from benchmark_tools.plot_ob_stratified_errors import ORDER, LABELS
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.run_simulation_methods import read_frozen

SOURCE_SHA = "1850a897bef81e7a5a2829dc6770388b4936d767b4e41f1e704ee8ffc4093811"
METHODS = ("orthohmm_high_sensitivity", "orthohmm_phylogeny_satellite_v2",
           "orthofinder_3_1_5_full", "orthofinder_3_1_5_sequence_only",
           "sonicparanoid_2_0_9", "proteinortho_6_3_6", "fastoma_0_3_5", "orthomcl_1_4")
NAMES = ("OH sensitive", "OH phylogenetic", "OF full", "OF sequence-only",
         "SonicParanoid", "ProteinOrtho", "FastOMA", "OrthoMCL")
METRICS = ("f_score", "precision", "recall")


def matrices(report):
    if (report["status"] != "all_method_orthobench_descriptive_strata"
            or report["new_intervals_calculated"] is not False):
        raise ValueError("Require descriptive export")
    rows = report["rows"]
    indexed = {(r["stratum"], r["method"]): r for r in rows}
    if len(rows) != 112 or set(indexed) != {(s, m) for s in ORDER for m in METHODS}:
        raise ValueError("Require exactly one row for every method/bin")
    values = np.full((3, 14, 8), np.nan)
    counts = []
    for i, label in enumerate(ORDER):
        first = indexed[label, METHODS[0]]
        counts.append(first["family_count"])
        for j, method in enumerate(METHODS):
            row = indexed[label, method]
            if (row["families"] != first["families"] or row["family_count"] != len(row["families"])
                    or len(set(row["families"])) != row["family_count"]):
                raise ValueError("Changed family membership")
            point = row["metrics_percent"]
            if not row["family_count"]:
                if point is not None or row["status"] != "empty_nonestimable":
                    raise ValueError("Empty bin must be missing")
            else:
                if set(point) != set(METRICS) or row["status"] != "descriptive":
                    raise ValueError("Invalid metrics or status")
                scores = [point[k] for k in METRICS]
                if not np.isfinite(scores).all() or any(not 0 <= v <= 100 for v in scores):
                    raise ValueError("Scores must be finite percentages")
                values[:, i, j] = scores
    return values, counts


def plot(report):
    values, counts = matrices(report)
    fig, axes = plt.subplots(1, 3, figsize=(22, 10), sharey=True)
    fig.subplots_adjust(left=.19, right=.94, bottom=.25, top=.86, wspace=.08)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#eeeeee")
    for k, ax in enumerate(axes):
        mesh = ax.imshow(values[k], vmin=0, vmax=100, cmap=cmap, aspect="auto")
        ax.set_title(("A  Weighted F1", "B  Precision", "C  Recall")[k], loc="left", fontsize=13)
        ax.set_xticks(range(8), NAMES, rotation=55, ha="right", fontsize=9)
        ax.set_yticks(range(14), [f"{label} (n={n})" for label, n in zip(LABELS, counts)], fontsize=10)
        ax.tick_params(length=0)
        for i in range(14):
            for j in range(8):
                value = values[k, i, j]
                ax.text(j, i, "NA" if np.isnan(value) else f"{value:.1f}", ha="center", va="center",
                        fontsize=8, color="#333333" if np.isnan(value) or value >= 55 else "white")
        for boundary in (2.5, 4.5, 7.5, 10.5):
            ax.axhline(boundary, color="white", linewidth=2)
        for spine in ax.spines.values():
            spine.set_visible(False)
    cax = fig.add_axes([.955, .25, .009, .61])
    fig.colorbar(mesh, cax=cax, label="Percent")
    fig.suptitle("OrthoBench: all-method descriptive strata", x=.03, ha="left", y=.96, fontsize=20)
    fig.text(.03, .915, "Eight retained methods | 14 frozen bins | 70 reference families | fixed method order", fontsize=12)
    for y, text in zip((.10, .068, .036), (
            "Development-exposed descriptive results; no new confidence intervals or significance tests. Empty bins are NA, not zero.",
            "Statistics recomputed from family-size-weighted TP/FP/FN, not mean family F1. Bins overlap across dimensions; n counts reference families.",
            "OH = OrthoHMM; OF = OrthoFinder 3.1.5. Descriptors are not validated domain, fragmentation or duplication-history annotations.")):
        fig.text(.03, y, text, fontsize=11)
    return fig


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = read_frozen(args.source, SOURCE_SHA)
    args.output.mkdir(parents=True, exist_ok=False)
    fig = plot(report)
    for extension in ("png", "pdf", "svg"):
        fig.savefig(args.output / ("ob_complete_strata." + extension), dpi=160)
    plt.close(fig)
    receipt = dict(source=record(args.source.resolve()), plotter=record(Path(__file__).resolve()),
        displayed_finite_values=int(np.isfinite(matrices(report)[0]).sum()),
        missing_values=int(np.isnan(matrices(report)[0]).sum()), new_intervals=False,
        outputs=[record(p.resolve()) for p in sorted(args.output.glob("ob_complete_strata.*"))])
    (args.output / "manifest.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt))
