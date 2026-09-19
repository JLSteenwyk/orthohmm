"""Render an explicitly identified SwissTrees factorial with all 42 endpoints."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from benchmark_tools.bootstrap_qfo_factorial import CELLS, PROTOCOL_SHA, contrasts
from benchmark_tools.bootstrap_qfo_corrected_factorial import CORRECTED_PROTOCOL_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen

METRICS = ("F1", "PPV", "TPR")


def validate(report, input_release="original"):
    if input_release not in ("original", "corrected"):
        raise ValueError("Unknown input release")
    status = "paired_qfo_factorial_swiss_intervals"
    if input_release == "corrected":
        status = "paired_corrected_qfo_factorial_swiss_intervals"
        if (report.get("input_release") != "QfO 2020_04 corrected UP000008143"
                or report.get("corrected_protocol", {}).get("sha256") != CORRECTED_PROTOCOL_SHA):
            raise ValueError("Changed corrected release identity or protocol")
    elif "corrected_protocol" in report or "input_release" in report:
        raise ValueError("Unexpected release annotation on historical results")
    expected = {"status": status, "replicates": 100000,
                "seed": 20260922, "alpha": .05, "multiplicity_endpoints": 42,
                "units": "raw 0-to-1 metric units", "quantile_method": "linear"}
    if any(report.get(k) != v for k, v in expected.items()):
        raise ValueError("Changed frozen uncertainty protocol")
    if report["protocol"]["sha256"] != PROTOCOL_SHA:
        raise ValueError("Changed protocol checksum")
    if len(report["families"]) != 18 or len(set(report["families"])) != 18:
        raise ValueError("Require 18 distinct families")
    if set(report["point_estimates"]) != set(CELLS):
        raise ValueError("Require all eight cells")
    scores = np.asarray([[report["point_estimates"][cell][m] for m in METRICS] for cell in CELLS])
    if not np.isfinite(scores).all() or np.any((scores < 0) | (scores > 1)):
        raise ValueError("Invalid raw metric units")
    if not np.allclose(scores[:, 0], 2 * scores[:, 1] * scores[:, 2] /
                       (scores[:, 1] + scores[:, 2]), atol=1e-12, rtol=0):
        raise ValueError("F1 differs from harmonic macro precision and recall")
    expected_contrasts = contrasts()
    if len(report["comparisons"]) != len(expected_contrasts):
        raise ValueError("Incomplete contrast inventory")
    for row, contrast in zip(report["comparisons"], expected_contrasts):
        if any(row.get(k) != v for k, v in contrast.items()):
            raise ValueError("Changed contrast identity")
        for j, metric in enumerate(METRICS):
            value = row["metrics"][metric]
            point = value["difference"]
            nominal, adjusted = (value[k] for k in ("paired_percentile_ci", "bonferroni_percentile_ci"))
            if (len(nominal) != 2 or len(adjusted) != 2
                    or not np.isfinite([point, *nominal, *adjusted]).all()
                    or not adjusted[0] <= nominal[0] <= nominal[1] <= adjusted[1]):
                raise ValueError("Invalid or nonnested intervals")
            if not np.isclose(point, np.asarray(contrast["weights"]) @ scores[:, j], atol=1e-12, rtol=0):
                raise ValueError("Effect differs from cell estimates")
            counts = [value[k] for k in ("family_wins", "family_ties", "family_losses")]
            if any(type(v) is not int or v < 0 for v in counts) or sum(counts) != 18:
                raise ValueError("Invalid family direction counts")
    return scores * 100


def plot(report, input_release="original"):
    scores = validate(report, input_release)
    fig, axes = plt.subplots(1, 4, figsize=(19, 9), gridspec_kw={"width_ratios": [1.1, 1, 1, 1]})
    fig.subplots_adjust(left=.055, right=.985, top=.79, bottom=.26, wspace=.7)
    fig.suptitle("QfO SwissTrees: HMM refinement, candidate expansion and reconciliation", x=.025, ha="left", y=.97, fontsize=17)
    release_label = ("Corrected-release inputs (984,137 genes)" if input_release == "corrected"
                     else "Original-release inputs (976,504 genes)")
    fig.text(.025, .915, release_label + " | 18 reference families | development-exposed analysis", fontsize=11)
    fig.text(.025, .865, "P = profile refinement; C = candidate expansion; R = reconciliation. Cells show P C R (0: off, 1: on).", fontsize=10)
    left = axes[0]
    left.imshow(scores, vmin=0, vmax=100, cmap="Greys", aspect="auto")
    left.set_xticks(range(3), ["F1", "Precision", "Recall"])
    left.set_yticks(range(8), [" ".join(c[-1] for c in cell.split("_")) for cell in CELLS])
    left.set_title("A  Observed scores (%)", loc="left", pad=12)
    for i in range(8):
        for j in range(3):
            left.text(j, i, f"{scores[i, j]:.2f}", ha="center", va="center", color="white" if scores[i, j] >= 55 else "black", fontsize=10)
    colors = ("#007d83", "#a66b0b", "#755297", "#333333")
    for k, ax in enumerate(axes[1:]):
        ax.axvline(0, color="#999999", linewidth=1)
        labels = []
        for i, row in enumerate(report["comparisons"]):
            value = row["metrics"][METRICS[k]]
            color = colors[min(i // 4, 3)]
            ax.plot(np.asarray(value["bonferroni_percentile_ci"]) * 100, [i, i], color=color, linewidth=1)
            ax.plot(np.asarray(value["paired_percentile_ci"]) * 100, [i, i], color=color, linewidth=4, solid_capstyle="butt")
            ax.plot(value["difference"] * 100, i, "o", color=color, markersize=4)
            labels.append(row["name"].replace("_by_", " x ").replace("_at_", " | ").replace("_", " "))
        ax.set_yticks(range(14), labels if k == 0 else [""] * 14, fontsize=9)
        ax.set_ylim(13.6, -.6)
        for boundary in (3.5, 7.5, 11.5):
            ax.axhline(boundary, color="#dddddd", linewidth=1)
        ax.set_xlabel("Effect (percentage points)", fontsize=9)
        ax.set_title(f"{'BCD'[k]}  {('F1', 'Precision', 'Recall')[k]} effects", loc="left", pad=12)
        ax.grid(axis="x", alpha=.2)
    for ax in axes:
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
    nonzero = sum(not (r["metrics"]["F1"]["bonferroni_percentile_ci"][0] <= 0 <=
                      r["metrics"]["F1"]["bonferroni_percentile_ci"][1]) for r in report["comparisons"])
    conclusion = (f"{nonzero}/14 adjusted F1 intervals exclude zero." if nonzero
                  else "All adjusted F1 intervals include zero.")
    release_note = ("No independent validation, selection adjustment, or intervals for other QfO endpoints or the secondary mean."
                    if input_release == "corrected" else
                    "No independent validation, equivalence claim, other-QfO-endpoint intervals, or corrected-input results are shown.")
    notes = [
        "Rows: P (teal), C (ochre), R (purple): on minus off, holding other factors fixed. Final rows: C x R interactions (gray).",
        "Thick: nominal 95% CI. Thin: Bonferroni CI across 42 endpoints. 100,000 paired family draws; seed 20260922.",
        "F1 is the harmonic mean of macro precision and recall, recomputed per draw. " + conclusion,
        "Profile-off retains initial HMM search. R changes group-derived pairs to native inferred pairs, not only group splitting.",
        release_note,
    ]
    for y, note in zip((.195, .155, .115, .075, .035), notes):
        fig.text(.025, y, note, fontsize=10)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--input-release", choices=("original", "corrected"), default="original")
    args = parser.parse_args()
    report = read_frozen(args.results, args.sha256)
    source, plotter = record(args.results), record(__file__)
    if args.output.exists():
        raise FileExistsError(args.output)
    fig = plot(report, args.input_release)
    args.output.mkdir(parents=True)
    outputs = []
    for extension in ("png", "pdf", "svg"):
        path = args.output / ("qfo_factorial_swiss." + extension)
        fig.savefig(path, dpi=180)
        outputs.append(record(path))
    plt.close(fig)
    check(source)
    check(plotter)
    manifest = {"source_results": source, "plotter": plotter, "outputs": outputs,
                "matplotlib": matplotlib.__version__, "endpoints_shown": 42,
                "input_release": args.input_release, "publication_ready": False,
                "helpers": [record(Path(__file__).with_name(name)) for name in (
                    "bootstrap_qfo_factorial.py", "bootstrap_qfo_corrected_factorial.py")]}
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
