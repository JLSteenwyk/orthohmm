"""Plot every prespecified SwissTrees comparator contrast from pinned intervals."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from benchmark_tools.audit_qfo_swiss_comparators import METHODS
from benchmark_tools.bootstrap_qfo_swiss_comparators import CONTRASTS, METRICS
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check

RESULT_SHA = "6317d0274142354dbc0745dc3c75a2532d0d39dedeba8f30b9aca8c6fee21005"
NAMES = ("OrthoHMM high sensitivity", "OrthoHMM phylogeny", "OrthoFinder MCL checkpoint*",
         "SonicParanoid", "Proteinortho", "FastOMA (supplied tree)*", "OrthoMCL",
         "OrthoHMM phylogeny\nminus high sensitivity")
POSITIONS = (0, 1, 2, 3, 4, 5, 6, 8)
COLORS = ("#00858a", "#a33b45", "#666666", "#303030", "#303030", "#666666", "#303030", "#a33b45")


def validate(report):
    expected = [(METHODS[a], METHODS[b]) for a, b in CONTRASTS]
    if (report["status"] != "paired_swiss_comparator_intervals" or
            [(r["candidate"], r["reference"]) for r in report["comparisons"]] != expected or
            report["replicates"] != 100000 or report["seed"] != 20260920 or
            report["multiplicity_endpoints"] != 24 or len(set(report["families"])) != 18):
        raise ValueError("Changed frozen contrast or resampling inventory")
    for row in report["comparisons"]:
        for metric in METRICS:
            values = row["metrics"][metric]
            nominal, adjusted = values["paired_percentile_ci"], values["bonferroni_percentile_ci"]
            if len(nominal) != 2 or len(adjusted) != 2 or not np.isfinite([values["difference"], *nominal, *adjusted]).all():
                raise ValueError("Invalid plotted interval")
            if not adjusted[0] <= nominal[0] <= nominal[1] <= adjusted[1]:
                raise ValueError("Invalid interval ordering")
            if min(*adjusted, values["difference"]) * 100 < -55 or max(*adjusted, values["difference"]) * 100 > 50:
                raise ValueError("Plotted evidence exceeds fixed axis limits")


def plot(report):
    validate(report)
    fig, axes = plt.subplots(1, 3, figsize=(14, 8.5))
    fig.subplots_adjust(left=.255, right=.98, bottom=.25, top=.80, wspace=.13)
    fig.suptitle("SwissTrees: paired method differences", x=.025, y=.965, ha="left", fontsize=18)
    fig.text(.025, .915, "18 curated families | 100,000 paired family resamples | 24 planned metric contrasts", fontsize=11)
    fig.text(.025, .865, "First seven rows: candidate minus full OrthoFinder 3.1.5. Last row: comparison of OrthoHMM configurations.", fontsize=10)
    for column, (ax, metric, title) in enumerate(zip(axes, METRICS, ("F1", "Precision", "Recall"))):
        ax.axvline(0, color="#999999", linewidth=1, linestyle="--")
        ax.axhline(7, color="#dddddd", linewidth=1)
        for y, color, row in zip(POSITIONS, COLORS, report["comparisons"]):
            values = row["metrics"][metric]
            adjusted = np.asarray(values["bonferroni_percentile_ci"]) * 100
            nominal = np.asarray(values["paired_percentile_ci"]) * 100
            ax.plot(adjusted, [y, y], color=color, linewidth=1.2)
            ax.plot(nominal, [y, y], color=color, linewidth=4, solid_capstyle="butt")
            ax.plot(values["difference"] * 100, y, "o", color=color, markersize=5)
        ax.set(xlim=(-55, 50), ylim=(8.7, -.7), xlabel="Difference (percentage points)")
        ax.set_xticks([-50, -25, 0, 25, 50])
        ax.set_yticks(POSITIONS, NAMES if column == 0 else [""] * 8, fontsize=10)
        ax.set_title(f"{'ABC'[column]}  {title}", loc="left", fontsize=12, pad=12)
        ax.grid(axis="x", alpha=.15)
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
    fig.text(.025, .17, "Thick lines: nominal 95% intervals. Thin lines: Bonferroni-adjusted intervals across 24 endpoints. Points: observed differences.", fontsize=10)
    fig.text(.025, .12, "Retrospective, development-exposed estimates. Shared evolutionary history and merged predictions can correlate families.", fontsize=10)
    fig.text(.025, .075, "* MCL checkpoint is diagnostic; FastOMA used a supplied OrthoFinder tree. OrthoHMM configurations differ beyond reconciliation.", fontsize=10)
    fig.text(.025, .03, "F1 is the harmonic mean of macro precision and recall. Intervals do not cover other QfO metrics or the six-metric secondary mean.", fontsize=10)
    return fig


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    identity = record(args.results)
    if identity["sha256"] != RESULT_SHA:
        raise ValueError("Changed frozen interval result")
    report = json.loads(args.results.read_text())
    if args.output.exists():
        raise FileExistsError(args.output)
    fig = plot(report)
    args.output.mkdir(parents=True)
    outputs = []
    for extension in ("png", "pdf", "svg"):
        path = args.output / ("swiss_comparator_intervals." + extension)
        fig.savefig(path, dpi=180)
        outputs.append(record(path))
    plt.close(fig)
    check(identity)
    manifest = {"results": identity, "plotter": record(__file__), "outputs": outputs,
                "matplotlib": matplotlib.__version__, "numpy": np.__version__,
                "display_conversion": "raw difference and interval bounds multiplied by 100",
                "contrasts": len(CONTRASTS), "metric_panels": list(METRICS)}
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
