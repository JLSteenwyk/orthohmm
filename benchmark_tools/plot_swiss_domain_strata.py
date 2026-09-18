"""Plot all frozen domain-stratified SwissTrees contrasts without recomputation."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from benchmark_tools.analyze_swiss_domain_strata import CONTRASTS, METHODS, METRICS, PRIMARY
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check

RESULT_SHA = "229269f3db0bfb0e99893de00a1c2467285a07d7ca148c4854b7661ca13a5891"
POSITIONS = (0, 1, 2, 4, 5, 6, 8, 9, 10)
LABELS = ("High sensitivity minus full OrthoFinder",
          "Phylogeny minus full OrthoFinder", "Phylogeny minus high sensitivity")
BIN_LABELS = ("Median Pfam types < 2 (12 families)",
              "Median Pfam types >= 2 (6 families)", "Interaction: higher minus lower")
COLORS = ("#00858a", "#a33b45", "#555555")


def rows(report):
    for contrast in report["contrasts"]:
        yield from (row["metrics"] for row in contrast["primary_strata"])
        yield contrast["interaction_higher_minus_lower"]


def validate(report):
    if (report["status"] != "frozen_swiss_domain_strata_evaluated" or
            report["replicates"] != 100000 or report["seed"] != 20260921 or
            report["multiplicity_endpoints"] != 27 or
            [(r["candidate"], r["reference"]) for r in report["contrasts"]] !=
            [(METHODS[a], METHODS[b]) for a, b in CONTRASTS]):
        raise ValueError("Changed frozen analysis inventory")
    low, high = (report["strata"][name] for name in PRIMARY)
    if len(low) != 12 or len(high) != 6 or len(set(low + high)) != 18:
        raise ValueError("Changed primary family bins")
    for contrast in report["contrasts"]:
        if [row["stratum"] for row in contrast["primary_strata"]] != list(PRIMARY):
            raise ValueError("Changed stratum order")
    for row in rows(report):
        for metric in METRICS:
            values = row[metric]
            nominal, adjusted = values["nominal95"], values["bonferroni27"]
            if len(nominal) != 2 or len(adjusted) != 2 or not np.isfinite([values["difference"], *nominal, *adjusted]).all():
                raise ValueError("Invalid interval")
            if not adjusted[0] <= nominal[0] <= nominal[1] <= adjusted[1]:
                raise ValueError("Invalid interval ordering")
            if min(values["difference"], *adjusted) < -.7 or max(values["difference"], *adjusted) > .7:
                raise ValueError("Evidence exceeds fixed axes")


def plot(report):
    validate(report)
    fig, axes = plt.subplots(1, 3, figsize=(16, 9))
    fig.subplots_adjust(left=.30, right=.98, bottom=.24, top=.79, wspace=.13)
    fig.suptitle("SwissTrees: domain-stratified method differences", x=.025, y=.96, ha="left", fontsize=18)
    fig.text(.025, .905, "18 curated families | 100,000 paired family resamples within bins | 27 planned endpoints", fontsize=11)
    fig.text(.025, .855, "OrthoHMM high sensitivity and phylogeny configurations; reference: full OrthoFinder 3.1.5.", fontsize=11)
    for column, (ax, metric, title) in enumerate(zip(axes, METRICS, ("F1", "Precision", "Recall"))):
        ax.axvline(0, color="#999999", linestyle="--", linewidth=1)
        for index, (y, row) in enumerate(zip(POSITIONS, rows(report))):
            values, color = row[metric], COLORS[index % 3]
            ax.plot(np.asarray(values["bonferroni27"]) * 100, [y, y], color=color, linewidth=1.2)
            ax.plot(np.asarray(values["nominal95"]) * 100, [y, y], color=color, linewidth=4, solid_capstyle="butt")
            ax.plot(values["difference"] * 100, y, "o", color=color, markersize=5)
        ax.set(xlim=(-70, 70), ylim=(10.7, -1.2), xlabel="Difference (percentage points)")
        ax.set_xticks([-60, -30, 0, 30, 60])
        ax.set_yticks(POSITIONS, BIN_LABELS * 3 if column == 0 else [""] * 9, fontsize=10)
        ax.set_title(f"{'ABC'[column]}  {title}", loc="left", fontsize=12, pad=12)
        ax.grid(axis="x", alpha=.15)
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
    for y, label in zip((-.8, 3.2, 7.2), LABELS):
        axes[0].text(-.02, y, label, transform=axes[0].get_yaxis_transform(),
                     ha="right", fontsize=10, fontweight="bold")
    fig.text(.025, .175, "Thick lines: nominal 95% intervals. Thin lines: Bonferroni-adjusted intervals across all 27 endpoints. Points: observed differences.", fontsize=10)
    fig.text(.025, .13, "All nine adjusted interaction intervals include zero; different within-bin results do not establish a difference between bins.", fontsize=10)
    fig.text(.025, .085, "Retrospective, development-exposed analysis; possible family dependence. Domain annotations do not establish causal mechanisms.", fontsize=10)
    fig.text(.025, .04, "OrthoHMM configurations differ beyond reconciliation. F1 is harmonic macro precision/recall. Repeat-type bins are descriptive only and not plotted.", fontsize=10)
    return fig


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    identity = record(args.results)
    if identity["sha256"] != RESULT_SHA:
        raise ValueError("Changed frozen result")
    if args.output.exists():
        raise FileExistsError(args.output)
    fig = plot(json.loads(args.results.read_text()))
    args.output.mkdir(parents=True)
    outputs = []
    for extension in ("png", "pdf", "svg"):
        path = args.output / ("swiss_domain_strata." + extension)
        fig.savefig(path, dpi=180)
        outputs.append(record(path))
    plt.close(fig)
    check(identity)
    manifest = {"results": identity, "plotter": record(__file__), "outputs": outputs,
                "matplotlib": matplotlib.__version__, "numpy": np.__version__,
                "display_conversion": "raw differences and bounds multiplied by 100",
                "endpoints": 27, "metric_panels": list(METRICS)}
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
