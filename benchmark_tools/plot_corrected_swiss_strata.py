"""Render all prespecified corrected composition-stratum endpoints."""

import argparse
import csv
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from benchmark_tools.bootstrap_corrected_swiss_strata import METHODS, CONTRASTS
from benchmark_tools.bootstrap_qfo_swiss_stages import METRICS
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen


def endpoints(report):
    if (report.get("status") != "corrected_swiss_primary_stratified_intervals"
            or report.get("scientific_inputs_admitted") is not True
            or report.get("uncertainty_admitted") is not True
            or report.get("publication_ready") is not False
            or report.get("replicates") != 100000 or report.get("seed") != 20260924
            or report.get("multiplicity_endpoints") != 27
            or report.get("interaction_direction") != "higher contrast minus lower contrast"
            or report.get("units") != "raw 0-to-1 differences"):
        raise ValueError("Require admitted frozen corrected strata analysis")
    bins = report["bins"]
    if set(bins) != {"lower", "higher", "missing"}:
        raise ValueError("Changed bin inventory")
    families = [f for b in bins.values() for f in b["families"]]
    if len(families) != 18 or len(set(families)) != 18:
        raise ValueError("Changed or overlapping family universe")
    for name, value in bins.items():
        expected_eligible = name != "missing" and len(value["families"]) >= 5
        if value["interval_eligible"] is not expected_eligible:
            raise ValueError("Invalid interval eligibility")
    expected = [(METHODS[a], METHODS[b]) for a, b in CONTRASTS]
    for comparisons in [bins[b]["comparisons"] for b in bins] + [report["interactions"]]:
        if [(r["candidate"], r["reference"]) for r in comparisons] != expected:
            raise ValueError("Changed contrast inventory")
    rows = []
    for index, (candidate, reference) in enumerate(expected):
        for name in ("lower", "higher", "interaction"):
            eligible = (bins["lower"]["interval_eligible"] and bins["higher"]["interval_eligible"]
                        if name == "interaction" else bins[name]["interval_eligible"])
            comparison = (report["interactions"][index] if name == "interaction"
                          else bins[name]["comparisons"][index])
            if set(comparison["metrics"]) != set(METRICS):
                raise ValueError("Changed metric inventory")
            for metric in METRICS:
                value = comparison["metrics"][metric]
                point, nominal, adjusted = (value[k] for k in
                    ("difference", "paired_percentile_ci", "bonferroni_percentile_ci"))
                point_available = (bool(bins["lower"]["families"]) and bool(bins["higher"]["families"])
                                   if name == "interaction" else bool(bins[name]["families"]))
                if (point is not None) != point_available:
                    raise ValueError("Point availability differs from family inventory")
                if eligible and (point is None or nominal is None or adjusted is None):
                    raise ValueError("Eligible endpoint lacks interval")
                if not eligible and (nominal is not None or adjusted is not None):
                    raise ValueError("Unavailable interval must remain missing")
                bound = 2 if name == "interaction" else 1
                values = [] if point is None else [point]
                if nominal is not None:
                    if len(nominal) != 2 or len(adjusted) != 2:
                        raise ValueError("Malformed interval")
                    values += [*nominal, *adjusted]
                    if not adjusted[0] <= nominal[0] <= nominal[1] <= adjusted[1]:
                        raise ValueError("Invalid interval ordering")
                if any(type(v) not in (int, float) or not math.isfinite(v) or abs(v) > bound for v in values):
                    raise ValueError("Invalid endpoint value")
                rows.append(dict(candidate=candidate, reference=reference, stratum=name, metric=metric,
                                 difference=point, nominal=nominal, adjusted=adjusted))
    return rows


def plot(report):
    rows = endpoints(report)
    extent = max([abs(v) * 100 for row in rows for v in
                  ([row["difference"]] if row["difference"] is not None else []) +
                  (row["adjusted"] or [])] + [5])
    limit = math.ceil((extent + 2) / 10) * 10
    positions = (0, 1, 2, 4, 5, 6, 8, 9, 10)
    colors = ("#00858a", "#a33b45", "#555555")
    labels = ("High sensitivity minus full OrthoFinder", "Phylogeny minus full OrthoFinder",
              "Phylogeny minus high sensitivity")
    names = (f"Lower entropy ({len(report['bins']['lower']['families'])} families)",
             f"Higher entropy ({len(report['bins']['higher']['families'])} families)",
             "Interaction: higher minus lower")
    fig, axes = plt.subplots(1, 3, figsize=(16, 9))
    fig.subplots_adjust(left=.30, right=.98, bottom=.23, top=.79, wspace=.14)
    fig.suptitle("Corrected QfO SwissTrees: composition-stratified differences",
                 x=.025, y=.96, ha="left", fontsize=17)
    fig.text(.025, .90, "100,000 paired family resamples within bins | 27 planned endpoints | full OrthoFinder 3.1.5", fontsize=11)
    missing = len(report["bins"]["missing"]["families"])
    fig.text(.025, .85, f"18 reference families; missing-composition bin: {missing} (descriptive, not a primary interval).", fontsize=11)
    for column, (ax, metric) in enumerate(zip(axes, METRICS)):
        selected = [r for r in rows if r["metric"] == metric]
        ax.axvline(0, color="#999999", linestyle="--", linewidth=1)
        for index, (y, row) in enumerate(zip(positions, selected)):
            color = colors[index % 3]
            for key, width in (("adjusted", 1.2), ("nominal", 4)):
                if row[key] is not None:
                    ax.plot([100 * v for v in row[key]], [y, y], color=color, linewidth=width)
            if row["difference"] is not None:
                ax.plot(row["difference"] * 100, y, "o", color=color, markersize=5)
            if row["adjusted"] is None:
                ax.text(.98, y, "CI unavailable", transform=ax.get_yaxis_transform(), ha="right", fontsize=8)
        ax.set(xlim=(-limit, limit), ylim=(10.7, -1.2), xlabel="Difference (percentage points)")
        ax.set_yticks(positions, names * 3 if column == 0 else [""] * 9, fontsize=10)
        title = {"F1": "F1", "PPV": "Precision", "TPR": "Recall"}[metric]
        ax.set_title(f"{'ABC'[column]}  {title}", loc="left", fontsize=12, pad=12)
        ax.grid(axis="x", alpha=.15)
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
    for y, label in zip((-.8, 3.2, 7.2), labels):
        axes[0].text(-.02, y, label, transform=axes[0].get_yaxis_transform(), ha="right", fontsize=10, fontweight="bold")
    fig.text(.025, .16, "Points: observed differences. Thick lines: nominal 95% intervals. Thin lines: Bonferroni intervals across 27 endpoints.", fontsize=10)
    fig.text(.025, .10, "Development-exposed, conditional family bootstrap; composition strata are not causal explanations.", fontsize=10)
    fig.text(.025, .04, "Configuration contrasts are not pure phylogeny ablations. F1 is harmonic macro precision/recall; missing intervals are never zero.", fontsize=10)
    return fig, rows


def render(source, sha, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    identity = record(source)
    report = read_frozen(source, sha)
    fig, rows = plot(report)
    output.mkdir(parents=True)
    try:
        for extension in ("png", "pdf", "svg"):
            fig.savefig(output / ("corrected_swiss_strata." + extension), dpi=180)
    finally:
        plt.close(fig)
    with (output / "endpoints.tsv").open("x") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["candidate", "reference", "stratum", "metric", "difference", "nominal_lower",
                         "nominal_upper", "adjusted_lower", "adjusted_upper"])
        for row in rows:
            writer.writerow([row[k] for k in ("candidate", "reference", "stratum", "metric", "difference")]
                            + (row["nominal"] or [None, None]) + (row["adjusted"] or [None, None]))
    check(identity)
    manifest = dict(results=identity, source=record(__file__), matplotlib=matplotlib.__version__,
                    endpoints=27, display_conversion="raw differences multiplied by 100; TSV remains raw",
                    outputs=[record(p) for p in sorted(output.iterdir())], publication_ready=False)
    with (output / "manifest.json").open("x") as stream:
        json.dump(manifest, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--results-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    render(args.results.resolve(), args.results_sha256, args.output.absolute())
