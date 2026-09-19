"""Render every planned corrected SwissTrees contrast, including missing rows."""

import argparse
import csv
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.reproduce_corrected_swiss_comparison import verify, METHODS, METRICS
from benchmark_tools.run_simulation_methods import read_frozen

LABELS = ("OrthoHMM sensitive", "OrthoHMM phylogenetic", "OrthoFinder full",
          "OrthoFinder sequence-only", "SonicParanoid", "Proteinortho", "FastOMA", "OrthoMCL")


def endpoints(result):
    verify(result)
    rows = []
    for contrast in result["comparisons"]:
        for metric in METRICS:
            values = contrast["metrics"][metric] if contrast["status"] == "estimated" else None
            rows.append({"candidate": contrast["candidate"], "reference": contrast["reference"],
                "metric": metric, "status": contrast["status"],
                "difference": values["difference"] if values else None,
                "nominal_low": values["paired_percentile_ci"][0] if values else None,
                "nominal_high": values["paired_percentile_ci"][1] if values else None,
                "adjusted_low": values["bonferroni_percentile_ci"][0] if values else None,
                "adjusted_high": values["bonferroni_percentile_ci"][1] if values else None,
                **{key: values[key] if values else None for key in ("family_wins", "family_ties", "family_losses")}})
    return rows


def run(path, digest, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    result = read_frozen(path, digest)
    source = record(__file__)
    inputs = [record(path), source, record(Path(__file__).with_name("reproduce_corrected_swiss_comparison.py"))]
    for item in inputs:
        check(item)
    rows = endpoints(result)
    output.mkdir(parents=True, exist_ok=False)
    with (output / "endpoints.tsv").open("x") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    labels = dict(zip(METHODS, LABELS))
    descriptions = [labels[r["candidate"]] + "\nminus " + labels[r["reference"]] for r in result["comparisons"]]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "svg.hashsalt": "corrected-swiss-comparison"})
    fig, axes = plt.subplots(1, 3, figsize=(15, 7), sharey=True)
    colors = ("#167d8d", "#bd4148", "#9a710c")
    for axis, metric, title, color in zip(axes, METRICS, ("F1", "Precision", "Recall"), colors):
        axis.axvline(0, color="#777777", linewidth=.8, linestyle="--")
        for y, contrast in enumerate(result["comparisons"]):
            if contrast["status"] != "estimated":
                axis.text(0, y, "Not admitted", ha="center", va="center", color="#777777", fontsize=9,
                          bbox={"facecolor": "white", "edgecolor": "none", "pad": 2})
                continue
            values = contrast["metrics"][metric]
            lo, hi = [100 * x for x in values["bonferroni_percentile_ci"]]
            nl, nh = [100 * x for x in values["paired_percentile_ci"]]
            axis.plot([lo, hi], [y, y], color=color, linewidth=1.2)
            axis.plot([nl, nh], [y, y], color=color, linewidth=4, solid_capstyle="butt")
            axis.plot(100 * values["difference"], y, "o", color=color, markersize=5)
        axis.set_title(title, fontweight="bold")
        axis.set_xlim(-60, 60)
        axis.set_xticks([-60, -30, 0, 30, 60])
        axis.set_xlabel("Difference (percentage points)")
        axis.grid(axis="y", color="#eeeeee", linewidth=.7)
        axis.spines[["top", "right", "left"]].set_visible(False)
        axis.tick_params(axis="y", length=0)
    axes[0].set_yticks(range(8), descriptions)
    axes[0].set_ylim(7.65, -.65)
    fig.suptitle("Corrected QfO: SwissTrees paired comparisons", x=.03, ha="left", fontsize=16)
    fig.text(.03, .92, "18 families; 100,000 shared draws. Thick: nominal 95% interval. Thin: 24-endpoint adjusted interval.", fontsize=10)
    fig.text(.03, .025, "Development-exposed evidence. Missing methods retained. Sequence-only OrthoFinder is a diagnostic; configuration contrasts are not pure ablations.", fontsize=9)
    fig.subplots_adjust(left=.24, right=.985, top=.85, bottom=.13, wspace=.16)
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(output / ("corrected_swiss_comparison." + suffix), dpi=180)
    plt.close(fig)
    lines = ["# Corrected SwissTrees Paired Contrasts", "",
        "Candidate minus reference, raw 0-to-1 units. Adjustment retains all 24 planned endpoints.", "",
        "| Candidate | Reference | Metric | Difference | Adjusted interval | Wins/ties/losses |",
        "| --- | --- | --- | ---: | --- | --- |"]
    for row in rows:
        cells = [labels[row["candidate"]], labels[row["reference"]], row["metric"]]
        if row["status"] == "estimated":
            cells.extend([f'{row["difference"]:+.6f}', f'[{row["adjusted_low"]:.6f}, {row["adjusted_high"]:.6f}]',
                          f'{row["family_wins"]}/{row["family_ties"]}/{row["family_losses"]}'])
        else:
            cells.extend(["Not admitted", "Not estimable", "Not estimable"])
        lines.append("| " + " | ".join(cells) + " |")
    (output / "endpoints.md").write_text("\n".join(lines) + "\n")
    for item in inputs:
        check(item)
    manifest = {"status": "corrected_swiss_comparison_rendered", "inputs": inputs,
        "endpoints": len(rows), "display_scale": 100, "tsv_scale": 1,
        "outputs": [record(p) for p in sorted(output.iterdir())], "publication_ready": False}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--results-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.results.resolve(), args.results_sha256, args.output.absolute())
