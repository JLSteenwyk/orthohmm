"""Plot descriptive overlap strata without implying independent confirmation."""

import argparse
import csv
from fractions import Fraction
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from benchmark_tools.plot_ygob_validation import METHODS, NAMES, METRICS
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.run_simulation_methods import read_frozen

STRATA = ("screen_positive", "screen_negative")
COLORS = ("#009E73", "#0072B2", "#D55E00", "#E69F00")
MARKERS = ("s", "o", "^", "v")
SIZES = ("reference_groups", "reference_genes", "truth_pairs", "zero_truth_pair_pillars")


def plotting_data(report):
    if (report.get("status") != "ygob_overlap_strata_descriptive"
            or report.get("independent_confirmation") is not False
            or report.get("confidence_intervals_computed") is not False
            or set(report["strata"]) != set(STRATA)):
        raise ValueError("Require descriptive two-stratum results without confirmation or CI claims")
    rows = []
    for stratum in STRATA:
        group = report["strata"][stratum]
        if set(group) != set(METHODS):
            raise ValueError("Require all four frozen methods")
        signature = None
        for method, label, color, marker in zip(METHODS, NAMES, COLORS, MARKERS):
            item = group[method]
            counts = item["counts"]
            if (any(type(counts[k]) is not int or counts[k] < 0 for k in ("tp", "fn"))
                    or type(counts["fp"]) not in (float, int) or not math.isfinite(counts["fp"])
                    or counts["fp"] < 0 or counts["fp"] * 2 != int(counts["fp"] * 2)):
                raise ValueError("Require nonnegative integer TP/FN and half-integer FP")
            sizes = {key: item[key] for key in SIZES}
            if any(type(value) is not int or value < 0 for value in sizes.values()):
                raise ValueError("Invalid reference sizes")
            if (counts["tp"] + counts["fn"] != sizes["truth_pairs"]
                    or sizes["zero_truth_pair_pillars"] > sizes["reference_groups"]):
                raise ValueError("Inconsistent truth counts or singleton count")
            if signature is not None and signature != sizes:
                raise ValueError("Reference sizes differ across methods")
            signature = sizes
            tp, fp, fn = (Fraction(counts[key]) for key in ("tp", "fp", "fn"))
            ratios = dict(f1=(2*tp, 2*tp+fp+fn), precision=(tp, tp+fp), recall=(tp, tp+fn))
            for metric in METRICS:
                numerator, denominator = ratios[metric]
                defined = denominator != 0
                score = float(numerator / denominator) if defined else 0.0
                actual = item["metrics"][metric]
                if (type(actual) not in (int, float) or not math.isfinite(actual)
                        or not math.isclose(actual, score, rel_tol=0, abs_tol=1e-12)
                        or item["defined"][metric] is not defined):
                    raise ValueError("Metric or definition flag differs from sufficient statistics")
                rows.append(dict(stratum=stratum, method=method, label=label, color=color,
                    marker=marker, metric=metric, score_percent=score*100, defined=defined,
                    **counts, **sizes))
    return rows


def plot(rows):
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.4), sharex=True, sharey=True)
    fig.subplots_adjust(left=.22, right=.98, top=.77, bottom=.19, hspace=.87, wspace=.17)
    fig.suptitle("YGOB overlap screen: group-recovery trade-offs", x=.035, ha="left", y=.965, fontsize=17)
    fig.text(.035, .914, "Secondary descriptive analysis | all four frozen methods | no new confidence intervals", fontsize=11)
    indexed = {(r["stratum"], r["method"], r["metric"]): r for r in rows}
    for i, stratum in enumerate(STRATA):
        exemplar = indexed[stratum, METHODS[0], METRICS[0]]
        heading = "Screen-positive" if i == 0 else "Screen-negative (no qualifying hit)"
        size_text = (f"{exemplar['reference_groups']:,} pillars; {exemplar['reference_genes']:,} reference genes; "
                     f"{exemplar['truth_pairs']:,} truth pairs; "
                     f"{exemplar['zero_truth_pair_pillars']:,} singleton pillars")
        top = axes[i, 0].get_position().y1
        fig.text(.22, top+.083, heading, fontsize=12, fontweight="bold")
        fig.text(.22, top+.053, size_text, fontsize=10)
        for j, metric in enumerate(METRICS):
            ax = axes[i, j]
            for y, method in enumerate(METHODS):
                row = indexed[stratum, method, metric]
                if row["defined"]:
                    value = row["score_percent"]
                    ax.scatter(value, y, marker=row["marker"], color=row["color"], s=65, zorder=3)
                    ax.annotate(f"{value:.2f}", (value, y), xytext=(-8 if value > 90 else 8, 0),
                                textcoords="offset points", ha="right" if value > 90 else "left",
                                va="center", fontsize=9)
                else:
                    ax.text(.04, y, "undefined", transform=ax.get_yaxis_transform(), va="center", fontsize=9)
            title = "F1" if metric == "f1" else metric.title()
            ax.set_title(f"{chr(65+i*3+j)}  {title}", loc="left", fontsize=11, pad=7)
            ax.set(xlim=(0, 100), ylim=(3.6, -.6), xticks=[0, 25, 50, 75, 100],
                   yticks=range(4), yticklabels=NAMES)
            ax.tick_params(axis="both", length=0, labelsize=9, labelbottom=True)
            ax.grid(axis="x", alpha=.2, linewidth=.7)
            ax.set_axisbelow(True)
            for spine in ax.spines.values():
                spine.set_visible(False)
            if i == 1:
                ax.set_xlabel("Score (%)", fontsize=10)
    fig.text(.035, .104, "No-hit status does not establish family independence or absence of remote homology.", fontsize=10)
    fig.text(.035, .063, "Original cross-pillar false-positive allocations retained; predictions were not subset-rescored.", fontsize=10)
    fig.text(.035, .025, "Curated group co-membership includes within-species pairs; these are not resolved pairwise orthology scores.", fontsize=9)
    return fig


def render(source, sha256, output):
    rows = plotting_data(read_frozen(source, sha256))
    output.mkdir(parents=True, exist_ok=False)
    table = output / "plotted_values.tsv"
    with table.open("x") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    with plt.rc_context({"font.family": "DejaVu Sans", "svg.hashsalt": "ygob-overlap-strata-v1"}):
        fig = plot(rows)
        try:
            outputs = [file_provenance(table)]
            for extension in ("png", "pdf", "svg"):
                path = output / ("ygob_overlap_strata." + extension)
                metadata = {"CreationDate": None, "ModDate": None} if extension == "pdf" else {"Date": None} if extension == "svg" else None
                fig.savefig(path, dpi=180, metadata=metadata)
                outputs.append(file_provenance(path))
        finally:
            plt.close(fig)
    manifest = dict(source_results=file_provenance(source), plotter=file_provenance(Path(__file__)),
        outputs=outputs, matplotlib=matplotlib.__version__, rows=len(rows),
        independent_confirmation=False, confidence_intervals_computed=False,
        statistic="Ratios of summed original pillar sufficient statistics; original allocated FP retained")
    if manifest["source_results"]["sha256"] != sha256:
        raise ValueError("Source changed during figure generation")
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    render(args.results.resolve(), args.sha256, args.output.absolute())
