"""Generate development-evidence figures from the audited comparison JSON."""

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.publication_comparison import METHODS


LABELS = (
    "OrthoHMM sensitive", "OrthoHMM satellite_v2", "OrthoFinder 3.1.5 full",
    "OrthoFinder 3.1.5 sequence", "SonicParanoid 2.0.9", "ProteinOrtho 6.3.6",
    "FastOMA 0.3.5 (supplied tree)", "OrthoMCL 1.4",
)
COLORS = ("#0072B2", "#009E73", "#D55E00", "#E69F00", "#CC79A7", "#6B6B6B", "#5647A6", "#111111")
MARKERS = ("o", "s", "^", "v", "D", "P", "X", "*")
QFO_AXES = {
    "VGNC F": ("TPR", "PPV"), "SwissTrees F": ("TPR", "PPV"),
    "TreeFam-A F": ("TPR", "PPV"), "EC": ("NR_ORTHOLOGS", "avg Schlicker"),
    "GO": ("NR_ORTHOLOGS", "avg Schlicker"), "FAS": ("NR_ORTHOLOGS", "FAS"),
}


def bounded(value, upper=1):
    if not isinstance(value, (int, float)) or not math.isfinite(value) or not 0 <= value <= upper:
        raise ValueError("Invalid plotted metric")
    return value


def plotting_data(report):
    rows = report["methods"]
    indexed = {r["key"]: r for r in rows}
    if len(indexed) != len(rows) or set(indexed) != {m[0] for m in METHODS}:
        raise ValueError("Comparison must contain exactly the retained method panel")
    data = []
    for (key, _, _, _), label, color, marker in zip(METHODS, LABELS, COLORS, MARKERS):
        row = indexed[key]
        ob = row["orthobench"]
        record = {"key": key, "label": label, "color": color, "marker": marker,
                  "orthobench": {k: bounded(ob[k + "_percent"], 100) for k in ("f_score", "precision", "recall")},
                  "three_kingdoms_f1": bounded(row["three_kingdoms"]["score"]["f_score"]),
                  "qfo_status": row["qfo"]["status"], "qfo": {}, "uncertainty": None}
        if row["qfo"]["status"] == "metrics_available":
            for metric, expected_axes in QFO_AXES.items():
                item = row["qfo"]["metric_details"][metric]
                if tuple(item["axes"][k] for k in ("x_axis", "y_axis")) != expected_axes:
                    raise ValueError(f"Unexpected QfO axes for {metric}")
                p = item["participant"]
                record["qfo"][metric] = {
                    "x": bounded(p["metric_x"], float("inf") if expected_axes[0] == "NR_ORTHOLOGS" else 1),
                    "y": bounded(p["metric_y"]),
                }
        if key.startswith("orthohmm_"):
            uncertainty = ob["uncertainty"]
            if uncertainty["versus"] != "orthofinder_3_1_5_full":
                raise ValueError("Unexpected paired comparator")
            for metric in ("f_score", "precision", "recall"):
                point = record["orthobench"][metric] - indexed["orthofinder_3_1_5_full"]["orthobench"][metric + "_percent"]
                values = uncertainty["metrics"][metric]
                if not math.isclose(point, values["difference_percentage_points"], abs_tol=1e-8):
                    raise ValueError("Paired estimate disagrees with point estimates")
                for interval in ("paired_percentile_ci", "bonferroni_percentile_ci"):
                    limits = values[interval]
                    if len(limits) != 2 or not all(math.isfinite(v) for v in limits) or limits[0] > limits[1]:
                        raise ValueError("Invalid paired confidence interval")
            record["uncertainty"] = uncertainty["metrics"]
        data.append(record)
    return data


def style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(alpha=0.18, linewidth=0.6)
    ax.set_axisbelow(True)


def legend(fig, data, y=0.015):
    handles = [Line2D([], [], color=r["color"], marker=r["marker"], linestyle="none",
                      markersize=7, label=r["label"]) for r in data]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, y),
               ncol=4, frameon=False, fontsize=9)


def overview(data):
    fig, axes = plt.subplots(1, 2, figsize=(12, 6.4), gridspec_kw={"width_ratios": [1, 1.35]})
    for r in data:
        axes[0].scatter(r["orthobench"]["recall"], r["orthobench"]["precision"],
                        c=r["color"], marker=r["marker"], s=70, edgecolors="white", linewidths=0.5)
    axes[0].set(xlim=(0, 103), ylim=(0, 103), xlabel="Weighted recall (%)", ylabel="Weighted precision (%)",
                title="A  OrthoBench: development-exposed")
    axes[0].set_aspect("equal", adjustable="box")
    for i, r in enumerate(data):
        axes[1].plot(100 * r["three_kingdoms_f1"], i, marker=r["marker"], color=r["color"], markersize=8)
        axes[1].annotate(f"{100 * r['three_kingdoms_f1']:.2f}", (100 * r["three_kingdoms_f1"], i),
                         xytext=(7, 0), textcoords="offset points", va="center", fontsize=9)
    axes[1].set(yticks=range(len(data)), yticklabels=[r["label"] for r in data], xlim=(0, 111),
                xticks=(0, 25, 50, 75, 100), xlabel="BUSCO-reference pair F1 (%)",
                title="B  Three Kingdoms: supplementary")
    axes[1].invert_yaxis()
    for ax in axes:
        style(ax)
    fig.suptitle("Audited historical accuracy | Work in progress", fontsize=14)
    fig.text(0.5, 0.14, "BUSCO-only scoring excludes false positives involving non-reference genes.", ha="center", fontsize=10)
    legend(fig, data)
    fig.subplots_adjust(left=0.07, right=0.96, bottom=0.25, top=0.86, wspace=0.85)
    return fig


def qfo_endpoints(data):
    fig, axes = plt.subplots(2, 3, figsize=(12, 8.5))
    for ax, (metric, axis_names) in zip(axes.flat, QFO_AXES.items()):
        for r in data:
            if metric in r["qfo"]:
                point = r["qfo"][metric]
                ax.scatter(point["x"], point["y"], color=r["color"], marker=r["marker"], s=65,
                           edgecolors="white", linewidths=0.4)
        if axis_names[0] == "NR_ORTHOLOGS":
            ax.set_xscale("symlog", linthresh=1)
            ax.set_xlabel("Assessed relations (log scale)")
            ax.set_ylabel("FAS similarity" if metric == "FAS" else "Mean Schlicker similarity")
        else:
            ax.set(xlim=(-0.03, 1.03), xlabel="Recall (TPR)", ylabel="Precision (PPV)")
        ax.set(ylim=(-0.03, 1.03), title=metric.removesuffix(" F"))
        style(ax)
    missing = ", ".join(r["label"] for r in data if not r["qfo"])
    fig.suptitle("QfO individual endpoints | Development-exposed", fontsize=14)
    note = "Error bars omitted: recorded uncertainty fields differ in definition across endpoints."
    if missing:
        note += "\nPending, not plotted: " + missing
    fig.text(0.5, 0.145, note, ha="center", fontsize=9)
    legend(fig, data)
    fig.subplots_adjust(left=0.08, right=0.97, bottom=0.25, top=0.9, hspace=0.45, wspace=0.4)
    return fig


def paired_differences(data):
    selected = [r for r in data if r["uncertainty"] is not None]
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.5), sharey=True)
    for ax, metric, label in zip(axes, ("f_score", "precision", "recall"), ("F1", "Precision", "Recall")):
        for i, r in enumerate(selected):
            item = r["uncertainty"][metric]
            ax.hlines(i, *item["bonferroni_percentile_ci"], color=r["color"], linewidth=1.2)
            ax.hlines(i, *item["paired_percentile_ci"], color=r["color"], linewidth=4)
            ax.plot(item["difference_percentage_points"], i, marker=r["marker"], color=r["color"],
                    markersize=8, markeredgecolor="white", markeredgewidth=0.7)
        ax.axvline(0, color="#444444", linestyle="--", linewidth=1)
        ax.set(title=label, xlabel="Difference (percentage points)", ylim=(-0.65, 1.65),
               yticks=range(len(selected)), yticklabels=[r["label"] for r in selected])
        style(ax)
    axes[0].invert_yaxis()
    fig.suptitle("OrthoBench paired differences versus OrthoFinder 3.1.5 full", fontsize=13)
    fig.text(0.5, 0.09, "Thick: nominal 95% interval. Thin: Bonferroni-adjusted across six contrasts/metrics.\n"
             "20,000 paired RefOG bootstrap replicates; development-exposed, not selection-adjusted.", ha="center", fontsize=10)
    fig.subplots_adjust(left=0.19, right=0.97, bottom=0.3, top=0.8, wspace=0.35)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("Refusing to overwrite existing figure bundle")
    data = plotting_data(json.loads(args.comparison.read_text()))
    args.output.mkdir(parents=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "pdf.fonttype": 42, "svg.hashsalt": "orthohmm-publication-v1"})
    outputs = []
    for name, builder in (("accuracy_overview", overview), ("qfo_endpoints", qfo_endpoints),
                          ("orthobench_paired_differences", paired_differences)):
        fig = builder(data)
        for suffix in ("png", "pdf", "svg"):
            path = args.output / f"{name}.{suffix}"
            metadata = {"CreationDate": None, "ModDate": None} if suffix == "pdf" else {"Date": None} if suffix == "svg" else None
            fig.savefig(path, dpi=200, facecolor="white", metadata=metadata)
            if suffix == "svg":
                path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines()) + "\n")
            outputs.append(file_provenance(path))
        plt.close(fig)
    payload = {"schema_version": 1, "publication_ready": False,
               "generated_at": datetime.now(timezone.utc).isoformat(),
               "source": file_provenance(Path(__file__)), "input": file_provenance(args.comparison),
               "command": [sys.executable, *sys.argv], "matplotlib_version": matplotlib.__version__,
               "plotted_data": data, "outputs": outputs}
    (args.output / "manifest.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
