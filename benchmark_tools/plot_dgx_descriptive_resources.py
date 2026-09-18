"""Plot every validated DGX repeat with explicit observational limitations."""

import argparse
import json
import math
from pathlib import Path
import statistics
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen

METHODS = (("orthohmm_high_sensitivity", "OrthoHMM high sensitivity", "#007d83", "o"),
           ("orthohmm_satellite_v2", "OrthoHMM inferred phylogeny", "#a66b0b", "s"),
           ("orthofinder_full", "OrthoFinder 3.1.5 full", "#755297", "^"))
METRICS = (("native_elapsed_seconds", 60., "A  Native elapsed time", "Minutes"),
           ("native_cpu_seconds", 3600., "B  Native CPU use", "CPU-hours"),
           ("cgroup_peak_bytes", 1024.**3, "C  Cgroup memory peak", "GiB (not process RSS)"))


def validate(report):
    if (report["status"] != "descriptive_resource_observations_not_controlled_comparison"
            or report["scientific_timings_admitted"] != 0 or report["publication_ready"] is not False
            or len(report["runs"]) != 27 or len(report["summaries"]) != 9):
        raise ValueError("Require complete descriptive, non-admitted panel")
    groups = {}
    indices = set()
    for row in report["runs"]:
        if (row["native_validation_status"] != "native_outputs_validated"
                or row["scientific_timing_admitted"] is not False or row["host_status"] != "inconclusive"
                or type(row["index"]) is not int or row["index"] not in range(27) or row["index"] in indices
                or type(row["repeat"]) is not int or row["repeat"] not in range(3)):
            raise ValueError("Invalid native run or changed admission/host status")
        for key, *_ in METRICS:
            if isinstance(row[key], bool) or not math.isfinite(row[key]) or row[key] <= 0:
                raise ValueError("Invalid resource value")
        if row["proteins"] != {4: 73266, 8: 165168, 12: 251378}.get(row["proteomes"]):
            raise ValueError("Input count differs from figure annotation")
        indices.add(row["index"])
        groups.setdefault((row["method"], row["proteomes"]), []).append(row)
    expected = {(method, size) for method, *_ in METHODS for size in (4, 8, 12)}
    if set(groups) != expected:
        raise ValueError("Changed method/size inventory")
    summaries = {(s["method"], s["proteomes"]): s for s in report["summaries"]}
    if set(summaries) != expected:
        raise ValueError("Missing or duplicate summary cells")
    for key, rows in groups.items():
        if len(rows) != 3 or {r["repeat"] for r in rows} != {0, 1, 2}:
            raise ValueError("Missing or duplicate repeat")
        for metric, *_ in METRICS:
            values = [r[metric] for r in rows]
            if summaries[key][metric] != {"median": statistics.median(values), "minimum": min(values), "maximum": max(values)}:
                raise ValueError("Summary differs from all three retained repeats")
    return groups


def plot(report):
    groups = validate(report)
    fig, axes = plt.subplots(1, 3, figsize=(15, 7.5))
    fig.subplots_adjust(left=.065, right=.985, bottom=.30, top=.76, wspace=.32)
    fig.suptitle("Dedicated DGX resource observations", x=.045, ha="left", y=.975, fontsize=19)
    fig.text(.045, .915, "27 validated native runs | 20 CPUs and 96 GiB per allocation | three repeats per method and size", fontsize=11)
    handles = []
    for ax, (metric, divisor, title, ylabel) in zip(axes, METRICS):
        for method_index, (method, label, color, marker) in enumerate(METHODS):
            for position, size in enumerate((4, 8, 12)):
                rows = sorted(groups[method, size], key=lambda r: r["repeat"])
                values = [r[metric] / divisor for r in rows]
                x = position + (method_index - 1) * .22
                ax.plot([x, x], [min(values), max(values)], color=color, linewidth=1.5, zorder=2)
                ax.plot([x - .045, x + .045], [statistics.median(values)] * 2, color=color, linewidth=2.5, zorder=3)
                points = ax.scatter([x - .035, x, x + .035], values, s=26, marker=marker,
                                    facecolors="none", edgecolors=color, linewidths=1.1, zorder=4, label=label)
                if ax is axes[0] and position == 0:
                    handles.append(points)
        ax.set_title(title, loc="left", fontsize=12, pad=12)
        ax.set_xticks(range(3), ["4", "8", "12"])
        ax.set_xlabel("Complete proteomes")
        ax.set_ylabel(ylabel)
        ax.set_xlim(-.5, 2.5)
        ax.set_ylim(bottom=0)
        ax.grid(axis="y", color="#dddddd", linewidth=.6)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)
    fig.legend(handles, [item[1] for item in METHODS], loc="upper left", bbox_to_anchor=(.04, .87),
               ncol=3, frameon=False, fontsize=10)
    notes = [
        "Open symbols: individual repeats. Thick ticks: median. Vertical bars: observed range, not confidence intervals.",
        "All host assessments remain inconclusive. Matched allocations do not establish controlled speed comparisons.",
        "Native elapsed/CPU values are from GNU time. Cgroup peak includes charged cache/kernel memory and may predate launch.",
        "Sampled aggregate RSS is not plotted: process-read gaps differ by method. No speedup ratios or complexity fit are implied.",
        "One nested input series (73,266 / 165,168 / 251,378 proteins); taxon composition and dataset size co-vary.",
    ]
    for y, note in zip((.225, .18, .135, .09, .045), notes):
        fig.text(.045, y, note, fontsize=10)
    return fig


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = read_frozen(args.results, args.sha256)
    source, plotter = record(args.results), record(__file__)
    fig = plot(report)
    args.output.mkdir(parents=True)
    outputs = []
    for extension in ("png", "pdf", "svg"):
        path = args.output / ("dgx_descriptive_resources." + extension)
        fig.savefig(path, dpi=180)
        outputs.append(record(path))
    plt.close(fig)
    check(source)
    check(plotter)
    manifest = {"source_results": source, "plotter": plotter, "outputs": outputs,
        "matplotlib": matplotlib.__version__, "runs_shown": 27, "panels": 3,
        "scientific_timings_admitted": 0, "publication_ready": False}
    with (args.output / "manifest.json").open("x") as stream:
        json.dump(manifest, stream, indent=2, sort_keys=True)
        stream.write("\n")
