"""Plot a self-contained generation-timing snapshot without inference imports."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from statistics import median

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


METHODS = ("original", "sensitive", "phylogeny")
LABELS = ("OrthoHMM1 (v0.1.0)", "OrthoHMM2 sensitive", "OrthoHMM2 + phylogeny")
COLORS = ("#666666", "#0072B2", "#009E73")


def summarize(data):
    runs = data["runs"]
    expected = [(3*r+j, r+1, METHODS[(r+j) % 3]) for r in range(3) for j in range(3)]
    if [(r["index"],r["repeat"],r["method"]) for r in runs] != expected:
        raise ValueError("Require every planned attempt in the fixed balanced order")
    if any(r["status"] != "validated" for r in runs):
        raise ValueError("Do not graph failed or censored runs as successful runtimes")
    summaries = []
    for method, label in zip(METHODS, LABELS):
        selected = [r for r in runs if r["method"] == method]
        row = dict(method=method, label=label)
        for metric in ("elapsed_s", "search_s"):
            values = [float(r[metric]) for r in selected]
            if any(not math.isfinite(v) or v <= 0 for v in values):
                raise ValueError("Timings must be finite and positive")
            row[metric] = dict(values=values, median=median(values), minimum=min(values), maximum=max(values))
        summaries.append(row)
    for row in summaries:
        row["original_to_modern_elapsed_ratio"] = summaries[0]["elapsed_s"]["median"] / row["elapsed_s"]["median"]
    return summaries


def export(snapshot, output):
    if output.exists():
        raise FileExistsError("Preserve previous figure attempts; use a new directory")
    data = json.loads(snapshot.read_text())
    summaries = summarize(data)
    if data["summary"] != summaries:
        raise ValueError("Supplied medians or ratios differ from raw measurements")
    output.mkdir(parents=True)
    plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":9,
                         "svg.fonttype":"none", "pdf.fonttype":42})
    figure, axes = plt.subplots(1, 2, figsize=(11.3, 4.7))
    figure.subplots_adjust(left=.225, right=.97, top=.83, bottom=.25, wspace=.32)
    for axis, metric, title in zip(axes, ("elapsed_s", "search_s"),
            ("A  End-to-end group inference", "B  Initial homology search")):
        maximum = max(row[metric]["maximum"] for row in summaries)
        for y, (row, color) in enumerate(zip(summaries, COLORS)):
            values = row[metric]
            axis.hlines(y, values["minimum"], values["maximum"], color=color, linewidth=3, alpha=.45)
            axis.scatter(values["values"], [y-.10,y,y+.10], color=color, s=25, zorder=3)
            axis.plot(values["median"], y, marker="|", markersize=17, color="#161616", zorder=4)
            axis.text(values["maximum"] + maximum*.035, y, f'{values["median"]:.1f} s', va="center", fontsize=9)
        axis.set(xlim=(0, maximum*1.30), ylim=(2.5,-.5), yticks=range(3),
                 yticklabels=LABELS if metric == "elapsed_s" else ["", "", ""],
                 xlabel="Elapsed time (seconds)", title=title)
        axis.spines[["right", "top", "left"]].set_visible(False)
        axis.tick_params(axis="y", length=0, pad=10)
        axis.grid(axis="x", color="#e7e7e7", linewidth=.6)
        axis.set_axisbelow(True)
    proteins = sum(r["proteins"] for r in data["input"])
    figure.suptitle(f'{len(data["input"])} bacterial proteomes | {proteins:,} proteins | 16 eligible physical cores', fontsize=11, y=.96)
    ratios = summaries[1:]
    figure.text(.225,.12,"Original/modern elapsed-median ratios: " + "; ".join(
        f'{"sensitive" if r["method"] == "sensitive" else "phylogeny"} {r["original_to_modern_elapsed_ratio"]:.2f}x' for r in ratios), fontsize=9)
    figure.text(.225,.065,"Three fresh technical repeats; ticks are medians, lines observed ranges. Shared host; accuracy not assessed.", fontsize=8, color="#555555")
    try:
        for extension in ("png", "pdf", "svg"):
            figure.savefig(output / ("fig6_generations."+extension), dpi=220)
    finally:
        plt.close(figure)
    (output / "timing.json").write_bytes(snapshot.read_bytes())
    with (output / "runs.tsv").open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=("index", "repeat", "method", "status", "elapsed_s", "search_s"), delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(data["runs"])
    records = []
    for p in sorted(output.iterdir()):
        content = p.read_bytes()
        records.append(dict(path=p.name, bytes=len(content), sha256=hashlib.sha256(content).hexdigest()))
    (output / "manifest.json").write_text(json.dumps(dict(outputs=records,
        source_snapshot_sha256=hashlib.sha256(snapshot.read_bytes()).hexdigest(),
        scope="Original descriptive shared-host figure; not an isolated speedup or an equal-accuracy comparison."), indent=2, sort_keys=True)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    export(args.snapshot, args.output)
