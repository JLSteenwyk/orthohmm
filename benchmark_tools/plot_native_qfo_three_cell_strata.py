"""Render verified fixed-stratum contrasts as descriptive publication panels."""

import argparse
import csv
import json
import math
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.export_native_qfo_three_cell_strata import record, require

PINS = {
    "report": ("native_qfo_three_cell_strata_20261007_v1/report.json", "55ed72ce771c756424877cfa412da329ff104efb7dca4b35113b0d80703a59e5"),
    "readback": ("native_qfo_three_cell_strata_readback_20261007_v2.json", "6f3eef63c0e96e1d85805fcfb70ed9b66fbb1db76a08ca2a105483422a0f6969"),
}
SELECTION = (
    ("sequence", "all", "All reference families"),
    ("sequence", "higher_entropy", "Composition: higher entropy"),
    ("sequence", "lower_entropy", "Composition: lower entropy"),
    ("sequence", "short_relative", "Length: short-relative"),
    ("sequence", "not_short_relative", "Length: not-short-relative"),
    ("domain", "median_pfam_types_at_least_two", "Pfam types: median >=2"),
    ("domain", "median_pfam_types_below_two", "Pfam types: median <2"),
    ("domain", "repeated_type_fraction_at_least_quarter", "Repeated Pfam fraction: >=1/4"),
    ("domain", "repeated_type_fraction_below_quarter", "Repeated Pfam fraction: <1/4"),
    ("duplication", "lower_duplication_fraction", "Reference duplication fraction: lower"),
    ("duplication", "upper_duplication_fraction", "Reference duplication fraction: upper"),
)
CONTRASTS = ("R_at_P0_C0", "C_at_P0_R0")
METRICS = ("F1", "PPV", "TPR")
COLORS = ("#202020", "#b44735", "#12786f")
MARKERS = ("o", "s", "^")
FIELDS = ("contrast", "suite", "stratum", "families", "label", "metric", "difference_pp")
CAPTION = (
    "Dots are descriptive changes, not confidence intervals. Axes use different ranges, both in percentage points.\n"
    "All 18 families are development-exposed; strata overlap. Empty bins and redundant all-family copies remain in the full tables.\n"
    "Length/composition, Pfam annotation and reference duplication descriptors are not causal biological validation."
)


def figure_data(report, readback):
    require(report["schema"] == "native_qfo_three_cell_strata_v1"
            and readback["schema"] == "native_qfo_three_cell_strata_rational_readback_v2"
            and readback["score_rows_checked"] == 60 and readback["differences_checked"] == 40
            and readback["family_rows_checked"] == 54 and readback["proteins_checked"] == 563,
            "Incomplete verified stratum input")
    for document in (report, readback):
        require(all(document[k] is False for k in ("publication_ready", "new_uncertainty",
                "new_accuracy_or_resource_admission", "independent_confirmation", "scientific_timings_admitted"))
                and document["new_bootstrap_draws"] == 0, "Inflated plotted scope")
    lookup = {(r["suite"], r["stratum"], r["contrast"]): r for r in report["differences"]}
    require(len(lookup) == len(report["differences"]) == 40, "Duplicate/missing contrast row")
    represented = {(s, n) for s, n, _ in SELECTION}
    excluded = []
    for suite, bins in report["bins"].items():
        for name, members in bins.items():
            if (suite, name) not in represented:
                reason = "empty_bin" if not members else "redundant_all_family_copy"
                require(not members or members == report["bins"]["sequence"]["all"],
                        "Nonempty nonredundant bin excluded")
                excluded.append(dict(suite=suite, stratum=name, reason=reason, families=len(members)))
    points = []
    for contrast in CONTRASTS:
        candidate = "p0_c0_r1" if contrast == CONTRASTS[0] else "p0_c1_r0"
        for suite, name, label in SELECTION:
            row = lookup[suite, name, contrast]
            members = report["bins"][suite][name]
            require(row["status"] == "descriptive" and row["families"] == len(members) > 0
                    and row["family_members"] == members and row["reference"] == "p0_c0_r0"
                    and row["candidate"] == candidate, "Changed plotted unit/contrast")
            for metric in METRICS:
                value = row[metric]
                require(type(value) in (int, float) and math.isfinite(value) and -1 <= value <= 1,
                        "Invalid plotted difference")
                points.append(dict(contrast=contrast, suite=suite, stratum=name, families=len(members),
                                   label=f"{label} (n={len(members)})", metric=metric, difference_pp=100 * value))
    require(len(points) == 66 and len(excluded) == 9, "Incomplete plotted/excluded inventory")
    return points, excluded


def render(points, output):
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 10,
                         "svg.fonttype": "none", "svg.hashsalt": "native-three-cell-strata", "pdf.fonttype": 42}):
        figure, axes = plt.subplots(1, 2, figsize=(14, 8), sharey=True)
        labels = []
        for axis, contrast, limits, title in zip(axes, CONTRASTS, ((-15, 47), (-6, 8)),
            ("A  Reconciliation minus baseline", "B  Candidate expansion minus baseline")):
            selected = [p for p in points if p["contrast"] == contrast]
            axis.axvline(0, color="#777777", linestyle="--", linewidth=.8)
            for i, metric in enumerate(METRICS):
                series = [p for p in selected if p["metric"] == metric]
                values = [p["difference_pp"] for p in series]
                require(all(limits[0] < v < limits[1] for v in values), "Clipped plotted point")
                axis.scatter(values, [y + (i - 1) * .17 for y in range(11)],
                             s=36, color=COLORS[i], marker=MARKERS[i], edgecolor="white", linewidth=.45, zorder=3)
                if not labels:
                    labels = [p["label"] for p in series]
            axis.set(xlim=limits, ylim=(10.6, -.65), xlabel="Difference (percentage points)")
            axis.set_yticks(range(11))
            axis.set_title(title, loc="left", fontsize=11, pad=12)
            axis.grid(axis="x", color="#e5e5e5", linewidth=.7)
            for y in (.5, 4.5, 8.5):
                axis.axhline(y, color="#dddddd", linewidth=.7)
            for spine in ("top", "right"):
                axis.spines[spine].set_visible(False)
        axes[0].set_yticklabels(labels)
        figure.suptitle("Native SwissTrees fixed-stratum contrasts", y=.98, fontsize=15)
        figure.text(.5, .938, "Initial HMM search on; downstream profile refinement off in all three cells",
                    ha="center", fontsize=10)
        legend = [Line2D([], [], marker=m, color=c, linestyle="none", markersize=6, label=label)
                  for m, c, label in zip(MARKERS, COLORS, ("F1", "Precision", "Recall"))]
        figure.legend(handles=legend, loc="upper center", bbox_to_anchor=(.68, .91), frameon=False, ncol=3)
        figure.text(.02, .055, CAPTION, va="bottom", fontsize=9, linespacing=1.5)
        figure.subplots_adjust(left=.30, right=.985, top=.83, bottom=.18, wspace=.16)
        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        bounds = figure.bbox
        for text in figure.findobj(match=matplotlib.text.Text):
            if text.get_visible() and text.get_text():
                box = text.get_window_extent(renderer)
                require(box.x0 >= bounds.x0 and box.y0 >= bounds.y0 and box.x1 <= bounds.x1
                        and box.y1 <= bounds.y1, "Figure text exceeds canvas")
        for suffix in ("png", "pdf", "svg"):
            figure.savefig(output / ("fixed_stratum_contrasts." + suffix), dpi=160,
                           facecolor="white", metadata={"Creator": "OrthoHMM benchmark workflow"} if suffix != "png" else None)
        plt.close(figure)


def build(repo, output):
    output = Path(output)
    require(not output.exists() and not output.is_symlink(), "Output already exists")
    refs, docs = {}, {}
    for key, (name, sha) in PINS.items():
        ref = record(Path(repo) / "benchmark_tools/results" / name)
        require(ref["sha256"] == sha, "Changed verified figure input")
        refs[key], docs[key] = ref, json.loads(Path(ref["path"]).read_text())
    require(docs["readback"]["report"] == refs["report"], "Wrong readback/report binding")
    checked = [*refs.values(), docs["report"]["source"], docs["readback"]["source"], record(__file__)]
    for ref in checked:
        require(record(ref["path"]) == ref, "Changed figure source")
    points, excluded = figure_data(docs["report"], docs["readback"])
    output.mkdir(parents=True)
    with (output / "plotted_rows.tsv").open("x", newline="") as stream:
        writer = csv.DictWriter(stream, FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(points)
    render(points, output)
    for ref in checked:
        require(record(ref["path"]) == ref, "Figure evidence changed during render")
    manifest = dict(schema="native_qfo_three_cell_strata_figure_v1", source=record(__file__), inputs=refs,
        checked_inputs=checked, plotted_rows=points, excluded_bins=excluded, caption=CAPTION,
        matplotlib_version=matplotlib.__version__, figure_inches=[14, 8], png_dpi=160,
        x_limits_pp=[[-15, 47], [-6, 8]], new_uncertainty=False, new_bootstrap_draws=0,
        independent_confirmation=False, new_accuracy_or_resource_admission=False, publication_ready=False,
        outputs=[record(output / name) for name in ("plotted_rows.tsv", "fixed_stratum_contrasts.png",
                                                   "fixed_stratum_contrasts.pdf", "fixed_stratum_contrasts.svg")])
    with (output / "manifest.json").open("x") as stream:
        json.dump(manifest, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = build(args.repo, args.output)
    print(json.dumps(dict(plotted_points=len(result["plotted_rows"]), excluded_bins=len(result["excluded_bins"]))))
