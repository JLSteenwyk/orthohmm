"""Plot frozen descriptive duplication strata without adding inference."""

import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from benchmark_tools.plot_swiss_descriptive_features import validate, METHODS, METRICS
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

TABLE_SHA = "e6ec609604d9ab134bb7f7231a155d8d6685c05b44870affd840097752720619"
HELPER_SHA = "d9947a6679117dee52e46153a8b27f3d6bde6c218ffee20ed958d6b54afc7e5e"
BINS = dict(all=18, lower_duplication_fraction=9, upper_duplication_fraction=9,
            missing_duplication_fraction=0)
STRATA = (("lower_duplication_fraction", "Lower annotation fraction (9 families)", "#007f87", "o"),
          ("upper_duplication_fraction", "Higher annotation fraction (9 families)", "#b1394b", "s"))


def render(root, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    table = root / "benchmark_tools/results/swiss_duplication_strata_20260923/scores.tsv"
    helper = Path(__file__).with_name("plot_swiss_descriptive_features.py")
    inputs = [record(table), record(helper), record(__file__)]
    if [r["sha256"] for r in inputs[:2]] != [TABLE_SHA, HELPER_SHA]:
        raise ValueError("Changed pinned source table or validator")
    with table.open() as stream:
        indexed = validate(list(csv.DictReader(stream, delimiter="\t")), BINS)
    output.mkdir(parents=True)
    fig, axes = plt.subplots(1, 3, figsize=(15, 9))
    fig.subplots_adjust(left=.23, right=.98, top=.78, bottom=.22, wspace=.13)
    fig.suptitle("QfO SwissTrees: mapped-tree duplication annotations", x=.025, y=.96, ha="left", fontsize=18)
    fig.text(.025, .91, "Descriptive differences versus full OrthoFinder 3.1.5; no subgroup intervals or significance tests.", fontsize=11)
    handles = [Line2D([], [], color=color, marker=marker, linestyle="None", label=label)
               for _, label, color, marker in STRATA]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(.02, .88), ncol=2, frameon=False, fontsize=10)
    endpoints = []
    for column, (ax, metric) in enumerate(zip(axes, METRICS)):
        ax.axvline(0, color="#888888", linestyle="--", linewidth=1)
        for y, method in enumerate(METHODS):
            for i, (stratum, _, color, marker) in enumerate(STRATA):
                row = indexed[method, stratum]
                value = None if row["delta_"+metric] == "NA" else float(row["delta_"+metric])
                endpoints.append(dict(method=method, stratum=stratum, metric=metric,
                                      difference=value, status=row["status"]))
                if value is not None:
                    ax.plot(value*100, y+(i-.5)*.18, marker=marker, color=color, markersize=5, linestyle="None")
            if method == "orthomcl_1_4":
                ax.text(.5, y, "Not yet admitted", transform=ax.get_yaxis_transform(), ha="center", va="center", fontsize=9)
        ax.set(xlim=(-45, 20), ylim=(7.6, -.6), xlabel="Difference (percentage points)")
        ax.set_yticks(range(8), list(METHODS.values()) if column == 0 else [""]*8, fontsize=10)
        ax.set_title(f"{'ABC'[column]}  " + dict(F1="F1", PPV="Precision", TPR="Recall")[metric], loc="left", fontsize=12)
        ax.grid(axis="x", alpha=.15)
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
    notes = ["Development-exposed associations, not causal effects. F1 is harmonic macro precision/recall; all panels share the same x-axis scale.",
        "*Sequence-only OrthoFinder uses group-clique pairs. FastOMA uses a supplied tree. Configuration contrasts are not pure ablations.",
        "Feature: explicit duplication annotations / informative mapped nodes; median 7/48. Not an evolutionary duplication rate.",
        "Default-S nodes are not explicit speciation evidence. Missing bin: 0 families. Full-panel and empty-bin scores remain in the source table."]
    for y, note in zip((.155, .115, .075, .035), notes):
        fig.text(.025, y, note, fontsize=9)
    if any(r["difference"] is not None and not -45 < 100*r["difference"] < 20 for r in endpoints):
        plt.close(fig)
        raise ValueError("Endpoint outside figure limits")
    try:
        for extension in ("png", "pdf", "svg"):
            fig.savefig(output / ("swiss_duplication_descriptive."+extension), dpi=180)
    finally:
        plt.close(fig)
    with (output / "endpoints.tsv").open("x") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(endpoints[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({k: "NA" if v is None else v for k, v in row.items()} for row in endpoints)
    for item in inputs:
        check(item)
    report = dict(inputs=inputs, outputs=[record(p) for p in sorted(output.iterdir())],
        endpoints=len(endpoints), unavailable_endpoints=sum(r["difference"] is None for r in endpoints),
        matplotlib_version=matplotlib.__version__, new_inferential_claims=False,
        publication_ready=False, display_conversion="Raw TSV differences multiplied by 100 for display")
    (output / "manifest.json").write_text(json.dumps(report, indent=2, sort_keys=True)+"\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    render(args.root.resolve(), args.output.absolute())
