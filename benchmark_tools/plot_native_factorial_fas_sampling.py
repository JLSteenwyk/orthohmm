"""Plot conditional expected-FAS ranges separately from descriptive native Z."""

from itertools import combinations
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from benchmark_tools.export_native_factorial_progress import require
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check

CELLS = ("p0_c0_r0", "p0_c0_r1", "p0_c1_r0", "p1_c0_r1")
REPORT_SHA = "785921f29f176d92417e9ef239c69dc896410e00ad764594d8ea7119b63d7be1"


def validate(report):
    require(report["schema"] == "native_factorial_conditional_fas_sampling_v1"
            and report["status"] == "conditional_design_ranges_complete"
            and report["target"] == "expected_native_post_attrition_ratio_under_fixed_design"
            and report["alpha"] == .05 and report["component_error"] == .05 / 12
            and report["components"] == 12 and report["joint_error_bound"] == .05,
            "Changed conditional sampling scope")
    for key in ("historical_scores_rerun", "observed_scores_changed", "biological_generalization_intervals",
                "unconditional_historical_interval_admission", "other_endpoint_uncertainty_admitted",
                "failed_factorial_cells_repaired", "publication_ready"):
        require(report[key] is False, "Changed scientific admission: " + key)
    methods, contrasts = report["methods"], report["contrasts"]
    require([r["cell"] for r in methods] == list(CELLS), "Changed cell inventory")
    for row in methods:
        bounds = row["interval"]["expected_native_mean_bounds"]
        z = row["observed_native_mean"]
        require(len(bounds) == 2 and all(math.isfinite(x) for x in [*bounds, z])
                and 0 <= bounds[0] <= bounds[1] <= 1 and 0 <= z <= 1,
                "Invalid plotted mean/range")
    pairs = list(combinations(methods, 2))
    require(len(contrasts) == 6, "Incomplete contrast inventory")
    for contrast, (left, right) in zip(contrasts, pairs):
        expected = [left["interval"]["expected_native_mean_bounds"][0] - right["interval"]["expected_native_mean_bounds"][1],
                    left["interval"]["expected_native_mean_bounds"][1] - right["interval"]["expected_native_mean_bounds"][0]]
        require((contrast["left"], contrast["right"]) == (left["cell"], right["cell"])
                and contrast["conditional_expected_difference_bounds"] == expected
                and contrast["observed_difference"] == left["observed_native_mean"] - right["observed_native_mean"]
                and contrast["zero_included"] is (expected[0] <= 0 <= expected[1]),
                "Changed plotted contrast arithmetic or zero status")
    return methods, contrasts


def label(cell):
    return cell.upper().replace("_", "/")


def range_and_observation(axis, bounds, observed, y, color):
    axis.hlines(y - .12, *bounds, color=color, linewidth=3)
    axis.vlines(bounds, y - .19, y - .05, color=color, linewidth=1)
    axis.plot(observed, y + .13, marker="o", linestyle="none", markersize=5,
              markerfacecolor="white", markeredgecolor="#333333", markeredgewidth=1)


def render(report, output):
    methods, contrasts = validate(report)
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 10,
                         "svg.hashsalt": "native-factorial-fas-sampling", "svg.fonttype": "none",
                         "pdf.fonttype": 42}):
        fig, axes = plt.subplots(1, 2, figsize=(14.5, 6.8))
        for y, row in enumerate(methods):
            range_and_observation(axes[0], row["interval"]["expected_native_mean_bounds"],
                                  row["observed_native_mean"], y, "#12786f")
        axes[0].set_yticks(range(4), [label(r["cell"]) for r in methods])
        axes[0].set(title="A  Expected repeated native FAS", xlabel="FAS score", ylim=(3.6, -.6))
        for y, row in enumerate(contrasts):
            range_and_observation(axes[1], row["conditional_expected_difference_bounds"],
                                  row["observed_difference"], y, "#12786f")
        axes[1].set_yticks(range(6), [label(r["left"]) + " minus " + label(r["right"]) for r in contrasts])
        axes[1].tick_params(axis="y", labelsize=8)
        axes[1].axvline(0, color="#777777", linestyle="--", linewidth=.8)
        axes[1].set(title="B  All six expected-ratio differences", xlabel="FAS difference", ylim=(5.6, -.6))
        for axis in axes:
            axis.spines[["top", "right"]].set_visible(False)
            axis.grid(axis="x", color="#e9e9e9", linewidth=.7)
            axis.tick_params(length=3)
            title = axis.get_title()
            axis.set_title("")
            axis.set_title(title, loc="left", fontsize=11, pad=12)
        fig.suptitle("Native FAS sampling uncertainty", x=.07, ha="left", fontsize=16, y=.97)
        fig.text(.07, .91, "Four QfO development-exposed ablations; initial HMM search remains on. P = refinement, C = expansion, R = reconciliation.", fontsize=10)
        handles = [plt.Line2D([], [], color="#12786f", linewidth=3),
                   plt.Line2D([], [], marker="o", linestyle="none", markerfacecolor="white", markeredgecolor="#333333")]
        fig.legend(handles, ("Simultaneous conditional 95% expected-ratio range", "Observed native Z (descriptive)"),
                   loc="upper left", bbox_to_anchor=(.063, .87), ncol=2, frameon=False, fontsize=9)
        fig.text(.07, .105, "Ranges account for unknown means in both strata and unknown numeric-return counts; all six differences use the same joint rectangle.", fontsize=9)
        fig.text(.07, .065, "Uniform without-replacement sampling and fixed pair outcomes are assumptions. These are not biological/generalization error bars.", fontsize=9)
        fig.text(.07, .025, "Circles are offset descriptive observations, not expected-ratio estimates or interval centers. Failed cells/timing remain unavailable.", fontsize=9)
        fig.subplots_adjust(left=.12, right=.98, top=.74, bottom=.22, wspace=.95)
        for suffix in ("png", "svg", "pdf"):
            fig.savefig(output / ("native_factorial_fas_sampling." + suffix), dpi=200,
                        metadata={"Creator": "OrthoHMM benchmark workflow"})
        plt.close(fig)


def run(repo, output):
    require(not output.exists() and not output.is_symlink(), "Output already exists")
    path = repo / "benchmark_tools/results/native_factorial_fas_sampling_20261010_v1/report.json"
    ref = record(path)
    require(ref["sha256"] == REPORT_SHA, "Changed retained conditional result")
    report = json.loads(path.read_text())
    validate(report)
    output.mkdir(parents=True)
    render(report, output)
    check(ref)
    result = dict(schema="native_factorial_fas_sampling_figure_v1", source=record(__file__), report=ref,
                  output_files=[record(output / ("native_factorial_fas_sampling." + s)) for s in ("png", "svg", "pdf")],
                  plotted_methods=report["methods"], plotted_contrasts=report["contrasts"],
                  new_scoring_or_sampling=False, biological_generalization_intervals=False,
                  publication_ready=False, matplotlib_version=matplotlib.__version__)
    with (output / "figure.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result
