"""Integrate checked native-cell comparator intervals without altering frozen drafts."""

import argparse
import csv
import io
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


PINS = {
    "result": ("native_qfo_comparator_uncertainty_20261009_v1/report.json",
               "aea49a83bcce9a80237fd1c8f04ada43efeb9bd41f6e17ba0aa0c7f911018c02"),
    "execution": ("native_qfo_comparator_execution_20261009_v1.json",
                  "b7a407f57c7952079a129b8814ce96fcabb53456c88ceb3d719e01ce4eafa3b6"),
    "parent": ("native_qfo_terminal_failures_20261009_v1_manuscript.md",
               "16cb869bb9fc67905de0e98436242d7140851e26c10ada91f1670f169b871018"),
}
CELLS = tuple(f"p{p}_c{c}_r{r}" for p in (0, 1) for c in (0, 1) for r in (0, 1))
OBSERVED = ("p0_c0_r0", "p0_c0_r1", "p0_c1_r0", "p1_c0_r1")
COMPARATORS = ("orthofinder_3_1_5_full", "orthofinder_3_1_5_sequence_only")
METRICS = ("F1", "PPV", "TPR")
METHOD_ANCHOR = "### Fixed-Stratum SwissTrees Error Analyses\n"
RESULT_ANCHOR = "### Terminal Native QfO Outcomes\n"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def endpoints(result):
    require(result["schema"] == "native_qfo_comparator_sensitivity_v1"
            and result["planned_endpoints"] == result["multiplicity_endpoints"] == 48
            and result["estimated_endpoints"] == 24
            and result["replicates"] == result["new_bootstrap_draws"] == 100000
            and result["seed"] == 20260920 and result["alpha"] == .05
            and result["quantile_method"] == "linear" and result["protocol_controls_match"] is True,
            "Changed comparator statistic/protocol scope")
    require(len(result["families"]) == len(set(result["families"])) == 18
            and result["observed_cells"] == list(OBSERVED)
            and all(result[key] is False for key in (
                "publication_ready", "independent_confirmation", "unobserved_cells_imputed",
                "retained_intervals_reused", "new_accuracy_or_resource_admission")), "Changed evidence scope")
    expected = [(cell, comparator) for cell in CELLS for comparator in COMPARATORS]
    require([(row["candidate"], row["reference"]) for row in result["comparisons"]] == expected,
            "Changed planned contrast inventory/order")
    rows = []
    for comparison in result["comparisons"]:
        available = comparison["candidate"] in OBSERVED
        if available:
            require(comparison["status"] == "conditional_native_comparison_estimated"
                    and comparison["missing_reason"] is None
                    and set(comparison["metrics"]) == set(METRICS)
                    and [r["family"] for r in comparison["family_differences"]] == result["families"],
                    "Changed estimated contrast scope")
        else:
            require(comparison["status"] == "native_cell_unavailable"
                    and comparison["metrics"] is comparison["family_differences"] is None
                    and isinstance(comparison["missing_reason"], str) and comparison["missing_reason"],
                    "Missing endpoint was imputed")
        for metric in METRICS:
            value = comparison["metrics"][metric] if available else None
            row = dict(candidate=comparison["candidate"], reference=comparison["reference"],
                metric=metric, status=comparison["status"], missing_reason=comparison["missing_reason"],
                **dict.fromkeys(("difference", "nominal_low", "nominal_high", "adjusted_low",
                                 "adjusted_high", "family_wins", "family_ties", "family_losses")))
            if available:
                intervals = [value[key] for key in ("paired_percentile_ci", "bonferroni_percentile_ci")]
                require(type(value["difference"]) in (int, float) and math.isfinite(value["difference"])
                        and -1 <= value["difference"] <= 1
                        and all(len(pair) == 2 and all(type(v) in (int, float) and math.isfinite(v)
                                    and -1 <= v <= 1 for v in pair) and pair[0] <= pair[1] for pair in intervals)
                        and intervals[1][0] <= intervals[0][0] <= intervals[0][1] <= intervals[1][1],
                        "Invalid plotted interval")
                outcomes = [value[k] for k in ("family_wins", "family_ties", "family_losses")]
                require(all(type(v) is int and v >= 0 for v in outcomes) and sum(outcomes) == 18,
                        "Wrong descriptive family outcomes")
                row.update(difference=value["difference"], nominal_low=intervals[0][0], nominal_high=intervals[0][1],
                    adjusted_low=intervals[1][0], adjusted_high=intervals[1][1],
                    **dict(zip(("family_wins", "family_ties", "family_losses"), outcomes)))
            rows.append(row)
    return rows


def sections(result, figure_directory):
    rows = endpoints(result)
    methods = """### Retrospective Native-Cell Comparator Protocol

A separately frozen retrospective analysis compares all eight planned native
cells with full OrthoFinder 3.1.5 and its sequence-only MCL checkpoint, for
SwissTrees F1, precision and recall. Four cells are admitted, leaving 24 of the
48 endpoints estimable; missing endpoints remain in the correction. The 18
common families contain 563 disjoint represented proteins and 10,765 reference
relations. Raw counts use the official halving/unit-prior conversion, equivalent
to PPV=(TP+2)/(TP+FP+4) and TPR=(TP+2)/(TP+FN+4). Every replicate recomputes
harmonic F1 from macro-family precision and recall, not mean family F1.

The calculation uses 100,000 shared multinomial family draws, PCG64 seed
20260920, nominal percentile quantiles 0.025/0.975 and separate 48-endpoint
Bonferroni quantiles 0.05/96 and 1-0.05/96 with linear interpolation. This is
new computation, not interval transplantation from the earlier 24-endpoint
comparator or 42-endpoint internal factorial analysis. Both remain unchanged.
Independent literal-count readback checked statistics, all 24 intervals,
family outcomes and all 48 output rows; deterministic replay is verification,
not another dataset. Development exposure, family exchangeability, merged-
prediction dependence and finite tail resolution limit this conditional
analysis; it is not selection-adjusted or independent confirmation.
[Frozen protocol](NATIVE_QFO_COMPARATOR_UNCERTAINTY_PROTOCOL_20261009.md),
[actual execution/readback](native_qfo_comparator_execution_20261009_v1.json).

"""
    lines = ["### Retrospective Native-Cell Comparator Sensitivity", "",
        "Native minus OrthoFinder SwissTrees F1 differences and adjusted intervals",
        "below use raw 0-to-1 units. Family outcomes are descriptive, not independent",
        "pair counts. All planned missing comparisons remain unavailable.", "",
        "| Native Cell | Comparator | F1 Difference | 48-Endpoint Adjusted Interval | Family Wins/Ties/Losses |",
        "| --- | --- | ---: | --- | --- |"]
    for comparator, label in zip(COMPARATORS, ("Full OrthoFinder", "Sequence-Only OrthoFinder")):
        for cell in OBSERVED:
            row = next(r for r in rows if (r["candidate"], r["reference"], r["metric"]) == (cell, comparator, "F1"))
            outcomes = "/".join(str(row[k]) for k in ("family_wins", "family_ties", "family_losses"))
            lines.append(f"| {cell.upper().replace('_', '/')} | {label} | {row['difference']:+.6f} | "
                f"[{row['adjusted_low']:.6f}, {row['adjusted_high']:.6f}] | {outcomes} |")
    lines.extend(["", "All four F1 point estimates are below full OrthoFinder. Only P0/C0/R0",
        "has an adjusted F1 interval excluding zero in that direction. Both",
        "unreconciled cells have lower precision with adjusted intervals below zero.",
        "Reconciled precision differences versus full OrthoFinder have intervals",
        "including zero. Relative to its sequence-only checkpoint, P0/C0/R1 and",
        "P1/C0/R1 have positive precision intervals and negative recall intervals;",
        "their positive F1 point differences have intervals including zero.",
        "P0/C0/R0 also has lower recall than this checkpoint. Eight of 24 estimated",
        "adjusted intervals exclude zero, conditionally within this separate family.",
        "This is not a whole-study significance count, equivalence test or general",
        "superiority finding. No default selection or initial-HMM causal effect follows.", "",
        f"**Native Comparator Figure.** The [six-panel figure]({figure_directory}/native_comparator_sensitivity.pdf)",
        "shows F1, precision and recall against full OrthoFinder (top) and its",
        "sequence-only checkpoint (bottom). Dots: point differences. Thick lines:",
        "nominal 95% intervals. Thin lines: 48-endpoint adjusted intervals.",
        "The figure alone uses percentage points; the table/TSV retain raw units.",
        "All eight cells appear in both rows; unavailable points are not zeros.",
        "[Complete 48-row table](native_qfo_comparator_uncertainty_20261009_v1/TABLE.md),",
        "[full-precision machine result](native_qfo_comparator_uncertainty_20261009_v1/report.json),",
        "[claim-to-evidence addendum](NATIVE_QFO_COMPARATOR_RESULT_20261009.md).", "",
        "P1/C0/R0 is absent from this supplied fresh-cell snapshot, not globally",
        "missing: the historical high-sensitivity method remains separate.",
        "The other absent cells failed before inference, during assessment or",
        "during native inference, respectively; no score is imputed. The recovered",
        "P0/C0/R1 timing remains failed/ineligible. R changes pair-output semantics",
        "as well as reconciliation; initial HMM search remains on in every cell.",
        "These intervals do not cover other QfO endpoints, the secondary mean,",
        "resource timings or missing outcomes. The terminal failure section below",
        "retains its earlier scope; its no-new-draw statement concerns that failure",
        "reporting update, not this separately specified retrospective calculation.", "",
        "The [terminal resource consolidation](NATIVE_FACTORIAL_TERMINAL_RESOURCE_RESULT_20261009.md)",
        "joins all 13 observed ablation attempts with explicit accuracy/failure and",
        "resource-scope flags. It does not repair failed measurements or establish",
        "isolated performance. Original shared-host contention limits remain.", "", ""])
    return methods, "\n".join(lines)


def manuscript(parent, result, figure_directory):
    require(parent.count(METHOD_ANCHOR) == parent.count(RESULT_ANCHOR) == 1, "Ambiguous manuscript anchors")
    methods, results = sections(result, figure_directory)
    require(methods not in parent and results not in parent, "Already integrated manuscript")
    revised = parent.replace(METHOD_ANCHOR, methods + METHOD_ANCHOR).replace(RESULT_ANCHOR, results + RESULT_ANCHOR)
    require(revised.replace(methods, "", 1).replace(results, "", 1) == parent, "Changed parent body")
    return revised, methods, results


def plot(rows, output):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
    fig, axes = plt.subplots(2, 3, figsize=(16, 11), sharey=True)
    labels = [cell.upper().replace("_", "/") for cell in CELLS]
    missing = {"p0_c1_r1": "Pre-native failure", "p1_c0_r0": "Not supplied",
               "p1_c1_r0": "Assessment OOM", "p1_c1_r1": "Inference SIGSEGV"}
    for i, comparator in enumerate(COMPARATORS):
        for j, (metric, title, color) in enumerate(zip(METRICS, ("F1", "Precision", "Recall"),
                                                    ("#167d8d", "#bd4148", "#9a710c"))):
            axis = axes[i, j]
            axis.axvline(0, color="#777777", linewidth=.8, linestyle="--")
            for y, cell in enumerate(CELLS):
                row = next(r for r in rows if (r["candidate"], r["reference"], r["metric"]) == (cell, comparator, metric))
                if row["difference"] is None:
                    axis.text(0, y, missing[cell], ha="center", va="center", fontsize=9, color="#666666",
                              bbox={"facecolor": "white", "edgecolor": "none", "pad": 2})
                    continue
                axis.plot([100*row[k] for k in ("adjusted_low", "adjusted_high")], [y,y], color=color, linewidth=1.2)
                axis.plot([100*row[k] for k in ("nominal_low", "nominal_high")], [y,y], color=color, linewidth=4,
                          solid_capstyle="butt")
                axis.plot(100*row["difference"], y, "o", color=color, markersize=5)
            axis.set_title(title + (" vs full OrthoFinder" if i == 0 else " vs sequence-only checkpoint"), fontsize=11)
            axis.set_xlim(-60, 60)
            axis.set_xticks((-60, -30, 0, 30, 60))
            axis.set_ylim(7.6, -.6)
            axis.set_yticks(range(8), labels)
            axis.set_xlabel("Difference (percentage points)")
            axis.grid(axis="y", color="#eeeeee", linewidth=.7)
            axis.spines[["top", "right", "left"]].set_visible(False)
            axis.tick_params(axis="y", length=0)
    fig.suptitle("Native QfO: SwissTrees comparator sensitivity", x=.035, ha="left", fontsize=16)
    fig.text(.035, .935, "18 families; 100,000 shared draws. Thick: nominal 95% CI. Thin: 48-endpoint adjusted CI.")
    fig.text(.035, .025, "Retrospective, development-exposed. Initial HMM search retained. Missing endpoints are not imputed. No equivalence or default selection.", fontsize=9)
    fig.subplots_adjust(left=.115, right=.985, top=.885, bottom=.09, wspace=.16, hspace=.30)
    for suffix in ("png", "pdf"):
        fig.savefig(output / ("native_comparator_sensitivity." + suffix), dpi=180)
    plt.close(fig)


def run(root, output, revised_path):
    root, output, revised_path = Path(root).resolve(), Path(output).absolute(), Path(revised_path).absolute()
    directory = root / "benchmark_tools/results"
    for path in (output, revised_path):
        if path.exists() or path.is_symlink():
            raise FileExistsError(path)
    require(output.parent.resolve() == revised_path.parent.resolve() == directory.resolve(),
            "Generate beside the parent manuscript to preserve local links")
    docs, inputs = {}, []
    for key, (name, digest) in PINS.items():
        path = directory / name
        ref = record(path)
        require(ref["sha256"] == digest, "Changed scientific integration input: " + key)
        inputs.append(ref)
        docs[key] = path.read_text() if key == "parent" else json.loads(path.read_text())
    receipt = docs["execution"]
    require(receipt["schema"] == "native_qfo_comparator_execution_and_independent_readback_v1"
            and receipt["production"]["tool_result"]["exit_code"] == receipt["readback"]["tool_result"]["exit_code"] == 0
            and receipt["readback"]["result"]["status"] == "independent_literal_count_and_interval_readback_passed"
            and receipt["publication_ready"] is False
            and receipt["readback"]["result"]["report"] == {k: inputs[0][k] for k in ("bytes", "sha256")},
            "Unreviewed/mixed comparator result")
    rows = endpoints(docs["result"])
    text, methods, results = manuscript(docs["parent"], docs["result"], output.name)
    source = record(__file__)
    inputs.append(source)
    for ref in inputs:
        check(ref)
    output.mkdir(parents=True)
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    (output / "endpoints.tsv").write_text(stream.getvalue())
    (output / "methods_section.md").write_text(methods)
    (output / "results_section.md").write_text(results)
    plot(rows, output)
    revised_path.write_text(text)
    for ref in inputs:
        check(ref)
    manifest = dict(schema="native_qfo_comparator_integration_v1", inputs=inputs,
        outputs=[record(p) for p in sorted(output.iterdir())], manuscript=record(revised_path),
        planned_endpoints=48, estimated_endpoints=24, unavailable_endpoints=24,
        parent_body_unchanged_except_insertions=True, plot_scale=100, table_scale=1,
        no_new_bootstrap=True, no_inference_or_scoring=True, publication_ready=False,
        visual_reviewed=False, manuscript_rendered=False)
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "output", "manuscript"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.root, args.output, args.manuscript), allow_nan=False))
