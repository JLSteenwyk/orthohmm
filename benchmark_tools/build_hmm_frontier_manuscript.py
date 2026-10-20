"""Export original manuscript figures from the completed HMM campaign.

The snapshot plotting mode needs only NumPy and Matplotlib, not raw proteomes
or the repository's historical scientific runners.
"""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from statistics import median
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Rectangle
import numpy as np


ARCHIVE = "benchmark_tools/results/hmm_campaign_final_20261020"
TIMING = "benchmark_tools/results/threadripper_shared_panel_snapshot_20261004_v27/panel.json"
PINS = {
    ARCHIVE + "/summary.json": "d6b6a0d43d661db6c5822b625abc576cc8cb1f094c6a304729f7b5d69eedd47d",
    ARCHIVE + "/scores.tsv": "4862010db93cf45cc6e9c1894e232d9cef709736ee6152b730a36fe1758eceb1",
    ARCHIVE + "/endpoint_diagnostics.tsv": "2be9fa4bfbd1589798261d8d7e00659729e8e400a77ff3efc14bc0d993df2468",
    TIMING: "1a586d5d85262267ddfbdd4fccc43a2d43208dae5b028368964d710d2b587187",
}
ENDPOINTS = {
    "OB": "OrthoBench: OrthoBench F1", "VGNC": "QfO: VGNC F",
    "SwissTrees": "QfO: SwissTrees F", "TreeFam": "QfO: TreeFam-A F",
    "EC": "QfO: EC", "GO": "QfO: GO", "FAS": "QfO: FAS",
    "TK": "ThreeKingdoms: Three Kingdoms F1",
}
PRIMARY = ("OB", "VGNC", "SwissTrees")
ALL_F1 = (*PRIMARY, "TreeFam", "TK")
METHODS = (
    ("profile_weighting:none", "OrthoHMM canonical", "OH-C", "#0072B2", "o"),
    ("profile_weighting:henikoff", "OrthoHMM Henikoff", "OH-W", "#117733", ">"),
    ("profile_weighting:henikoff_cluster", "OrthoHMM clustered Henikoff", "OH-N", "#44AA99", "<"),
    ("projection_recovery:projection_off", "OrthoHMM projection repair", "OH-R", "#882255", "h"),
    ("orthohmm_high_sensitivity", "OrthoHMM sensitive (retained)", "OH-S", "#56B4E9", "D"),
    ("orthohmm_phylogeny_satellite_v2", "OrthoHMM phylogeny (retained)", "OH-P", "#009E73", "s"),
    ("orthofinder_3_1_5_full", "OrthoFinder 3.1.5 full", "OF-F", "#D55E00", "^"),
    ("orthofinder_3_1_5_sequence_only", "OrthoFinder sequence checkpoint", "OF-S", "#E69F00", "v"),
    ("sonicparanoid_2_0_9", "SonicParanoid 2.0.9", "SP", "#CC79A7", "P"),
    ("proteinortho_6_3_6", "ProteinOrtho 6.3.6", "PO", "#666666", "X"),
    ("fastoma_0_3_5", "FastOMA 0.3.5", "FA", "#332288", "d"),
    ("orthomcl_1_4", "OrthoMCL 1.4", "OM", "#111111", "*"),
)


def record(path):
    path = Path(path)
    data = path.read_bytes()
    return dict(path=str(path), bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def write_json(path, value):
    with Path(path).open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def number(value):
    if value is None or value == "NA":
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError("Nonfinite plotted evidence")
    return result


def dominates(left, right, directions):
    if len(left) != len(right) or len(left) != len(directions) or not directions:
        raise ValueError("Dominance requires equal nonempty dimensions")
    if any(d not in (1, -1) for d in directions):
        raise ValueError("Directions must be +1 (maximize) or -1 (minimize)")
    if any(number(v) is None for v in (*left, *right)):
        raise ValueError("Missing objectives cannot participate in dominance")
    differences = [(float(a) - float(b)) * d for a, b, d in zip(left, right, directions)]
    return all(v >= 0 for v in differences) and any(v > 0 for v in differences)


def frontier(methods, objectives):
    keys = [m["key"] for m in methods]
    if len(keys) != len(set(keys)) or not objectives or len(set(objectives)) != len(objectives):
        raise ValueError("Unique method identities and objectives required")
    vectors = []
    for method in methods:
        try:
            vector = [number(method["scores"][k]) for k in objectives]
        except KeyError as error:
            raise ValueError("Incomplete objective vector") from error
        if any(v is None for v in vector):
            raise ValueError("Incomplete objective vector")
        vectors.append(vector)
    result = []
    for i, method in enumerate(methods):
        dominators = [other["key"] for j, other in enumerate(methods)
                      if i != j and dominates(vectors[j], vectors[i], [1] * len(objectives))]
        result.append(dict(key=method["key"], objectives=list(objectives),
                           nondominated=not dominators, dominators=dominators))
    return result


def uncertainty_row(label, values):
    intervals = {k: values[k] for k in ("paired_percentile_ci", "bonferroni_percentile_ci")}
    if any(len(v) != 2 or not all(math.isfinite(x) for x in v) or v[0] > v[1]
           for v in intervals.values()):
        raise ValueError("Invalid original uncertainty interval")
    return dict(label=label, difference=values["difference_percentage_points"], **intervals)


def validate_timing(panel):
    rows = panel["runs"]
    if len(rows) != 27 or len({r["index"] for r in rows}) != 27:
        raise ValueError("Timing panel must retain all 27 identities")
    if (sum(r["comparative_timing_eligible"] is True for r in rows) != 24
            or sum(r["resources"] is not None for r in rows) != 25):
        raise ValueError("Retain 24 eligible and 25 measured attempts")
    for row in rows:
        if row["resources"] is None:
            if row["comparative_timing_eligible"] or row["status"] != "reviewed_pre_native_abort":
                raise ValueError("Missing resources cannot be an eligible observation")
        elif any(number(v) is None or number(v) <= 0 for v in row["resources"].values()):
            raise ValueError("Measured resources must be finite and positive")
    for cell in panel["cells"]:
        selected = [r for r in rows if (r["method"], r["proteomes"]) == (cell["method"], cell["proteomes"])]
        eligible = [r for r in selected if r["comparative_timing_eligible"]]
        if len(selected) != 3 or len(eligible) != cell["eligible_repeats"]:
            raise ValueError("Changed timing-cell coverage")
        for metric, summary in cell["resources"].items():
            expected = None
            if len(eligible) == 3:
                values = [r["resources"][metric] for r in eligible]
                expected = dict(median=median(values), minimum=min(values), maximum=max(values))
            if summary != expected:
                raise ValueError("Cell summaries require exactly three eligible repeats")


def collect(repo):
    repo = repo.resolve()
    refs = []
    for path, digest in PINS.items():
        ref = record(repo / path)
        if ref["sha256"] != digest:
            raise ValueError("Changed admitted input: " + path)
        refs.append(ref)
    source = json.loads((repo / ARCHIVE / "summary.json").read_text())
    if (source["configurations"] != 29 or source["endpoint_slots"] != 232
            or source["scored_endpoint_observations"] != 125 or source["defaults_promoted"]):
        raise ValueError("Require the complete unchanged campaign")
    rows = list(csv.DictReader((repo / ARCHIVE / "endpoint_diagnostics.tsv").open(), delimiter="\t"))
    index = {(r["method"], r["dataset"] + ": " + r["endpoint"]): r for r in rows}
    if len(rows) != 232 or len(index) != 232:
        raise ValueError("Incomplete or duplicate diagnostics")
    originals = {m["key"]: m for m in [*source["context"], *source["methods"]]}
    baseline = originals["profile_weighting:none"]["scores"]
    aliases = ("canonical_collection:canonical", "projection_recovery:legacy_off")
    if any(originals[k]["scores"] != baseline for k in aliases):
        raise ValueError("Repeated canonical controls no longer match")
    methods = []
    for key, label, short, color, marker in METHODS:
        original = originals[key]
        if set(original["scores"]) != set(ENDPOINTS.values()):
            raise ValueError("A main configuration lacks a complete vector")
        method = dict(key=key, label=label, short=short, color=color, marker=marker,
                      scores={}, diagnostics={}, evidence=original["status"])
        for endpoint, name in ENDPOINTS.items():
            value = original["scores"][name]
            diag = index[key, name]
            if number(diag["score"]) != value or not 0 <= value <= 1:
                raise ValueError("Plot/diagnostic score ownership mismatch")
            method["scores"][endpoint] = value
            method["diagnostics"][endpoint] = {k: diag[k] for k in (
                "precision", "recall", "assessed_relations", "numeric_scored_pairs",
                "semantics", "native_genes", "native_species", "reference_gene_coverage")}
            if endpoint in ALL_F1:
                for k in ("precision", "recall"):
                    n = number(diag[k])
                    if n is None or not 0 <= n <= 1:
                        raise ValueError("Reference endpoint requires recorded precision/recall")
        methods.append(method)
    # This existing reader checks the complete frozen timing inventory and summaries.
    sys.path.insert(0, str(repo))
    from benchmark_tools.results.plot_shared_prepared_resources_20261004 import validate
    from benchmark_tools.prepare_ob_candidate_neighborhood import check
    timing = json.loads((repo / TIMING).read_text())
    validate(timing)
    validate_timing(timing)
    for ref in [timing["source"], timing["plan"], *timing["evidence"]]:
        check(ref)
    if timing["reviewed_attempts"] != 27 or timing["eligible_attempts"] != 24:
        raise ValueError("Require the complete final descriptive resource panel")
    standalone = []
    for key, label in (("henikoff", "Henikoff"), ("henikoff_cluster", "Clustered Henikoff")):
        standalone.append(uncertainty_row(label, source["uncertainty"]["weighting"]["comparisons"][key]["metrics"]["f_score"]))
    standalone.append(uncertainty_row("Projection repair", source["uncertainty"]["projection_recovery"]["versus_legacy"]["comparisons"]["projection_off"]["metrics"]["f_score"]))
    interaction = []
    for contrast in source["interaction_uncertainty"]["contrasts"]:
        candidate = contrast["candidate"]
        values = contrast["result"]["comparisons"][candidate]["metrics"]["f_score"]
        interaction.append(uncertainty_row(candidate + " vs " + contrast["control"], values))
    scopes = {"primary_three_F1": list(PRIMARY), "all_five_F1": list(ALL_F1)}
    scopes.update({"pair_" + a + "_" + b: [a, b] for a, b in (
        ("OB", "VGNC"), ("OB", "SwissTrees"), ("VGNC", "SwissTrees"))})
    return dict(schema_version=1, status="descriptive_posthoc_manuscript_analysis",
        methods=methods, sources=refs, collapsed_aliases={"profile_weighting:none": list(aliases)},
        all_scores_tsv=(repo / ARCHIVE / "scores.tsv").read_text(),
        endpoint_names=ENDPOINTS, scopes=scopes,
        pareto={key: frontier(methods, objectives) for key, objectives in scopes.items()},
        timing={k: timing[k] for k in ("runs", "cells", "planned_attempts", "reviewed_attempts",
            "resource_reviewed_attempts", "pre_native_aborted_indices",
            "eligible_attempts", "excluded_attempts", "primary_scopes", "execution_scope",
            "uncontended_timing", "scientific_timings_admitted", "limitations")},
        uncertainty=dict(standalone=standalone, interaction=interaction, replicates=100000,
            scope="Original phase-specific paired/Bonferroni intervals; development-exposed, not campaign-wide or independent confirmation."),
        defaults_promoted=False, independent_confirmation_established=False,
        limitations=["Point-estimate dominance is not statistical superiority.",
            "Functional means use different conditional populations and are excluded from accuracy frontiers.",
            "Timing revisions differ from current experimental configurations; no joint accuracy-runtime frontier.",
            "TreeFam supporting; Three Kingdoms limited marker co-membership, not genome-wide orthology truth."])


def setup():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 11,
        "axes.labelsize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "axes.spines.top": False, "axes.spines.right": False,
        "savefig.facecolor": "white", "figure.facecolor": "white"})


def style(ax):
    ax.grid(color="#e5e5e5", linewidth=.5)
    ax.set_axisbelow(True)


def legend(fig, methods, bottom=.02):
    handles = [Line2D([], [], marker=m["marker"], color=m["color"], linestyle="none",
                      markersize=6, label=m["short"] + "  " + m["label"]) for m in methods]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, bottom),
               frameon=False, ncol=3, fontsize=7.5, columnspacing=1.5, handletextpad=.4)


def point(ax, method, x, y, **kwargs):
    ax.scatter(x, y, marker=method["marker"], color=method["color"],
               s=55, edgecolors="white", linewidths=.5, zorder=5, **kwargs)


def front_line(ax, methods, coordinates):
    selected = []
    for method in methods:
        x, y = coordinates(method)
        if not any(dominates(coordinates(other), (x, y), (1, 1)) for other in methods):
            selected.append((x, y))
    selected = sorted(set(selected))
    if len(selected) > 1:
        ax.step(*zip(*selected), where="pre", color="#444444", linestyle="--", linewidth=.9, zorder=2)


def workflow(data):
    fig, ax = plt.subplots(figsize=(11.3, 5.3))
    ax.set(xlim=(0, 12), ylim=(0, 6)); ax.axis("off")
    labels = ["Protein\nsequences", "Built-in HMM\nhomology search", "Seed families +\nprofile recruitment",
              "Candidate-family\nexpansion", "Gene/species trees\n+ reconciliation", "Groups +\northolog pairs"]
    for i, label in enumerate(labels):
        x = .1 + i * 2
        ax.add_patch(Rectangle((x, 4.25), 1.8, 1.05, facecolor="#edf3f4", edgecolor="#50757b", linewidth=.9))
        ax.text(x + .9, 4.78, label, ha="center", va="center", fontsize=9)
        if i < 5:
            ax.add_patch(FancyArrowPatch((x + 1.8, 4.78), (x + 2, 4.78), arrowstyle="-|>", mutation_scale=12, color="#555555"))
    ax.text(.1, 5.75, "A  HMM-centered inference with optional phylogenetic refinement", weight="bold", fontsize=12)
    ax.text(.1, 3.65, "B  Four separately controlled interventions", weight="bold", fontsize=12)
    interventions = [("Sequence weighting", "Relative residue weights;\nactual versus cluster effective N"),
        ("Alignment-aware recruitment", "Aligned residues, fragments\nand fusion/fission asymmetry"),
        ("Local joint-tree rescue", "Bounded neighboring families;\nreconciliation, not automatic union"),
        ("Selective HMMER validation", "Ambiguous candidates; separate\nmeasurement from reselection")]
    for i, (title, body) in enumerate(interventions):
        x = .1 + i * 3
        ax.add_patch(Rectangle((x, 1.45), 2.8, 1.55, facecolor="#fafafa", edgecolor="#bbbbbb", linewidth=.8))
        ax.text(x + .14, 2.64, title, fontsize=9, weight="bold")
        ax.text(x + .14, 2.1, body, fontsize=8.5, va="center", linespacing=1.6)
    ax.text(.1, .65, "Individual ablations  >  fixed advancement gates  >  justified interaction  >  full-vector reporting", fontsize=10)
    ax.text(.1, .1, "HMM support nominates candidates; it does not alone establish orthology. Benchmark labels are evaluation-only.", fontsize=9, color="#555555")
    fig.subplots_adjust(left=.02, right=.99, top=.98, bottom=.04)
    return fig


def orthobench(data):
    methods = data["methods"]
    fig, axes = plt.subplots(1, 2, figsize=(11.3, 6.8), gridspec_kw={"width_ratios": [1.2, 1]})
    for i, m in enumerate(methods):
        score = m["scores"]["OB"]
        axes[0].barh(i, score, color=m["color"], height=.62)
        axes[0].text(score + .006, i, f"{score:.4f}", va="center", fontsize=8)
        point(axes[1], m, number(m["diagnostics"]["OB"]["recall"]), number(m["diagnostics"]["OB"]["precision"]))
    axes[0].set(yticks=range(len(methods)), yticklabels=[m["short"] for m in methods],
        xlim=(0, .85), xlabel="Weighted group F1", title="A  OrthoBench group recovery")
    axes[0].invert_yaxis()
    front_line(axes[1], methods, lambda m: (number(m["diagnostics"]["OB"]["recall"]), number(m["diagnostics"]["OB"]["precision"])))
    axes[1].set(xlim=(0, 1.02), ylim=(0, 1.02), xlabel="Weighted recall", ylabel="Weighted precision", title="B  Precision-recall trade-off")
    for m in methods:
        if m["short"] in ("OH-C", "OF-F"):
            label = "OH-C/P" if m["short"] == "OH-C" else m["short"]
            axes[1].annotate(label, (number(m["diagnostics"]["OB"]["recall"]), number(m["diagnostics"]["OB"]["precision"])), xytext=(5, 8), textcoords="offset points", fontsize=8)
    for ax in axes: style(ax)
    legend(fig, methods)
    fig.subplots_adjust(left=.07, right=.98, top=.94, bottom=.26, wspace=.30)
    return fig


def pareto_plot(data):
    methods = data["methods"]
    fig, axes = plt.subplots(2, 2, figsize=(11.3, 9))
    for ax, (a, b), letter in zip(axes.flat, (("OB", "VGNC"), ("OB", "SwissTrees"), ("VGNC", "SwissTrees")), "ABC"):
        for m in methods: point(ax, m, m["scores"][a], m["scores"][b])
        front_line(ax, methods, lambda m: (m["scores"][a], m["scores"][b]))
        ax.set(xlim=(0, 1.02), ylim=(0, 1.02), xlabel=a + " F1", ylabel=b + " F1", title=letter + "  Pairwise point-estimate frontier")
        for m in methods:
            if m["short"] in ("OF-F", "OH-C"):
                offset = (-45, 10) if m["short"] == "OF-F" else (-45, -23)
                if letter == "A" and m["short"] == "OF-F":
                    offset = (-45, -13)
                ax.annotate(m["short"], (m["scores"][a], m["scores"][b]), xytext=offset, textcoords="offset points", fontsize=8,
                    arrowprops=dict(arrowstyle="-", linewidth=.5, color=m["color"]))
        style(ax)
    ax = axes[1, 1]
    selected = [m for m in methods if m["short"] in ("OH-C", "OH-W", "OH-N", "OH-R", "OH-P")]
    for m in selected:
        point(ax, m, m["scores"]["OB"], m["scores"]["VGNC"])
        dy = {"OH-W": -13, "OH-N": 9}.get(m["short"], 5)
        ax.annotate(m["short"], (m["scores"]["OB"], m["scores"]["VGNC"]), xytext=(7, dy), textcoords="offset points", fontsize=8)
    front_line(ax, methods, lambda m: (m["scores"]["OB"], m["scores"]["VGNC"]))
    ax.set(xlim=(.74095, .7422), ylim=(.90065, .9023), xlabel="OB F1 (magnified)", ylabel="VGNC F1 (magnified)", title="D  Small development differences")
    ax.ticklabel_format(useOffset=False, axis="both"); style(ax)
    members = {r["key"] for r in data["pareto"]["primary_three_F1"] if r["nondominated"]}
    fig.suptitle("Primary three-objective frontier: " + ", ".join(m["short"] for m in methods if m["key"] in members), fontsize=10, y=.992)
    legend(fig, methods)
    fig.subplots_adjust(left=.075, right=.98, top=.95, bottom=.20, hspace=.40, wspace=.32)
    return fig


def resources(data):
    fig, axes = plt.subplots(1, 3, figsize=(11.3, 5.8))
    panel = data["timing"]
    mapping = (("orthohmm_high_sensitivity", "OH-S"), ("orthohmm_satellite_v2", "OH-P"), ("orthofinder_3_1_5_full", "OF-F"))
    styles = {m["short"]: m for m in data["methods"]}
    metrics = (("wall_seconds", 60, "A  Native elapsed time", "Minutes"),
        ("cpu_seconds", 3600, "B  Native CPU use", "CPU-hours"),
        ("peak_memory_bytes", 1024 ** 3, "C  Native-step memory peak", "GiB (includes launcher)"))
    for ax, (metric, divisor, title, ylabel) in zip(axes, metrics):
        for index, (key, short) in enumerate(mapping):
            m = styles[short]
            for sizeindex, size in enumerate((4, 8, 12)):
                x = sizeindex + (index - 1) * .24
                for row in panel["runs"]:
                    if (row["method"] != key or row["proteomes"] != size
                            or row["resources"] is None or row["resources"][metric] is None):
                        continue
                    eligible = row["comparative_timing_eligible"]
                    ax.scatter(x + (row["repeat"] - 1) * .045, row["resources"][metric] / divisor,
                        marker=m["marker"] if eligible else "x", color=m["color"] if eligible else "#777777", s=30, zorder=4)
                cell = next(c for c in panel["cells"] if c["method"] == key and c["proteomes"] == size)
                summary = cell["resources"][metric]
                if summary is not None:
                    ax.plot([x, x], [summary["minimum"] / divisor, summary["maximum"] / divisor], color=m["color"], linewidth=1.4)
                    ax.plot([x - .06, x + .06], [summary["median"] / divisor] * 2, color=m["color"], linewidth=2)
                ax.text(x, 1.035, f"{cell['eligible_repeats']}/3", transform=ax.get_xaxis_transform(), ha="center", fontsize=7, color=m["color"])
        ax.set(xticks=(0, 1, 2), xticklabels=(4, 8, 12), xlabel="Complete proteomes", ylabel=ylabel, xlim=(-.5, 2.5), ylim=(0, None))
        ax.set_title(title, pad=30); style(ax)
    handles = [Line2D([], [], color=styles[short]["color"], marker=styles[short]["marker"], linestyle="none", label=styles[short]["label"]) for _, short in mapping]
    handles.append(Line2D([], [], color="#777777", marker="x", linestyle="none", label="Excluded measured attempt"))
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, .07), ncol=2, frameon=False, fontsize=8)
    fig.text(.5, .01, "Shared host; 32 native CPUs, 128 GiB/run. Median/range only for 3/3 eligible repeats. No uncontended speed ranking.", ha="center", fontsize=8)
    fig.subplots_adjust(left=.07, right=.98, top=.80, bottom=.28, wspace=.40)
    return fig


def qfo(data):
    methods = data["methods"]
    fig, axes = plt.subplots(2, 3, figsize=(11.3, 8.6))
    for ax, endpoint, letter in zip(axes[0], ("VGNC", "SwissTrees", "TreeFam"), "ABC"):
        for m in methods:
            d = m["diagnostics"][endpoint]
            point(ax, m, number(d["recall"]), number(d["precision"]))
        front_line(ax, methods, lambda m: (number(m["diagnostics"][endpoint]["recall"]), number(m["diagnostics"][endpoint]["precision"])))
        ax.set(xlim=(0, 1.03), ylim=(0, 1.03), xlabel="Recall", ylabel="Precision", title=letter + "  " + endpoint + (" (supporting)" if endpoint == "TreeFam" else "")); style(ax)
    for ax, endpoint, letter in zip(axes[1], ("EC", "GO", "FAS"), "DEF"):
        for m in methods:
            count = number(m["diagnostics"][endpoint]["assessed_relations"])
            if count is None or count <= 0:
                raise ValueError("Functional plot requires an explicit assessed/eligible population")
            point(ax, m, count, m["scores"][endpoint])
        ax.set_xscale("log")
        ax.set(xlabel="Eligible pair pool (not sample n)" if endpoint == "FAS" else "Assessed relations", ylabel="Conditional mean similarity",
               ylim=(0, 1.03), title=letter + "  " + endpoint + " conservation"); style(ax)
    legend(fig, methods)
    fig.subplots_adjust(left=.07, right=.98, top=.95, bottom=.21, hspace=.42, wspace=.36)
    return fig


def heatmap(data):
    methods = data["methods"]
    matrix = np.array([[m["scores"][k] for k in ENDPOINTS] for m in methods])
    fig, ax = plt.subplots(figsize=(11.3, 6.5))
    image = ax.imshow(matrix, cmap="cividis", vmin=0, vmax=1, aspect="auto")
    for i in range(len(methods)):
        for j in range(8):
            ax.text(j, i, f"{matrix[i, j]:.4f}", ha="center", va="center", fontsize=8, color="white" if matrix[i, j] < .55 else "#111111")
    ax.set(xticks=range(8), xticklabels=("OB F1", "VGNC F1", "Swiss F1", "TreeFam F1", "EC mean", "GO mean", "FAS mean", "TK F1"),
        yticks=range(len(methods)), yticklabels=[m["short"] + "  " + m["label"] for m in methods])
    ax.tick_params(length=0); fig.colorbar(image, ax=ax, shrink=.8, label="Score (not a shared accuracy scale)")
    fig.subplots_adjust(left=.30, right=.98, bottom=.10, top=.97)
    return fig


def uncertainty(data):
    fig, axes = plt.subplots(1, 2, figsize=(11.3, 5.5), gridspec_kw={"width_ratios": [1, 1.35]})
    for ax, name, letter in zip(axes, ("standalone", "interaction"), "AB"):
        rows = data["uncertainty"][name]
        for i, row in enumerate(rows):
            ax.plot(row["bonferroni_percentile_ci"], [i, i], color="#b6b6b6", linewidth=4, solid_capstyle="butt")
            ax.plot(row["paired_percentile_ci"], [i, i], color="#0072B2", linewidth=2)
            ax.scatter(row["difference"], i, color="#111111", s=22, zorder=4)
        ax.axvline(0, color="#555555", linestyle="--", linewidth=.8)
        ax.set(yticks=range(len(rows)), yticklabels=[r["label"].replace("repair_henikoff_cluster", "N + repair").replace("repair_henikoff", "W + repair").replace("legacy_henikoff_cluster", "N alone").replace("legacy_henikoff", "W alone").replace("legacy_none", "canonical").replace("repair_none", "repair alone") for r in rows],
            xlabel="OrthoBench F1 difference (percentage points)", title=letter + ("  Standalone versus canonical" if name == "standalone" else "  Registered interaction contrasts"))
        ax.invert_yaxis(); style(ax)
    fig.legend(handles=[Line2D([], [], color="#0072B2", linewidth=2, label="Paired percentile interval"),
        Line2D([], [], color="#b6b6b6", linewidth=4, label="Original phase-adjusted interval")], loc="lower center", bbox_to_anchor=(.5, .01), ncol=2, frameon=False, fontsize=8)
    fig.subplots_adjust(left=.16, right=.98, top=.90, bottom=.21, wspace=.85)
    return fig


def export(data, output):
    if output.exists():
        raise FileExistsError("Use a fresh output directory; preserve previous figure attempts")
    if len(data["methods"]) != 12 or len({m["key"] for m in data["methods"]}) != 12:
        raise ValueError("Require the complete fixed 12-configuration panel")
    for key, objectives in data["scopes"].items():
        if frontier(data["methods"], objectives) != data["pareto"][key]:
            raise ValueError("Snapshot frontier differs from exact recomputation")
    validate_timing(data["timing"])
    setup(); output.mkdir(parents=True)
    write_json(output / "analysis.json", data)
    with (output / "scores_all_29.tsv").open("x") as stream:
        stream.write(data["all_scores_tsv"])
    with (output / "pareto_membership.tsv").open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=("scope", "key", "objectives", "nondominated", "dominators"), delimiter="\t")
        writer.writeheader()
        for scope, values in data["pareto"].items():
            for row in values:
                writer.writerow({"scope": scope, **row, "objectives": json.dumps(row["objectives"]), "dominators": json.dumps(row["dominators"])})
    with (output / "plot_coordinates.tsv").open("x", newline="") as stream:
        columns = ("key", "short", "endpoint", "score", "precision", "recall", "assessed_relations", "numeric_scored_pairs", "semantics")
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", extrasaction="ignore"); writer.writeheader()
        for m in data["methods"]:
            for endpoint in ENDPOINTS:
                writer.writerow(dict(key=m["key"], short=m["short"], endpoint=endpoint, score=m["scores"][endpoint], **m["diagnostics"][endpoint]))
    with (output / "timing_attempts.tsv").open("x", newline="") as stream:
        columns = ("index", "job_id", "method", "proteomes", "repeat", "status", "comparative_timing_eligible", "wall_seconds", "cpu_seconds", "peak_memory_bytes")
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", extrasaction="ignore"); writer.writeheader()
        for row in data["timing"]["runs"]:
            writer.writerow({**row, **(row["resources"] or {})})
    figures = (("fig1_workflow", workflow), ("fig2_orthobench", orthobench),
        ("fig3_pareto", pareto_plot), ("fig4_resources", resources), ("fig5_qfo", qfo),
        ("figS1_all_scores", heatmap), ("figS2_uncertainty", uncertainty))
    for name, draw in figures:
        fig = draw(data)
        try:
            for extension in ("png", "pdf", "svg"):
                fig.savefig(output / (name + "." + extension), dpi=220)
        finally:
            plt.close(fig)
    refs = [record(p) for p in sorted(output.iterdir()) if p.is_file()]
    write_json(output / "build_manifest.json", dict(status="original_figures_and_exact_frontiers_exported",
        builder=record(__file__), sources=data["sources"], outputs=refs,
        numpy_version=np.__version__, matplotlib_version=matplotlib.__version__,
        scope="No new inference, scoring, resampling or default promotion; rendered layout review is separate."))
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--repo", type=Path)
    source.add_argument("--snapshot", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = collect(args.repo) if args.repo else json.loads(args.snapshot.read_text())
    export(data, args.output)
    print(json.dumps({"output": str(args.output), "primary_frontier": [r["key"] for r in data["pareto"]["primary_three_F1"] if r["nondominated"]]}))
