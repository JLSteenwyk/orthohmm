"""Plot verified native SwissTrees stage separation and frozen-bin trade-offs."""

import argparse
from collections import Counter
import csv
import json
import math
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.export_native_factorial_progress import load, require
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

CELLS = ("p0_c0_r0", "p0_c0_r1")
BINS = ("all", "lower_entropy", "higher_entropy", "short_relative", "not_short_relative")
LABELS = ("All families (18)", "Lower entropy (9)", "Higher entropy (9)", "Short-relative (7)", "Other families (11)")
COLORS = ("#12786f", "#b64b37")
SUPPORT = ("no_direct_hit", "one_direction", "both_directions")
KINDS = {
    "search": ("native_qfo_swiss_direct_search_support_v1", "trace_native_qfo_swiss_search_support.py",
               "native_qfo_swiss_search_support_code_readback_v1", "readback_native_qfo_swiss_search_support.py"),
    "reconciliation": ("native_qfo_swiss_reconciliation_localization_v1", "trace_native_qfo_swiss_reconciliation.py",
                       "native_qfo_swiss_reconciliation_newick_readback_v1", "readback_native_qfo_swiss_reconciliation.py"),
    "strata": ("native_qfo_swiss_sequence_strata_v1", "export_native_qfo_swiss_sequence_strata.py",
               "native_qfo_swiss_sequence_strata_readback_v1", "readback_native_qfo_swiss_sequence_strata.py"),
}


def paired_reports(kind, report_ref, readback_ref, evidence):
    report, observed = load(report_ref["path"], report_ref["sha256"], evidence)
    readback, verified = load(readback_ref["path"], readback_ref["sha256"], evidence)
    schema, source, readback_schema, readback_source = KINDS[kind]
    require(observed == report_ref and verified == readback_ref and report["schema"] == schema
            and readback["schema"] == readback_schema and readback["report"] == observed
            and report["source"] == record(Path(__file__).with_name(source))
            and readback["source"] == record(Path(__file__).with_name(readback_source))
            and all(item["publication_ready"] is False and item["independent_confirmation"] is False
                    for item in (report, readback)), "Wrong diagnostic/readback source or scope: " + kind)
    for index, item in enumerate((report, readback)):
        flags = (("new_uncertainty", "new_accuracy_or_resource_admission") if kind == "strata" else
                 ("new_scoring_or_admission", "uncertainty_admitted", "scientific_timings_admitted"))
        if kind == "reconciliation" and index == 1:
            flags = ("new_scoring_or_admission", "uncertainty_admitted")
        require(all(item[k] is False for k in flags), "Diagnostic claims new inference: " + kind)
        evidence.append(item["source"])
    return report


def figure_data(search, reconciliation, strata):
    cases = search["cases"]
    require(len(cases) == search["changed_pairs"] == reconciliation["changed_pairs_traced"] == 2023,
            "Changed-pair universe differs")
    require(len({(r["family"], r["protein_a"], r["protein_b"]) for r in cases}) == len(cases)
            and all(r["before"] in ("TP", "FP") and r["after"] == ("FN" if r["before"] == "TP" else "TN")
                    for r in cases), "Invalid pair identities/truth transitions")
    counts, roots = Counter(), Counter()
    for row in cases:
        require(type(row["same_root_hog"]) is bool and row["selected_directed_score_multisets_identical"] is True
                and set(row["direct_search"]) == set(CELLS), "Wrong native search case scope")
        roots[row["before"], row["same_root_hog"]] += 1
        for cell in CELLS:
            view = row["direct_search"][cell]
            forward, reverse = view["gene_a_to_b"], view["gene_b_to_a"]
            category = "both_directions" if forward and reverse else "one_direction" if forward or reverse else "no_direct_hit"
            require(view["support"] == category, "Incorrect direct-search classification")
            require(all(type(h["row"]) is int and h["row"] >= 0 and type(h["score"]) in (int, float)
                        and math.isfinite(h["score"]) for hits in (forward, reverse) for h in hits), "Invalid hit record")
            counts[row["before"], cell, category] += 1
        require(all(sorted(h["score"] for h in row["direct_search"][CELLS[0]][direction]) ==
                    sorted(h["score"] for h in row["direct_search"][CELLS[1]][direction])
                    for direction in ("gene_a_to_b", "gene_b_to_a")), "Selected search scores differ")
    expected_summary = [dict(before=label, cell=cell, support=category, pairs=counts[label, cell, category])
        for label in ("TP", "FP") for cell in CELLS for category in SUPPORT]
    require(search["summary"] == expected_summary, "Search summary differs from cases")
    require(all(roots[label, same] == reconciliation["summary"].get(label + ("_same_root" if same else "_different_root"), 0)
                for label in ("TP", "FP") for same in (True, False)), "Reconciliation RootHOG counts differ")
    stage_rows = []
    for label in ("TP", "FP"):
        total = sum(counts[label, CELLS[0], category] for category in SUPPORT)
        require(total == (334 if label == "TP" else 1689), "Changed excluded-positive cohort size")
        for category in SUPPORT:
            require(counts[label, CELLS[0], category] == counts[label, CELLS[1], category], "Different search support views")
            stage_rows.append(dict(panel="A", label=label, category=category, count=counts[label, CELLS[0], category],
                                   total=total, percentage=100 * counts[label, CELLS[0], category] / total))
        for same in (True, False):
            stage_rows.append(dict(panel="B", label=label, category="same_root" if same else "different_root",
                                   count=roots[label, same], total=total, percentage=100 * roots[label, same] / total))
    all_rows = {(r["cell"], r["stratum"]): r for r in strata["rows"]}
    require(len(all_rows) == len(strata["rows"]) == 22 and len(strata["differences"]) == 11, "Wrong native strata inventory")
    deltas = {r["stratum"]: r for r in strata["differences"]}
    require(len(deltas) == 11 and set(BINS) <= set(deltas)
            and all((cell, name) in all_rows for cell in CELLS for name in BINS), "Missing or duplicate plotted strata")
    plotted, changes = [], []
    for name, size in zip(BINS, (18, 9, 9, 7, 11)):
        for cell in CELLS:
            row = all_rows[cell, name]
            require(row["status"] == "descriptive" and row["families"] == len(row["family_members"]) == size
                    and all(type(row[k]) in (int, float) and math.isfinite(row[k]) and 0 < row[k] <= 1
                            for k in ("F1", "PPV", "TPR")), "Invalid plotted native stratum")
            require(math.isclose(row["F1"], 2 * row["PPV"] * row["TPR"] / (row["PPV"] + row["TPR"]),
                                 rel_tol=0, abs_tol=1e-12), "Plotted F1 arithmetic differs")
            plotted.append(dict(cell=cell, stratum=name, families=size, F1=row["F1"], PPV=row["PPV"], TPR=row["TPR"]))
        left, right = (all_rows[cell, name] for cell in CELLS)
        require(left["family_members"] == right["family_members"], "Different native stratum membership")
        for metric in ("PPV", "TPR"):
            delta = right[metric] - left[metric]
            require(deltas[name][metric] == delta, "Plotted stratum difference differs")
            changes.append(dict(stratum=name, metric=metric, difference_pp=100 * delta))
    return stage_rows, plotted, changes


def render(stage, scores, changes, output):
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 10, "svg.hashsalt": "native-swiss-mechanism",
                         "svg.fonttype": "none", "pdf.fonttype": 42}):
        fig, axes = plt.subplots(2, 2, figsize=(13.2, 8.8))
        for axis, panel, title, categories, colors, names in (
            (axes[0, 0], "A", "A  Direct significant search support", SUPPORT,
             ("#b6bfc4", "#e4ad48", "#26859a"), ("No direct hit", "One direction", "Both directions")),
            (axes[0, 1], "B", "B  Final membership of excluded pairs", ("same_root", "different_root"),
             ("#576774", "#c7cdd0"), ("Same RootHOG", "Different RootHOG"))):
            for y, label in enumerate(("TP", "FP")):
                left = 0
                for category, color in zip(categories, colors):
                    row = next(r for r in stage if (r["panel"], r["label"], r["category"]) == (panel, label, category))
                    axis.barh(y, row["percentage"], left=left, height=.45, color=color, edgecolor="white", linewidth=1)
                    if row["percentage"] >= 8:
                        axis.text(left + row["percentage"] / 2, y, str(row["count"]), ha="center", va="center",
                                  color="white" if color in ("#26859a", "#576774") else "#222222", fontsize=10)
                    left += row["percentage"]
            axis.set(xlim=(0, 100), ylim=(1.6, -.6), xlabel="Percentage of removed reference pairs")
            axis.set_yticks((0, 1), ("Removed TP (334)", "Removed FP (1,689)"))
            axis.legend([plt.Rectangle((0, 0), 1, 1, color=c) for c in colors], names,
                        frameon=False, fontsize=9, loc="lower left", bbox_to_anchor=(-.02, 1.04),
                        borderaxespad=0, ncol=len(names))
            axis.set_title(title, loc="left", fontsize=11, pad=40)
            axis.set_xticks((0, 25, 50, 75, 100))
        values = {(r["cell"], r["stratum"]): r for r in scores}
        axis = axes[1, 0]
        for y, name in enumerate(BINS):
            points = [100 * values[cell, name]["F1"] for cell in CELLS]
            axis.plot(points, (y, y), color="#a8a8a8", linewidth=1.4)
            for point, color in zip(points, COLORS):
                axis.scatter(point, y, color=color, s=55, edgecolor="white", linewidth=.5, zorder=3)
        axis.set(xlim=(50, 100), ylim=(4.55, -.55), xlabel="F1 (%)")
        axis.set_yticks(range(5), LABELS)
        axis.set_title("C  Frozen sequence-strata F1", loc="left", fontsize=11, pad=40)
        axis.grid(axis="x", color="#e8e8e8", linewidth=.7)
        axis.legend([plt.Line2D([], [], marker="o", linestyle="none", color=c) for c in COLORS],
                    ("R-off: group-clique", "R-on: resolved pairs"), frameon=False, fontsize=9,
                    loc="lower left", bbox_to_anchor=(-.02, 1.04), borderaxespad=0, ncol=2)
        axis = axes[1, 1]
        axis.axvline(0, color="#777777", linewidth=.8)
        for metric, offset, color in (("PPV", -.13, COLORS[0]), ("TPR", .13, COLORS[1])):
            points = [next(r["difference_pp"] for r in changes if r["stratum"] == name and r["metric"] == metric)
                      for name in BINS]
            axis.barh([y + offset for y in range(5)], points, height=.22, color=color, label="Precision" if metric == "PPV" else "Recall")
        axis.set(xlim=(-16, 42), ylim=(4.55, -.55), xlabel="R-on minus R-off (percentage points)")
        axis.set_yticks(range(5), LABELS)
        axis.set_title("D  Precision-recall trade-offs", loc="left", fontsize=11, pad=40)
        axis.grid(axis="x", color="#e8e8e8", linewidth=.7)
        axis.legend(frameon=False, fontsize=9, loc="lower left", bbox_to_anchor=(-.02, 1.04), borderaxespad=0, ncol=2)
        for axis in axes.flat:
            axis.spines[["top", "right"]].set_visible(False)
            axis.set_axisbelow(True)
            axis.tick_params(length=3)
        fig.suptitle("Native SwissTrees: search support and reconciliation trade-offs", x=.045, ha="left", fontsize=15, y=.98)
        fig.text(.045, .932, "Initial HMM search on; P = 0, C = 0. A-B: all 2,023 changed pairs. C-D: unchanged input-only family bins.", fontsize=10)
        fig.text(.045, .065, "A: selected direct-search scores identical in both views; homology support is not orthology confidence. Small one-direction counts: TP 6, FP 20.", fontsize=8.8)
        fig.text(.045, .041, "18 development-exposed families; descriptive overlapping bins, no subgroup intervals or significance claims. Entropy is not divergence; relative shortness is not fragmentation.", fontsize=8.8)
        fig.text(.045, .017, "Not selected-default superiority or a total-HMM ablation. Failed R-on native timing remains ineligible; no inference timing is displayed.", fontsize=8.8)
        fig.subplots_adjust(left=.155, right=.985, top=.79, bottom=.17, hspace=.82, wspace=.70)
        for suffix in ("png", "pdf", "svg"):
            fig.savefig(output / ("native_swiss_mechanism." + suffix), dpi=200,
                        metadata={"Creator": "OrthoHMM benchmark workflow"})
        plt.close(fig)


def run(inputs, output):
    require(set(inputs) == set(KINDS) and not output.exists() and not output.is_symlink(), "Invalid inputs or occupied output")
    evidence = [record(Path(__file__).parent / "results/NATIVE_QFO_SWISS_MECHANISM_FIGURE_PROTOCOL_20261006.md")]
    reports = {kind: paired_reports(kind, *refs, evidence) for kind, refs in inputs.items()}
    search, reconciliation, strata = (reports[k] for k in ("search", "reconciliation", "strata"))
    require(search["localization"] == inputs["reconciliation"][0], "Search/localization binding differs")
    transition, _ = load(reconciliation["transition"]["path"], reconciliation["transition"]["sha256"], evidence)
    binding, _ = load(strata["binding"]["path"], strata["binding"]["sha256"], evidence)
    require([r["cell"] for r in transition["cells"]] == list(CELLS)
            and all(r["admission"] == binding["bound_cells"][r["cell"]]["admission"] for r in transition["cells"]),
            "Different native scientific admissions")
    check(reconciliation["pair_ledger"])
    evidence.append(reconciliation["pair_ledger"])
    with open(reconciliation["pair_ledger"]["path"], newline="") as stream:
        localized = list(csv.DictReader(stream, delimiter="\t"))
    require([(r["family"], r["protein_a"], r["protein_b"], r["before"], r["after"], r["same_root_hog"]) for r in localized] ==
            [(r["family"], r["protein_a"], r["protein_b"], r["before"], r["after"], str(r["same_root_hog"])) for r in search["cases"]],
            "Different localized case identities or RootHOG flags")
    stage, scores, changes = figure_data(search, reconciliation, strata)
    for ref in evidence:
        check(ref)
    output.mkdir(parents=True)
    for name, rows in (("stages.tsv", stage), ("strata.tsv", scores), ("differences.tsv", changes)):
        with (output / name).open("x", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    render(stage, scores, changes, output)
    for ref in evidence:
        check(ref)
    result = dict(schema="native_qfo_swiss_mechanism_figure_v1", inputs=inputs, evidence=evidence,
        source=record(__file__), outputs=[record(p) for p in sorted(output.iterdir())],
        panels=4, changed_pairs=2023, stage_rows=len(stage), strata_rows=len(scores), difference_rows=len(changes),
        scientific_evidence_replayed=False, raw_search_rescanned=False, new_uncertainty=False,
        new_scoring_or_admission=False, scientific_timings_admitted=False, publication_ready=False,
        visual_review_complete=False, python_version=sys.version, matplotlib_version=matplotlib.__version__,
        limitations=["Direct source/report/readback bindings reused; no new array/tree/raw scientific replay.",
                     "Complete changed-pair cohort and frozen development-exposed descriptive strata, not biological causality.",
                     "No subgroup intervals, new uncertainty, selected-default comparison or total-HMM ablation.",
                     "Rendering is shared-host postprocessing with unknown potentially tool-dependent contention."])
    with (output / "manifest.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for kind in KINDS:
        for suffix in ("report", "readback"):
            parser.add_argument("--" + kind + "-" + suffix, nargs=2, required=True, metavar=("PATH", "SHA256"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    inputs = {}
    for kind in KINDS:
        refs = []
        for suffix in ("report", "readback"):
            path, digest = getattr(args, kind + "_" + suffix)
            ref = record(path)
            require(ref["sha256"] == digest, "Changed supplied figure input")
            refs.append(ref)
        inputs[kind] = refs
    result = run(inputs, args.output)
    print(json.dumps({k: result[k] for k in ("panels", "changed_pairs", "stage_rows", "strata_rows", "difference_rows")}))
