"""Present three admitted native QfO cells without imputing missing contrasts."""

import argparse
import csv
import json
import math
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools import plot_native_qfo_scientific_scores as previous
from benchmark_tools.export_native_factorial_progress import finite, load, require
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

CELLS = ("p0_c0_r0", "p0_c0_r1", "p0_c1_r0")
CONTRASTS = ("C_at_P0_R0", "R_at_P0_C0")
COLORS = ("#12786f", "#b64b37", "#4653a6")
ENDPOINTS = ("VGNC", "SwissTrees", "TreeFam-A", "GO", "EC", "FAS")
STEM = "native_qfo_three_cell"


def figure_data(snapshot, binding, readback):
    require(snapshot["schema"] == "native_qfo_scientific_reporting_snapshot_v1"
            and snapshot["publication_ready"] is False and len(snapshot["rows"]) == 7,
            "Require partial native scientific snapshot")
    rows = {row["cell"]: row for row in snapshot["rows"]}
    require(len(rows) == 7 and set(CELLS) <= set(rows)
            and sum(row["accuracy_admitted"] is True for row in rows.values()) == 3,
            "Require exactly three distinct admitted cells")
    scores, coverage = [], []
    for cell in CELLS:
        row = rows[cell]
        require(row["accuracy_admitted"] is True and row["status"] in (
            "supplied_native_admission", "supplied_recovered_scientific_admission")
            and set(row["scores"]) == set(ENDPOINTS), "Unadmitted or incomplete plotted cell")
        if cell == CELLS[1]:
            require(row["status"] == "supplied_recovered_scientific_admission"
                    and row["resources"] is None and row["timing_admitted"] is False
                    and row["timing_eligible"] is False, "Recovered plot relabels failed timing")
        for endpoint in ENDPOINTS:
            scores.append(dict(cell=cell, endpoint=endpoint,
                statistic="F1" if endpoint in ENDPOINTS[:3] else "similarity",
                value=finite(row["scores"][endpoint], endpoint, high=1)))
        for endpoint in ENDPOINTS[:3]:
            detail = row["endpoint_details"][endpoint]
            p, r = (finite(detail[k], k, high=1) for k in ("precision", "recall"))
            require(math.isclose(row["scores"][endpoint], 2*p*r/(p+r) if p+r else 0,
                                 rel_tol=0, abs_tol=1e-12), "Native precision/recall do not reproduce F1")
        n, covered, count = (row[k] for k in ("input_accessions", "relation_accessions", "submitted_pairs"))
        require(type(n) is int and n > 0 and type(covered) is int and 0 <= covered <= n
                and type(count) is int and count >= 0 and covered <= 2*count
                and (count == 0 or covered >= 2)
                and finite(row["relation_coverage"], "coverage", high=1) == covered/n,
                "Changed all-input coverage denominator")
        coverage.append(dict(cell=cell, input_accessions=n, relation_accessions=covered,
            submitted_pairs=count, relation_coverage=covered/n,
            precision=row["endpoint_details"]["SwissTrees"]["precision"],
            recall=row["endpoint_details"]["SwissTrees"]["recall"],
            prediction_semantics=row["prediction_semantics"]))
    for cell, row in rows.items():
        if cell not in CELLS:
            require(row["accuracy_admitted"] is False and all(v is None for v in row["scores"].values()),
                    "Missing plotted-panel score must remain null")
    require(binding["schema"] == "native_qfo_retained_swiss_uncertainty_binding_v1"
            and binding["publication_ready"] is False and binding["independent_confirmation"] is False
            and binding["multiplicity_endpoints"] == 42 and binding["replicates_reused"] == 100000
            and type(binding["new_bootstrap_draws"]) is int and binding["new_bootstrap_draws"] == 0,
            "Changed uncertainty scope")
    require(readback["schema"] == "native_qfo_candidate_swiss_rational_readback_v1"
            and readback["families_checked"] == 18 and readback["candidate_pair_labels_matched"] == 10765
            and all(readback[k] is False for k in ("publication_ready", "independent_confirmation",
                "new_accuracy_or_resource_admission", "scientific_timings_admitted"))
            and type(readback["new_bootstrap_draws"]) is int and readback["new_bootstrap_draws"] == 0,
            "Changed candidate readback scope")
    intervals = []
    require(len(binding["contrasts"]) == 14, "Changed contrast inventory")
    for name, candidate in zip(CONTRASTS, (CELLS[2], CELLS[1])):
        matches = [r for r in binding["contrasts"] if r["name"] == name]
        require(len(matches) == 1, "Missing or duplicated plotted contrast")
        effect = matches[0]
        require(effect["status"] == "native_records_matched" and not effect["missing_cells"]
                and not effect["differing_cells"] and effect["candidate"] == candidate
                and effect["reference"] == CELLS[0] and set(effect["metrics"]) == {"F1", "PPV", "TPR"},
                "Unmatched or changed contrast")
        for metric in ("F1", "PPV", "TPR"):
            value = effect["metrics"][metric]
            delta = finite(value["difference"], "difference", low=-1, high=1)
            nominal, adjusted = value["paired_percentile_ci"], value["bonferroni_percentile_ci"]
            require(len(nominal) == len(adjusted) == 2, "Changed interval shape")
            for bound in [*nominal, *adjusted]:
                finite(bound, "interval", low=-1, high=1)
            require(adjusted[0] <= nominal[0] <= nominal[1] <= adjusted[1], "Invalid interval order/nesting")
            if metric == "F1":
                expected = rows[candidate]["scores"]["SwissTrees"] - rows[CELLS[0]]["scores"]["SwissTrees"]
            else:
                key = "precision" if metric == "PPV" else "recall"
                expected = rows[candidate]["endpoint_details"]["SwissTrees"][key] - rows[CELLS[0]]["endpoint_details"]["SwissTrees"][key]
            require(abs(delta-expected) <= 1e-7, "Interval point differs from native decimal endpoint")
            intervals.append(dict(contrast=name, metric=metric, difference_pp=100*delta,
                nominal_low_pp=100*nominal[0], nominal_high_pp=100*nominal[1],
                adjusted_low_pp=100*adjusted[0], adjusted_high_pp=100*adjusted[1]))
    require({r["name"] for r in binding["contrasts"] if r["status"] == "native_records_matched"} == set(CONTRASTS)
            and all(r["metrics"] is None for r in binding["contrasts"] if r["name"] not in CONTRASTS),
            "Unavailable contrasts must remain null")
    return rows, scores, coverage, intervals


def render(rows, intervals, output):
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 9,
        "svg.hashsalt": "native-qfo-three-cell", "svg.fonttype": "none", "pdf.fonttype": 42}):
        fig, axes = plt.subplots(2, 3, figsize=(13.8, 9.8))
        for axis, endpoints, title in ((axes[0, 0], ENDPOINTS[:3], "A  Orthology F1"),
                                       (axes[0, 1], ENDPOINTS[3:], "B  Similarities (not F1)")):
            for offset, cell, color in zip((-.14, 0, .14), CELLS, COLORS):
                axis.scatter([rows[cell]["scores"][e] for e in endpoints],
                    [y+offset for y in range(3)], color=color, s=45, edgecolor="white", linewidth=.5)
            axis.set(xlim=(-.02, 1.02), ylim=(2.45, -.45), xlabel="Native endpoint score")
            axis.set_yticks(range(3), endpoints)
            axis.set_xticks([0, .25, .5, .75, 1])
            axis.grid(axis="x", color="#e9e9e9", linewidth=.7)
            axis.set_title(title, loc="left", fontsize=11, pad=12)
        axis = axes[0, 2]
        values = [100*rows[cell]["relation_coverage"] for cell in CELLS]
        axis.barh(range(3), values, height=.45, color=COLORS)
        for y, value in enumerate(values):
            axis.text(value+1.5, y, f"{value:.2f}%", va="center", fontsize=9)
        axis.set(xlim=(0, 100), ylim=(2.45, -.45), xlabel="Inputs in any submitted pair (%)")
        axis.set_yticks(range(3), [c.upper().replace("_", "") for c in CELLS])
        axis.set_title("C  All-input relation coverage", loc="left", fontsize=11, pad=12)
        axis.grid(axis="x", color="#e9e9e9", linewidth=.7)
        axis.set_axisbelow(True)
        for axis, name, color, title in ((axes[1, 0], CONTRASTS[0], COLORS[2], "D  SwissTrees: candidate expansion"),
            (axes[1, 1], CONTRASTS[1], COLORS[1], "E  SwissTrees: reconciliation")):
            selected = [r for r in intervals if r["contrast"] == name]
            axis.axvline(0, color="#777777", linestyle="--", linewidth=.8)
            for y, row in enumerate(selected):
                axis.plot([row["adjusted_low_pp"], row["adjusted_high_pp"]], [y, y], color="#333333", linewidth=1.5)
                axis.plot([row["nominal_low_pp"], row["nominal_high_pp"]], [y, y], color="#333333", linewidth=4)
                axis.plot(row["difference_pp"], y, "o", color=color, markersize=6)
            bounds = [r[k] for r in selected for k in ("adjusted_low_pp", "adjusted_high_pp")]
            axis.set(xlim=(min(0, *bounds)-3, max(0, *bounds)+3), ylim=(2.45, -.45),
                     xlabel="Difference from P0C0R0 (pp)")
            axis.set_yticks(range(3), ("F1", "Precision", "Recall"))
            axis.set_title(title, loc="left", fontsize=11, pad=12)
            axis.grid(axis="y", color="#e9e9e9", linewidth=.7)
        axis = axes[1, 2]
        for cell, color in zip(CELLS, COLORS):
            detail = rows[cell]["endpoint_details"]["SwissTrees"]
            axis.scatter(detail["recall"], detail["precision"], s=65, color=color, edgecolor="white", linewidth=.5)
        axis.set(xlim=(0, 1.02), ylim=(0, 1.02), xlabel="Recall", ylabel="Precision")
        axis.set_title("F  SwissTrees precision-recall", loc="left", fontsize=11, pad=12)
        axis.grid(color="#e9e9e9", linewidth=.7)
        for axis in axes.flat:
            axis.spines[["top", "right"]].set_visible(False)
            axis.tick_params(length=3)
        fig.suptitle("Native OrthoHMM ablations: three admitted QfO cells", x=.065, ha="left", fontsize=16, y=.98)
        fig.text(.065, .942, "3/7 fresh cells admitted; four scores unavailable. Initial HMM search on; P = 0 in all displayed cells.")
        handles = [plt.Line2D([], [], marker="o", linestyle="none", color=c, markersize=7) for c in COLORS]
        fig.legend(handles, ("P0C0R0: clique pairs", "P0C0R1: inferred pairs (recovered science)",
            "P0C1R0: expanded clique pairs"), loc="upper left", bbox_to_anchor=(.058, .925),
            frameon=False, ncol=3, columnspacing=2)
        fig.text(.065, .125, "D/E: thick = nominal 95% paired-family interval; thin = 42-endpoint adjusted interval. Both adjusted F1 intervals include zero.")
        fig.text(.065, .098, "18 development-exposed SwissTrees families; retained draws reused after exact native family matches. Conditional exchangeability/percentile limits.")
        fig.text(.065, .071, "No intervals for other QfO endpoints. Coverage is not accuracy. Unseeded FAS sampling; the project-defined secondary mean is not plotted.")
        fig.text(.065, .044, "Failed R-on timing remains excluded. No timing comparison, independent confirmation, full-factorial result or OrthoFinder-superiority claim.")
        fig.subplots_adjust(left=.105, right=.98, top=.83, bottom=.23, hspace=.75, wspace=.55)
        for suffix in ("png", "pdf", "svg"):
            fig.savefig(output / f"{STEM}.{suffix}", dpi=200, metadata={"Creator": "OrthoHMM benchmark workflow"})
        plt.close(fig)


def run(snapshot_path, snapshot_sha, binding_path, binding_sha, readback_path, readback_sha, output, validation_python):
    require(not output.exists() and not output.is_symlink(), "Output already exists")
    evidence = []
    snapshot, snapshot_ref = load(snapshot_path, snapshot_sha, evidence)
    binding, binding_ref = load(binding_path, binding_sha, evidence)
    readback, readback_ref = load(readback_path, readback_sha, evidence)
    require(binding["snapshot"] == snapshot_ref and readback["binding"] == binding_ref,
            "Different figure snapshots or readback binding")
    require(readback["source"] == record(Path(__file__).with_name("readback_native_qfo_candidate_swiss.py")),
            "Changed candidate readback source")
    validation = previous.replay_binding(validation_python, snapshot_path, snapshot_sha, binding)
    rows, scores, coverage, intervals = figure_data(snapshot, binding, readback)
    evidence.extend([*binding["evidence"], binding["source"], *readback["checked_inputs"], readback["source"],
        validation["python"], record(previous.__file__), record(__file__)])
    output.mkdir(parents=True, exist_ok=False)
    for name, values in (("scores.tsv", scores), ("coverage.tsv", coverage), ("swiss_intervals.tsv", intervals)):
        with (output / name).open("x", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(values[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(values)
    render(rows, intervals, output)
    for ref in evidence:
        check(ref)
    result = dict(schema="native_qfo_three_cell_figure_v1", snapshot=snapshot_ref, swiss_binding=binding_ref,
        candidate_readback=readback_ref, plotted_cells=list(CELLS), admitted_cells_in_snapshot=3,
        unavailable_score_cells=[r["cell"] for r in snapshot["rows"] if r["accuracy_admitted"] is False],
        plotted_score_endpoints=18, plotted_coverage_endpoints=3, plotted_swiss_contrast_endpoints=6,
        evidence=evidence, outputs=[record(p) for p in sorted(output.iterdir())], source=record(__file__),
        validation=validation, observer_command=sys.orig_argv, python_version=sys.version,
        matplotlib_version=matplotlib.__version__, new_bootstrap_draws=0, new_scoring_or_admission=False,
        scientific_timings_admitted=False, publication_ready=False, visual_review_complete=False)
    with (output / "manifest.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("snapshot", "swiss-binding", "candidate-readback", "output", "validation-python"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("snapshot-sha256", "swiss-binding-sha256", "candidate-readback-sha256"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    result = run(args.snapshot, args.snapshot_sha256, args.swiss_binding, args.swiss_binding_sha256,
                 args.candidate_readback, args.candidate_readback_sha256, args.output, args.validation_python)
    print(json.dumps({k: result[k] for k in ("plotted_cells", "plotted_score_endpoints", "plotted_swiss_contrast_endpoints")}))
