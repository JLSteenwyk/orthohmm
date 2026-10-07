"""Plot four admitted native QfO cells and three conditional SwissTrees contrasts."""

import argparse
import csv
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from benchmark_tools.export_native_factorial_progress import finite, load, require
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

CELLS = ("p0_c0_r0", "p0_c0_r1", "p0_c1_r0", "p1_c0_r1")
CONTRASTS = (("C_at_P0_R0", CELLS[2], CELLS[0]),
             ("R_at_P0_C0", CELLS[1], CELLS[0]),
             ("P_at_C0_R1", CELLS[3], CELLS[1]))
ENDPOINTS = ("VGNC", "SwissTrees", "TreeFam-A", "GO", "EC", "FAS")
COLORS = ("#12786f", "#b64b37", "#4653a6", "#9a6c00")
STEM = "native_qfo_four_cell"


def figure_data(snapshot, binding, readback):
    require(snapshot["schema"] == "allocated_native_qfo_scientific_reporting_snapshot_v1"
            and snapshot["publication_ready"] is False and len(snapshot["rows"]) == 7
            and snapshot["supplied_admissions"] == 4 and snapshot["supplied_allocated_admissions"] == 1,
            "Require four-cell allocated reporting snapshot")
    rows = {row["cell"]: row for row in snapshot["rows"]}
    require(len(rows) == 7 and set(CELLS) <= set(rows)
            and sum(row["accuracy_admitted"] is True for row in rows.values()) == 4,
            "Require exactly four distinct admitted cells")
    scores, coverage = [], []
    for cell, status in zip(CELLS, ("supplied_native_admission", "supplied_recovered_scientific_admission",
                                  "supplied_native_admission", "supplied_allocated_native_admission")):
        row = rows[cell]
        require(row["accuracy_admitted"] is True and row["status"] == status
                and set(row["scores"]) == set(ENDPOINTS), "Changed plotted admission or endpoint scope")
        if cell == CELLS[1]:
            require(row["resources"] is None and row["timing_admitted"] is False
                    and row["timing_eligible"] is False, "Recovered plot repairs failed timing")
        if cell == CELLS[3]:
            require(row["scientific_timings_admitted"] is False, "Allocated plot claims isolated timing admission")
        for endpoint in ENDPOINTS:
            scores.append(dict(cell=cell, endpoint=endpoint, statistic="F1" if endpoint in ENDPOINTS[:3]
                               else "similarity", value=finite(row["scores"][endpoint], endpoint, high=1)))
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
            recall=row["endpoint_details"]["SwissTrees"]["recall"], prediction_semantics=row["prediction_semantics"]))
    require(all(row["accuracy_admitted"] is False and all(v is None for v in row["scores"].values())
                for cell, row in rows.items() if cell not in CELLS), "Unavailable scores must remain null")
    require(binding["schema"] == "allocated_native_qfo_retained_swiss_uncertainty_binding_v1"
            and all(binding[k] is False for k in ("publication_ready", "independent_confirmation",
                "new_accuracy_or_resource_admission")) and binding["multiplicity_endpoints"] == 42
            and binding["replicates_reused"] == 100000 and binding["seed_reused"] == 20260922
            and type(binding["new_bootstrap_draws"]) is int and binding["new_bootstrap_draws"] == 0,
            "Changed uncertainty scope")
    require(readback["schema"] == "allocated_native_qfo_profile_swiss_rational_readback_v1"
            and readback["cells"] == [CELLS[1], CELLS[3]] and readback["families_checked"] == 18
            and readback["profile_pair_labels_matched"] == 10765
            and readback["prior_matched_contrasts_checked"] == 2
            and readback["prior_matched_contrasts_unchanged"] is True
            and all(readback[k] is False for k in ("publication_ready", "independent_confirmation",
                "new_accuracy_or_resource_admission", "scientific_timings_admitted"))
            and type(readback["new_bootstrap_draws"]) is int and readback["new_bootstrap_draws"] == 0,
            "Changed profile readback scope")
    effects = {row["name"]: row for row in binding["contrasts"]}
    require(len(effects) == len(binding["contrasts"]) == 14
            and {name for name, row in effects.items() if row["status"] == "native_records_matched"}
                == {name for name, _, _ in CONTRASTS}, "Changed available contrast inventory")
    require(all(row["metrics"] is None for row in effects.values() if row["status"] != "native_records_matched"),
            "Unavailable contrasts must remain null")
    require(readback["contrast"] == effects["P_at_C0_R1"], "Readback profile contrast differs")
    intervals = []
    for name, candidate, reference in CONTRASTS:
        effect = effects[name]
        require(not effect["missing_cells"] and not effect["differing_cells"] and effect["candidate"] == candidate
                and effect["reference"] == reference and set(effect["metrics"]) == {"F1", "PPV", "TPR"},
                "Changed conditional contrast")
        for metric in ("F1", "PPV", "TPR"):
            value = effect["metrics"][metric]
            delta = finite(value["difference"], "difference", low=-1, high=1)
            nominal, adjusted = value["paired_percentile_ci"], value["bonferroni_percentile_ci"]
            require(len(nominal) == len(adjusted) == 2, "Changed interval shape")
            for bound in [*nominal, *adjusted]:
                finite(bound, "interval", low=-1, high=1)
            require(adjusted[0] <= nominal[0] <= nominal[1] <= adjusted[1], "Invalid interval order/nesting")
            if metric == "F1":
                expected = rows[candidate]["scores"]["SwissTrees"] - rows[reference]["scores"]["SwissTrees"]
            else:
                key = "precision" if metric == "PPV" else "recall"
                expected = rows[candidate]["endpoint_details"]["SwissTrees"][key] - rows[reference]["endpoint_details"]["SwissTrees"][key]
            require(abs(delta - expected) <= 1e-7, "Interval point differs from native decimal endpoint")
            intervals.append(dict(contrast=name, metric=metric, difference_pp=100*delta,
                nominal_low_pp=100*nominal[0], nominal_high_pp=100*nominal[1],
                adjusted_low_pp=100*adjusted[0], adjusted_high_pp=100*adjusted[1]))
    return rows, scores, coverage, intervals


def replay_snapshot(python, path, digest):
    root = Path(__file__).resolve().parent.parent
    script = (
        "import json,sys;sys.path.insert(0,sys.argv[3]);"
        "from benchmark_tools.audit_allocated_native_qfo_swiss_counts import replay_snapshot;"
        "from benchmark_tools.prepare_ob_candidate_neighborhood import check;"
        "e=[];s,r=replay_snapshot(sys.argv[1],sys.argv[2],e);"
        "[check(ref) for ref in e];"
        "print(json.dumps(dict(exact_snapshot_replay=True,admitted_cells=s['supplied_admissions'],"
        "python_version=sys.version)))"
    )
    command = [str(python.absolute()), "-I", "-B", "-c", script, str(path.absolute()), digest, str(root)]
    env = {k: v for k, v in os.environ.items() if k not in (
        "PYTHONPATH", "PYTHONHOME", "PYTHONUSERBASE", "LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT")}
    env.update(PYTHONNOUSERSITE="1", PYTHONDONTWRITEBYTECODE="1", PYTHONHASHSEED="0",
               OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    result = subprocess.run(command, env=env, capture_output=True, text=True, timeout=300, check=False)
    require(result.returncode == 0, "Scientific reporting replay failed: " + result.stderr[-4000:])
    value = json.loads(result.stdout)
    require(value["exact_snapshot_replay"] is True and value["admitted_cells"] == 4, "Wrong reporting replay scope")
    return dict(command=command, returncode=result.returncode, stdout=result.stdout, stderr=result.stderr,
        python=record(python), result=value,
        replay_helper=record(Path(__file__).with_name("audit_allocated_native_qfo_swiss_counts.py")))


def render(rows, intervals, output):
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 9,
        "svg.hashsalt": "native-qfo-four-cell", "svg.fonttype": "none", "pdf.fonttype": 42}):
        fig, axes = plt.subplots(2, 3, figsize=(14.8, 10.6))
        for axis, endpoints, title in ((axes[0, 0], ENDPOINTS[:3], "A  Orthology F1"),
                                      (axes[0, 1], ENDPOINTS[3:], "B  Similarities (not F1)")):
            for offset, cell, color in zip((-.18, -.06, .06, .18), CELLS, COLORS):
                axis.scatter([rows[cell]["scores"][e] for e in endpoints],
                    [y+offset for y in range(3)], color=color, s=42, edgecolor="white", linewidth=.5)
            axis.set(xlim=(-.02, 1.02), ylim=(2.45, -.45), xlabel="Native endpoint score")
            axis.set_yticks(range(3), endpoints)
            axis.set_xticks([0, .25, .5, .75, 1])
            axis.grid(axis="x", color="#e9e9e9", linewidth=.7)
            axis.set_title(title, loc="left", fontsize=11, pad=12)
        axis = axes[0, 2]
        values = [100*rows[cell]["relation_coverage"] for cell in CELLS]
        axis.barh(range(4), values, height=.48, color=COLORS)
        for y, value in enumerate(values):
            axis.text(value+1.5, y, f"{value:.2f}%", va="center")
        axis.set(xlim=(0, 100), ylim=(3.55, -.55), xlabel="Inputs in any submitted pair (%)")
        axis.set_yticks(range(4), [c.upper().replace("_", "") for c in CELLS])
        axis.set_title("C  All-input relation coverage", loc="left", fontsize=11, pad=12)
        axis.grid(axis="x", color="#e9e9e9", linewidth=.7)
        axis.set_axisbelow(True)
        for axis, (name, candidate, reference), color, title in zip(axes[1], CONTRASTS,
            (COLORS[2], COLORS[1], COLORS[3]), ("D  Candidate expansion (P0/R0)",
            "E  Reconciliation (P0/C0)", "F  Profile refinement (C0/R1)")):
            selected = [row for row in intervals if row["contrast"] == name]
            axis.axvline(0, color="#777777", linestyle="--", linewidth=.8)
            for y, row in enumerate(selected):
                axis.plot([row["adjusted_low_pp"], row["adjusted_high_pp"]], [y, y], color="#333333", linewidth=1.5)
                axis.plot([row["nominal_low_pp"], row["nominal_high_pp"]], [y, y], color="#333333", linewidth=4)
                axis.plot(row["difference_pp"], y, "o", color=color, markersize=6)
            bounds = [row[key] for row in selected for key in ("adjusted_low_pp", "adjusted_high_pp")]
            padding = max(.5, (max(0, *bounds) - min(0, *bounds)) * .08)
            axis.set(xlim=(min(0, *bounds)-padding, max(0, *bounds)+padding), ylim=(2.45, -.45),
                xlabel=candidate.upper().replace("_", "") + " - " + reference.upper().replace("_", "") + " (pp)")
            axis.set_yticks(range(3), ("F1", "Precision", "Recall"))
            axis.set_title(title, loc="left", fontsize=11, pad=12)
            axis.grid(axis="y", color="#e9e9e9", linewidth=.7)
        for axis in axes.flat:
            axis.spines[["top", "right"]].set_visible(False)
            axis.tick_params(length=3)
        fig.suptitle("Native OrthoHMM ablations: four admitted QfO cells", x=.065, ha="left", fontsize=16, y=.98)
        fig.text(.065, .944, "4/7 fresh cells admitted; three scores unavailable. Initial HMM search remains on in every cell.")
        handles = [plt.Line2D([], [], marker="o", linestyle="none", color=color, markersize=7) for color in COLORS]
        fig.legend(handles, ("P0C0R0: clique pairs", "P0C0R1: inferred pairs (recovered science)",
            "P0C1R0: expanded clique pairs", "P1C0R1: profile-refined inferred pairs"),
            loc="upper left", bbox_to_anchor=(.058, .922), frameon=False, ncol=2, columnspacing=3)
        fig.text(.065, .126, "D-F: thick = nominal 95% paired-family interval; thin = 42-endpoint adjusted interval. All three adjusted F1 intervals include zero.")
        fig.text(.065, .099, "18 development-exposed SwissTrees families; retained draws reused after exact native family matches. Conditional percentile/exchangeability limits.")
        fig.text(.065, .072, "No intervals for other endpoints. Coverage is not accuracy. FAS sampling is unseeded; the project-defined secondary mean is not plotted.")
        fig.text(.065, .045, "Failed reference timing remains excluded. No isolated timing, independent confirmation, full-factorial or OrthoFinder-superiority claim.")
        fig.subplots_adjust(left=.105, right=.98, top=.80, bottom=.23, hspace=.80, wspace=.57)
        for suffix in ("png", "pdf", "svg"):
            fig.savefig(output / f"{STEM}.{suffix}", dpi=200, metadata={"Creator": "OrthoHMM benchmark workflow"})
        plt.close(fig)


def run(snapshot_path, snapshot_sha, binding_path, binding_sha, readback_path, readback_sha, output, validation_python):
    require(not output.exists() and not output.is_symlink(), "Output already exists")
    evidence = []
    snapshot, snapshot_ref = load(snapshot_path, snapshot_sha, evidence)
    binding, binding_ref = load(binding_path, binding_sha, evidence)
    readback, readback_ref = load(readback_path, readback_sha, evidence)
    require(snapshot["source"] == record(Path(__file__).with_name("export_allocated_native_qfo_scientific_scores.py"))
        and binding["source"] == record(Path(__file__).with_name("bind_allocated_native_qfo_swiss_uncertainty.py"))
        and readback["source"] == record(Path(__file__).with_name("readback_allocated_native_qfo_profile_swiss.py")),
        "Changed reporting/binding/readback source")
    require(binding["snapshot"] == snapshot_ref and readback["binding"] == binding_ref, "Different figure bindings")
    rows, scores, coverage, intervals = figure_data(snapshot, binding, readback)
    validation = replay_snapshot(validation_python, snapshot_path, snapshot_sha)
    evidence.extend([*snapshot["evidence"], *snapshot["outputs"], snapshot["source"], *binding["evidence"],
        *binding["helpers"], binding["source"], *readback["checked_inputs"], readback["source"],
        validation["python"], validation["replay_helper"], record(__file__)])
    for ref in evidence:
        check(ref)
    output.mkdir(parents=True, exist_ok=False)
    for name, values in (("scores.tsv", scores), ("coverage.tsv", coverage), ("swiss_intervals.tsv", intervals)):
        with (output / name).open("x", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(values[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(values)
    render(rows, intervals, output)
    for ref in evidence:
        check(ref)
    result = dict(schema="native_qfo_four_cell_figure_v1", source=record(__file__), snapshot=snapshot_ref,
        swiss_binding=binding_ref, profile_readback=readback_ref, plotted_cells=list(CELLS),
        admitted_cells_in_snapshot=4, unavailable_score_cells=[r["cell"] for r in snapshot["rows"]
            if r["accuracy_admitted"] is False], plotted_score_endpoints=24, plotted_coverage_endpoints=4,
        plotted_swiss_contrast_endpoints=9, evidence=evidence,
        outputs=[record(path) for path in sorted(output.iterdir())], validation=validation,
        observer_command=sys.orig_argv, python_version=sys.version, matplotlib_version=matplotlib.__version__,
        new_bootstrap_draws=0, new_scoring_or_admission=False, scientific_timings_admitted=False,
        independent_confirmation=False, publication_ready=False, visual_review_complete=False)
    with (output / "manifest.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("snapshot", "swiss-binding", "profile-readback", "output", "validation-python"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("snapshot-sha256", "swiss-binding-sha256", "profile-readback-sha256"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    result = run(args.snapshot, args.snapshot_sha256, args.swiss_binding, args.swiss_binding_sha256,
                 args.profile_readback, args.profile_readback_sha256, args.output, args.validation_python)
    print(json.dumps({key: result[key] for key in ("plotted_cells", "plotted_score_endpoints", "plotted_swiss_contrast_endpoints")}))


if __name__ == "__main__":
    main()
