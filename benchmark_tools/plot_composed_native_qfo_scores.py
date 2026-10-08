"""Plot admitted composed-native QfO scores with explicit missing outcomes."""

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
import numpy as np

from benchmark_tools.bootstrap_qfo_factorial import CELLS as FACTORIAL_CELLS, contrasts
from benchmark_tools.export_native_factorial_progress import finite, load, require
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

CELLS = ("p0_c0_r0", "p0_c0_r1", "p0_c1_r0", "p1_c0_r1", "p1_c1_r1")
ENDPOINTS = ("VGNC", "SwissTrees", "TreeFam-A", "GO", "EC", "FAS")
COLORS = ("#12786f", "#b64b37", "#4653a6", "#9a6c00", "#7d438d")
STEM = "composed_native_qfo"
LABELS = ("Orthology F1", "Similarities (not F1)", "All-input relation coverage",
    "Conditional SwissTrees contrasts", "5/7 native cells admitted", "two scores unavailable",
    "Initial HMM search remains on", "42-endpoint adjusted intervals",
    "18 development-exposed SwissTrees families", "Coverage is not accuracy",
    "secondary mean is not plotted", "Failed timing remains excluded", "No isolated timing")


def figure_data(snapshot, uncertainty):
    require(snapshot["schema"] == "composed_native_qfo_scientific_reporting_v1"
        and snapshot["publication_ready"] is False and snapshot["supplied_admissions"] == 5
        and snapshot["historical_four_cell_rows_preserved"] is True,
        "Require the five-cell composed reporting snapshot")
    rows = {row["cell"]: row for row in snapshot["rows"]}
    require(len(snapshot["rows"]) == len(rows) == 7
        and [row["index"] for row in snapshot["rows"]] == list(range(6, 13))
        and {cell for cell, row in rows.items() if row["accuracy_admitted"] is True} == set(CELLS),
        "Changed admitted or missing cell inventory")
    require(rows["p0_c1_r1"]["status"] == "no_supplied_native_admission"
        and rows["p1_c1_r0"]["status"] == "retained_composed_scoring_failure"
        and rows["p1_c1_r0"]["scoring_status"] == "OUT_OF_MEMORY",
        "Missing inference or scoring failure was relabeled")
    require(rows[CELLS[1]]["resources"] is None
        and rows[CELLS[1]]["timing_admitted"] is False
        and rows[CELLS[1]]["timing_eligible"] is False, "Failed timing was repaired")
    require(all(rows[cell]["scientific_timings_admitted"] is False for cell in CELLS[3:]),
        "Allocated/composed figure claims isolated timing admission")
    scores, coverage, status = [], [], []
    for row in snapshot["rows"]:
        admitted = row["accuracy_admitted"]
        require(type(admitted) is bool and set(row["scores"]) == set(ENDPOINTS), "Changed endpoint inventory")
        if admitted:
            for endpoint in ENDPOINTS[:3]:
                detail = row["endpoint_details"][endpoint]
                p, r = (finite(detail[key], key, high=1) for key in ("precision", "recall"))
                value = finite(row["scores"][endpoint], endpoint, high=1)
                require(math.isclose(value, 2*p*r/(p+r) if p+r else 0, rel_tol=0, abs_tol=1e-12),
                    "Precision/recall do not reproduce native F1")
        else:
            require(all(value is None for value in row["scores"].values()) and row["secondary_mean"] is None,
                "Missing scores or mean were imputed")
        for endpoint in ENDPOINTS:
            value = finite(row["scores"][endpoint], endpoint, high=1) if admitted else None
            scores.append(dict(cell=row["cell"], endpoint=endpoint,
                statistic="F1" if endpoint in ENDPOINTS[:3] else "similarity", value=value,
                accuracy_admitted=admitted))
        status.append(dict(index=row["index"], cell=row["cell"], status=row["status"],
            accuracy_admitted=admitted, scoring_status=row.get("scoring_status", "")))
        n, covered, pairs = (row[key] for key in ("input_accessions", "relation_accessions", "submitted_pairs"))
        if n is None:
            require(covered is None and pairs is None and row["relation_coverage"] is None and not admitted,
                "Partially missing coverage")
        else:
            require(type(n) is int and n > 0 and type(covered) is int and 0 <= covered <= n
                and type(pairs) is int and pairs >= 0 and covered <= 2*pairs
                and (pairs == 0 or covered >= 2) and row["relation_coverage"] == covered/n,
                "Coverage count or denominator differs")
        coverage.append(dict(cell=row["cell"], input_accessions=n, relation_accessions=covered,
            submitted_pairs=pairs, relation_coverage=row["relation_coverage"],
            accuracy_admitted=admitted, prediction_semantics=row["prediction_semantics"]))
    require(uncertainty["schema"] == "composed_native_qfo_swiss_bootstrap_v1"
        and uncertainty["status"] == "paired_composed_native_qfo_swiss_intervals"
        and uncertainty["observed_cells"] == list(CELLS)
        and uncertainty["replicates"] == uncertainty["new_bootstrap_draws"] == 100000
        and uncertainty["seed"] == 20260922 and uncertainty["alpha"] == .05
        and uncertainty["multiplicity_endpoints"] == 42 and uncertainty["quantile_method"] == "linear"
        and len(uncertainty["families"]) == len(set(uncertainty["families"])) == 18
        and all(uncertainty[key] is False for key in ("retained_intervals_reused", "unobserved_cells_imputed",
            "independent_confirmation", "new_accuracy_or_resource_admission", "publication_ready")),
        "Changed native uncertainty scope")
    require(set(uncertainty["point_estimates"]) == set(uncertainty["bound_cells"]) == set(CELLS),
        "Changed observed native uncertainty cells")
    for cell in CELLS:
        original, bound = rows[cell], uncertainty["bound_cells"][cell]
        require(bound["admission"] == original["admission"] and bound["index"] == original["index"]
            and bound["native_job_id"] == original["native_job_id"], "Count/admission binding differs")
        for metric, key in (("F1", None), ("PPV", "precision"), ("TPR", "recall")):
            point = finite(uncertainty["point_estimates"][cell][metric], metric, high=1)
            value = original["scores"]["SwissTrees"] if key is None else original["endpoint_details"]["SwissTrees"][key]
            require(abs(point-value) <= 1e-7, "Bootstrap point differs from admitted native endpoint")
    definitions = contrasts()
    require(len(uncertainty["comparisons"]) == len(definitions), "Changed contrast inventory")
    intervals = []
    for effect, definition in zip(uncertainty["comparisons"], definitions):
        require(all(effect[key] == value for key, value in definition.items()), "Changed contrast definition")
        needed = [cell for cell, weight in zip(FACTORIAL_CELLS, definition["weights"]) if weight]
        absent = [cell for cell in needed if cell not in CELLS]
        require(effect["missing_cells"] == absent, "Changed missing contrast cells")
        if absent:
            require(effect["status"] == "native_counts_unavailable" and effect["metrics"] is None
                and effect["family_differences"] is None, "Unavailable interval was imputed")
            continue
        require(effect["status"] == "native_counts_resampled" and set(effect["metrics"]) == {"F1", "PPV", "TPR"},
            "Changed available interval scope")
        for metric in ("F1", "PPV", "TPR"):
            value = effect["metrics"][metric]
            difference = finite(value["difference"], "difference", low=-2, high=2)
            nominal, adjusted = value["paired_percentile_ci"], value["bonferroni_percentile_ci"]
            require(len(nominal) == len(adjusted) == 2, "Changed interval shape")
            for limit in [*nominal, *adjusted]: finite(limit, "interval", low=-2, high=2)
            require(adjusted[0] <= nominal[0] <= nominal[1] <= adjusted[1], "Invalid interval nesting")
            expected = sum(weight*uncertainty["point_estimates"][cell][metric]
                for cell, weight in zip(FACTORIAL_CELLS, definition["weights"]) if weight)
            require(abs(difference-expected) <= 1e-12, "Contrast point differs from native arithmetic")
            intervals.append(dict(contrast=definition["name"], metric=metric, difference_pp=100*difference,
                nominal_low_pp=100*nominal[0], nominal_high_pp=100*nominal[1],
                adjusted_low_pp=100*adjusted[0], adjusted_high_pp=100*adjusted[1]))
    return rows, scores, coverage, status, intervals


def replay_snapshot(python, path, digest):
    script = ("import json,sys;sys.path.insert(0,sys.argv[3]);"
        "from benchmark_tools.audit_composed_native_qfo_swiss_counts import replay_snapshot;"
        "from benchmark_tools.prepare_ob_candidate_neighborhood import check;"
        "e=[];s,r=replay_snapshot(sys.argv[1],sys.argv[2],e);[check(x) for x in e];"
        "print(json.dumps(dict(exact_snapshot_replay=True,admitted_cells=s['supplied_admissions'])))")
    command = [str(python.absolute()), "-I", "-B", "-c", script, str(path.absolute()), digest,
        str(Path(__file__).resolve().parent.parent)]
    env = {key: value for key, value in os.environ.items() if key not in (
        "PYTHONPATH", "PYTHONHOME", "PYTHONUSERBASE", "LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT")}
    env.update(PYTHONNOUSERSITE="1", PYTHONDONTWRITEBYTECODE="1", PYTHONHASHSEED="0",
        OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    observed = subprocess.run(command, env=env, capture_output=True, text=True, timeout=300, check=False)
    require(observed.returncode == 0, "Scientific replay failed: " + observed.stderr[-4000:])
    value = json.loads(observed.stdout)
    require(value == dict(exact_snapshot_replay=True, admitted_cells=5), "Changed scientific replay scope")
    return dict(command=command, returncode=observed.returncode, result=value, python=record(python),
        replay_helper=record(Path(__file__).with_name("audit_composed_native_qfo_swiss_counts.py")))


def render(rows, intervals, output):
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 9, "svg.fonttype": "none",
        "pdf.fonttype": 42, "svg.hashsalt": STEM}):
        fig = plt.figure(figsize=(15, 12))
        grid = fig.add_gridspec(2, 3, height_ratios=(1, 2.2), hspace=.35, wspace=.45)
        for column, endpoints, title in ((0, ENDPOINTS[:3], "A  Orthology F1"),
            (1, ENDPOINTS[3:], "B  Similarities (not F1)")):
            axis = fig.add_subplot(grid[0, column])
            for offset, cell, color in zip(np.linspace(-.22, .22, 5), CELLS, COLORS):
                axis.scatter([rows[cell]["scores"][e] for e in endpoints],
                    [y+offset for y in range(3)], color=color, s=40, edgecolor="white", linewidth=.5)
            axis.set(xlim=(-.02, 1.02), ylim=(2.5, -.5), xlabel="Native endpoint score")
            axis.set_yticks(range(3), endpoints)
            axis.set_title(title, loc="left", fontsize=11)
            axis.grid(axis="x", color="#e9e9e9")
        axis = fig.add_subplot(grid[0, 2])
        values = [100*rows[cell]["relation_coverage"] for cell in CELLS]
        axis.barh(range(5), values, color=COLORS, height=.5)
        for y, value in enumerate(values): axis.text(value+1.5, y, f"{value:.2f}%", va="center")
        axis.set(xlim=(0, 100), ylim=(4.6, -.6), xlabel="Inputs in any submitted pair (%)")
        axis.set_yticks(range(5), [cell.upper().replace("_", "") for cell in CELLS])
        axis.set_title("C  All-input relation coverage", loc="left", fontsize=11)
        axis = fig.add_subplot(grid[1, :])
        axis.axvline(0, color="#777777", linewidth=.8, linestyle="--")
        for y, row in enumerate(intervals):
            axis.plot([row["adjusted_low_pp"], row["adjusted_high_pp"]], [y, y], color="#444444", linewidth=1.5)
            axis.plot([row["nominal_low_pp"], row["nominal_high_pp"]], [y, y], color="#444444", linewidth=4)
            axis.plot(row["difference_pp"], y, "o", color="#202020", markersize=6)
        axis.set_yticks(range(len(intervals)), [row["contrast"]+" / "+row["metric"] for row in intervals])
        axis.set(ylim=(len(intervals)-.5, -.5), xlabel="Native metric difference (percentage points)")
        axis.set_title("D  Conditional SwissTrees contrasts", loc="left", fontsize=11)
        axis.grid(axis="y", color="#e9e9e9")
        for axis in fig.axes: axis.spines[["top", "right"]].set_visible(False)
        fig.suptitle("Native OrthoHMM QfO ablations", x=.06, ha="left", fontsize=16, y=.98)
        fig.text(.06, .95, "5/7 native cells admitted; two scores unavailable. Initial HMM search remains on in every cell.")
        handles = [plt.Line2D([], [], marker="o", linestyle="none", color=color) for color in COLORS]
        fig.legend(handles, [cell.upper().replace("_", "") for cell in CELLS], loc="upper left",
            bbox_to_anchor=(.052, .935), frameon=False, ncol=5)
        fig.text(.06, .895, "P: profile refinement; C: candidate expansion; R: phylogenetic pair inference (0 off, 1 on).")
        footers = ("Thick = nominal 95% paired-family intervals; thin = 42-endpoint adjusted intervals. No assumed direction or zero inclusion.",
            "18 development-exposed SwissTrees families; new native-count draws. Conditional percentile/exchangeability limitations remain.",
            "Coverage is not accuracy. Failed11 scoring retains coverage but no admitted score; see status and coverage tables.",
            "No intervals for other endpoints; the project-defined secondary mean is not plotted. Initial HMM search is retained.",
            "Failed timing remains excluded. No isolated timing, independent confirmation, full-factorial or superiority claim.")
        for y, text in zip((.135, .111, .087, .063, .039), footers): fig.text(.06, y, text)
        fig.subplots_adjust(left=.22, right=.98, top=.85, bottom=.21)
        for suffix in ("png", "pdf", "svg"):
            fig.savefig(output / f"{STEM}.{suffix}", dpi=200,
                metadata={"Creator": "OrthoHMM benchmark workflow"})
        plt.close(fig)


def run(snapshot_path, snapshot_sha, intervals_path, intervals_sha, output, validation_python):
    require(not output.exists() and not output.is_symlink(), "Output already exists")
    evidence = []
    snapshot, snapshot_ref = load(snapshot_path, snapshot_sha, evidence)
    uncertainty, uncertainty_ref = load(intervals_path, intervals_sha, evidence)
    require(snapshot["source"] == record(Path(__file__).with_name("export_composed_native_qfo_scientific_scores.py"))
        and uncertainty["source"] == record(Path(__file__).with_name("bootstrap_composed_native_qfo_swiss.py"))
        and uncertainty["snapshot"] == snapshot_ref, "Changed source or cross-snapshot intervals")
    rows, scores, coverage, status, intervals = figure_data(snapshot, uncertainty)
    validation = replay_snapshot(validation_python, snapshot_path, snapshot_sha)
    evidence.extend([*snapshot["evidence"], *snapshot["outputs"], snapshot["source"],
        *uncertainty["evidence"], *uncertainty["helpers"], uncertainty["source"],
        validation["python"], validation["replay_helper"], record(__file__)])
    for ref in evidence: check(ref)
    output.mkdir(parents=True, exist_ok=False)
    for name, values in (("scores.tsv", scores), ("coverage.tsv", coverage),
        ("status.tsv", status), ("swiss_intervals.tsv", intervals)):
        with (output / name).open("x", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(values[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(values)
    render(rows, intervals, output)
    for ref in evidence: check(ref)
    result = dict(schema="composed_native_qfo_figure_v1", source=record(__file__), snapshot=snapshot_ref,
        native_intervals=uncertainty_ref, plotted_cells=list(CELLS), unavailable_score_cells=[r["cell"]
            for r in snapshot["rows"] if not r["accuracy_admitted"]], plotted_score_endpoints=30,
        plotted_coverage_endpoints=5, plotted_swiss_contrast_endpoints=len(intervals),
        table_score_rows=len(scores), table_coverage_rows=len(coverage), table_status_rows=len(status),
        evidence=evidence, outputs=[record(path) for path in sorted(output.iterdir())], validation=validation,
        python_version=sys.version, matplotlib_version=matplotlib.__version__,
        new_bootstrap_draws=0, consumed_native_bootstrap_draws=100000,
        new_scoring_or_admission=False, scientific_timings_admitted=False,
        independent_confirmation=False, publication_ready=False, visual_review_complete=False,
        timing_disclosure=snapshot["timing_disclosure"])
    with (output / "manifest.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("snapshot", "native-intervals", "output", "validation-python"):
        parser.add_argument("--"+name, type=Path, required=True)
    for name in ("snapshot-sha256", "native-intervals-sha256"):
        parser.add_argument("--"+name, required=True)
    args = parser.parse_args()
    result = run(args.snapshot, args.snapshot_sha256, args.native_intervals, args.native_intervals_sha256,
        args.output, args.validation_python)
    print(json.dumps(dict(plotted_cells=result["plotted_cells"], scores=result["plotted_score_endpoints"])))


if __name__ == "__main__":
    main()
