"""Plot admitted P0/C0 QfO ablations and guarded retained SwissTrees intervals."""

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools import export_native_qfo_scientific_scores as reporter
from benchmark_tools.export_native_factorial_progress import finite, load, require
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

CELLS = ("p0_c0_r0", "p0_c0_r1")
COLORS = ("#12786f", "#b64b37")
VALIDATION_PYTHON_SHA = "8b1cd756be711ef53f35cb6c954472fdfc52094c4619a01f553f64354587388b"
VALIDATION_VERSIONS = dict(python="3.10.13", biopython="1.87", numpy="2.2.6", psutil="7.2.2")
WORKER = """import json, sys
from importlib.metadata import version
from pathlib import Path
from benchmark_tools import bind_native_qfo_swiss_uncertainty as binder
request = json.loads(sys.argv[1])
versions = dict(python=".".join(map(str, sys.version_info[:3])),
                **{name: version(name) for name in ("biopython", "numpy", "psutil")})
if versions != request["versions"]:
    raise ValueError("Validation environment versions differ")
result = binder.bind(Path(request["snapshot"]), request["snapshot_sha"], request["audits"],
                     request["retained_counts"], request["bootstrap"])
print(json.dumps(dict(binding=result, versions=versions, executable=sys.executable,
                     prefix=sys.prefix), sort_keys=True, allow_nan=False))
"""


def replay_binding(validation_python, snapshot_path, snapshot_sha, binding):
    # Keep the venv entry point: resolving its symlink selects a different environment.
    invocation = str(validation_python.absolute())
    python_ref = record(invocation)
    require(python_ref["sha256"] == VALIDATION_PYTHON_SHA, "Validation Python binary differs")
    audits = {r["count_audit"]["path"]: r["count_audit"]["sha256"] for r in binding["bound_cells"].values()}
    request = dict(snapshot=str(snapshot_path.absolute()), snapshot_sha=snapshot_sha,
        audits=list(audits.items()), retained_counts=binding["retained_counts"]["path"],
        bootstrap=binding["bootstrap"]["path"], versions=VALIDATION_VERSIONS)
    environment = dict(os.environ)
    for name in ("PYTHONHOME", "PYTHONPATH", "LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT"):
        environment.pop(name, None)
    environment.update(PYTHONNOUSERSITE="1", PYTHONDONTWRITEBYTECODE="1", PYTHONHASHSEED="0",
                       OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    command = [invocation, "-B", "-c", WORKER, json.dumps(request, sort_keys=True)]
    process = subprocess.run(command, cwd=Path(__file__).resolve().parent.parent,
        env=environment, capture_output=True, text=True, timeout=180, check=False)
    require(process.returncode == 0, "Validation worker failed: " + process.stderr)
    require(not process.stderr, "Validation worker emitted unexpected stderr")
    response = json.loads(process.stdout)
    require(response["versions"] == VALIDATION_VERSIONS
            and response["executable"] == invocation
            and response["prefix"] == str(validation_python.absolute().parent.parent),
            "Validation environment identity differs")
    require(response["binding"] == binding, "Retained native uncertainty binding replay differs")
    check(python_ref)
    return dict(python=python_ref, invocation=invocation, command=command,
        versions=response["versions"], prefix=response["prefix"],
        stdout_sha256=hashlib.sha256(process.stdout.encode()).hexdigest(),
        exact_binding_replay=True, new_scoring_or_admission=False)


def figure_data(snapshot, binding):
    require(snapshot["schema"] == "native_qfo_scientific_reporting_snapshot_v1"
            and snapshot["publication_ready"] is False and len(snapshot["rows"]) == 7,
            "Require incomplete-publication native scientific snapshot")
    rows = {row["cell"]: row for row in snapshot["rows"]}
    require(len(rows) == 7 and set(CELLS) <= set(rows), "Missing or duplicate native cells")
    admitted = sum(row["accuracy_admitted"] is True for row in rows.values())
    scores = []
    for cell in CELLS:
        row = rows[cell]
        require(row["accuracy_admitted"] is True and row["status"] in (
            "supplied_native_admission", "supplied_recovered_scientific_admission"), "Unadmitted plotted cell")
        if row["status"] == "supplied_recovered_scientific_admission":
            require(row["resources"] is None and row["timing_admitted"] is False
                    and row["timing_eligible"] is False, "Recovered plot relabels failed timing")
        require(set(row["scores"]) == set(reporter.original.ENDPOINTS), "Incomplete plotted endpoints")
        for endpoint, value in row["scores"].items():
            scores.append(dict(cell=cell, endpoint=endpoint, statistic="F1" if endpoint in
                reporter.original.F1_ENDPOINTS else "similarity", value=finite(value, endpoint, high=1)))
        detail = row["endpoint_details"]["SwissTrees"]
        p, r = (finite(detail[name], name, high=1) for name in ("precision", "recall"))
        require(math.isclose(row["scores"]["SwissTrees"], reporter.original.harmonic_mean(p, r),
                             rel_tol=0, abs_tol=1e-12), "Plotted precision/recall do not reproduce F1")
    require(binding["schema"] == "native_qfo_retained_swiss_uncertainty_binding_v1"
            and binding["snapshot"]["sha256"]
            and binding["publication_ready"] is False and binding["new_bootstrap_draws"] == 0
            and binding["multiplicity_endpoints"] == 42 and binding["replicates_reused"] == 100000,
            "Changed retained uncertainty scope")
    matches = [row for row in binding["contrasts"] if row["name"] == "R_at_P0_C0"]
    require(len(matches) == 1 and matches[0]["status"] == "native_records_matched"
            and not matches[0]["missing_cells"] and not matches[0]["differing_cells"]
            and matches[0]["candidate"] == CELLS[1] and matches[0]["reference"] == CELLS[0],
            "Require matched native R contrast")
    effect, intervals = matches[0], []
    require(set(effect["metrics"]) == {"F1", "PPV", "TPR"}, "Incomplete SwissTrees contrast")
    for metric in ("F1", "PPV", "TPR"):
        item = effect["metrics"][metric]
        d = finite(item["difference"], metric + " difference", low=-1, high=1)
        nominal, adjusted = item["paired_percentile_ci"], item["bonferroni_percentile_ci"]
        require(len(nominal) == len(adjusted) == 2, "Invalid interval shape")
        for v in [*nominal, *adjusted]:
            finite(v, metric + " interval", low=-1, high=1)
        require(adjusted[0] <= nominal[0] <= nominal[1] <= adjusted[1], "Invalid or nonnested intervals")
        if metric == "F1":
            expected = rows[CELLS[1]]["scores"]["SwissTrees"] - rows[CELLS[0]]["scores"]["SwissTrees"]
        else:
            key = "precision" if metric == "PPV" else "recall"
            expected = rows[CELLS[1]]["endpoint_details"]["SwissTrees"][key] - rows[CELLS[0]]["endpoint_details"]["SwissTrees"][key]
        require(abs(d - expected) <= 1e-7, "Interval point differs from native decimal endpoints")
        intervals.append(dict(metric=metric, difference_pp=100*d, nominal_low_pp=100*nominal[0],
            nominal_high_pp=100*nominal[1], adjusted_low_pp=100*adjusted[0], adjusted_high_pp=100*adjusted[1]))
    return rows, scores, intervals, admitted


def render(rows, intervals, admitted, output):
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 10,
        "svg.hashsalt": "native-qfo-p0c0", "svg.fonttype": "none", "pdf.fonttype": 42}):
        fig, axes = plt.subplots(2, 2, figsize=(12.5, 8))
        for axis, endpoints, title in ((axes[0, 0], ("VGNC", "SwissTrees", "TreeFam-A"), "A  Orthology F1"),
                                       (axes[0, 1], ("GO", "EC", "FAS"), "B  Similarity endpoints (not F1)")):
            for y, endpoint in enumerate(endpoints):
                values = [rows[cell]["scores"][endpoint] for cell in CELLS]
                axis.plot(values, [y, y], color="#a8a8a8", linewidth=1.5, zorder=1)
                for value, color in zip(values, COLORS):
                    axis.scatter(value, y, color=color, s=65, edgecolor="white", linewidth=.7, zorder=2)
            axis.set(title=title, xlim=(-.02, 1.02), ylim=(2.5, -.5), xlabel="Native score")
            axis.set_yticks(range(3), endpoints)
            axis.set_xticks([0, .25, .5, .75, 1])
            axis.grid(axis="x", color="#e9e9e9", linewidth=.7)
        axis = axes[1, 0]
        for cell, color in zip(CELLS, COLORS):
            detail = rows[cell]["endpoint_details"]["SwissTrees"]
            axis.scatter(detail["recall"], detail["precision"], color=color, s=90, edgecolor="white", linewidth=.7)
        axis.set(title="C  SwissTrees precision-recall", xlim=(0, 1.02), ylim=(0, 1.02),
                 xlabel="Recall", ylabel="Precision")
        axis.grid(color="#e9e9e9", linewidth=.7)
        axis = axes[1, 1]
        axis.axvline(0, color="#777777", linestyle="--", linewidth=.8)
        for y, row in enumerate(intervals):
            axis.plot([row["adjusted_low_pp"], row["adjusted_high_pp"]], [y, y], color="#333333", linewidth=1.5)
            axis.plot([row["nominal_low_pp"], row["nominal_high_pp"]], [y, y], color="#333333", linewidth=4)
            axis.plot(row["difference_pp"], y, "o", color=COLORS[1], markersize=6)
        bounds = [r[k] for r in intervals for k in ("difference_pp", "adjusted_low_pp", "adjusted_high_pp")]
        axis.set(title="D  SwissTrees: R-on minus R-off", xlim=(min(0, *bounds)-5, max(0, *bounds)+5),
                 ylim=(2.5, -.5), xlabel="Difference (percentage points)")
        axis.set_yticks(range(3), ("F1", "Precision", "Recall"))
        axis.grid(axis="y", color="#e9e9e9", linewidth=.7)
        for axis in axes.flat:
            axis.spines[["top", "right"]].set_visible(False)
            axis.tick_params(length=3)
            axis.set_title(axis.get_title(), loc="left", fontsize=11, pad=12)
            axis.set_title("")
        fig.suptitle("Native OrthoHMM ablation: QfO development set", x=.065, ha="left", fontsize=16, y=.98)
        fig.text(.065, .937, f"Initial HMM search on; P = 0, C = 0. Showing two cells; {admitted}/7 fresh cells admitted.")
        handles = [plt.Line2D([], [], marker="o", linestyle="none", color=color, markersize=7) for color in COLORS]
        fig.legend(handles, ("R-off: group-clique pairs", "R-on: inferred pairs"), loc="upper left",
                   bbox_to_anchor=(.055, .915), frameon=False, ncol=2)
        fig.text(.065, .065, "D: thick = nominal 95% interval; thin = 42-endpoint adjusted interval. Retained draws reused only after exact family matches.", fontsize=9)
        fig.text(.065, .040, "18 development-exposed SwissTrees families; conditional exchangeability and percentile-coverage limits. F1 adjusted interval includes zero.", fontsize=9)
        fig.text(.065, .015, "No intervals for other QfO endpoints; unseeded FAS sample. Failed R-on inference timing excluded. Not a selected-default comparison.", fontsize=9)
        fig.subplots_adjust(left=.10, right=.985, top=.81, bottom=.17, hspace=.65, wspace=.42)
        for suffix in ("png", "pdf", "svg"):
            fig.savefig(output / f"native_qfo_p0c0.{suffix}", dpi=200, metadata={"Creator": "OrthoHMM benchmark workflow"})
        plt.close(fig)


def run(snapshot_path, snapshot_sha, binding_path, binding_sha, output, validation_python):
    require(not output.exists() and not output.is_symlink(), "Output already exists")
    evidence = []
    snapshot, snapshot_ref = load(snapshot_path, snapshot_sha, evidence)
    binding, binding_ref = load(binding_path, binding_sha, evidence)
    require(binding["snapshot"] == snapshot_ref, "Plot inputs refer to different snapshots")
    validation = replay_binding(validation_python, snapshot_path, snapshot_sha, binding)
    rows, scores, intervals, admitted = figure_data(snapshot, binding)
    evidence.extend([*binding["evidence"], binding["source"], validation["python"], record(__file__)])
    output.mkdir(parents=True, exist_ok=False)
    for name, values in (("scores.tsv", scores), ("swiss_intervals.tsv", intervals)):
        with (output / name).open("x") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(values[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(values)
    render(rows, intervals, admitted, output)
    for ref in evidence:
        check(ref)
    result = dict(schema="native_qfo_p0c0_figure_v1", snapshot=snapshot_ref, swiss_binding=binding_ref,
        plotted_cells=list(CELLS), admitted_cells_in_snapshot=admitted, plotted_score_endpoints=12,
        plotted_swiss_contrast_endpoints=3, evidence=evidence, outputs=[record(p) for p in sorted(output.iterdir())],
        source=record(__file__), validation=validation, observer_command=sys.orig_argv, python_version=sys.version,
        matplotlib_version=matplotlib.__version__, new_bootstrap_draws=0, new_scoring_or_admission=False,
        scientific_timings_admitted=False, publication_ready=False, visual_review_complete=False)
    with (output / "manifest.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("snapshot", "swiss-binding", "output", "validation-python"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("snapshot-sha256", "swiss-binding-sha256"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    result = run(args.snapshot, args.snapshot_sha256, args.swiss_binding, args.swiss_binding_sha256,
                 args.output, args.validation_python)
    print(json.dumps(dict(plotted_cells=result["plotted_cells"], score_endpoints=12, swiss_endpoints=3)))
