"""Plot every prespecified null tail from retained scores, without native scoring."""

import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_frozen_null_scores import (
    BANDS, ENDPOINTS, LENGTHS, PAIRS_PER_SEED, REGIMES, REVISION, SEEDS, THRESHOLDS, tail,
)
from benchmark_tools.run_simulation_methods import read_frozen

LABELS = ("BLOSUM background", "Uniform residues", "Half background + half glutamine")
COLORS = ("#b34e26", "#007d83")
MAX_BYTES = 4 * 1024 * 1024


def validate(report):
    if (report.get("status") != "prespecified_frozen_null_panel_completed"
            or report.get("revision") != REVISION
            or report.get("independent_pairs") != 90000
            or report.get("native_score_evaluations") != 180000
            or report.get("tail_endpoints") != ENDPOINTS
            or any(report.get(k) is not False for k in (
                "prefilter_executed", "native_pipeline_rerun", "calibration_established",
                "benchmark_scores_or_defaults_changed", "publication_ready"))):
        raise ValueError("Require the completed frozen diagnostic and its unchanged scope")
    rows = {(r["regime"], r["length"], r["seed"]): r for r in report["rows"]}
    expected = {(r, length, seed) for r in REGIMES for length in LENGTHS for seed in SEEDS}
    if len(report["rows"]) != 90 or set(rows) != expected:
        raise ValueError("Changed seed-cell inventory")
    summaries = {(s["regime"], s["length"]): s for s in report["summaries"]}
    if len(report["summaries"]) != 9 or set(summaries) != {(r, l) for r in REGIMES for l in LENGTHS}:
        raise ValueError("Changed composition-length inventory")
    plotted = []
    for regime in REGIMES:
        for length in LENGTHS:
            summary = summaries[regime, length]
            scores = {}
            for band in BANDS:
                pooled = []
                for seed in SEEDS:
                    row = rows[regime, length, seed]
                    values = row["scores"][str(band)]
                    if (set(row["scores"]) != {"0", "64"} or len(values) != PAIRS_PER_SEED
                            or any(type(v) is not int or v < 0 for v in values)):
                        raise ValueError("Changed nonnegative integer score inventory")
                    pooled.extend(values)
                scores[band] = np.asarray(pooled)
                saved = summary["bands"][str(band)]
                if [s["threshold"] for s in saved] != list(THRESHOLDS):
                    raise ValueError("Missing or reordered planned tail endpoints")
                for threshold, original in zip(THRESHOLDS, saved):
                    computed = tail(scores[band], length, threshold)
                    if set(original) != set(computed):
                        raise ValueError("Changed endpoint fields")
                    for key, value in computed.items():
                        a, b = np.asarray(original[key]), np.asarray(value)
                        if a.shape != b.shape or not np.allclose(a, b, rtol=1e-12, atol=1e-10):
                            raise ValueError("Retained endpoint disagrees with raw scores: " + key)
                    low, high = computed["bonferroni_clopper_pearson"]
                    plotted.append(dict(regime=regime, length=length, band=band,
                        threshold=threshold, trials=computed["trials"], hits=computed["hits"],
                        fraction=computed["fraction"], adjusted_low=low, adjusted_high=high,
                        minimum_integer_score=computed["minimum_integer_score"],
                        poisson_model_tail_reference=computed["poisson_model_tail_reference"],
                        zero_hit_upper_limit=computed["hits"] == 0))
            if ((scores[64] > scores[0]).any()
                    or summary["band_changed_scores"] != np.count_nonzero(scores[64] != scores[0])):
                raise ValueError("Changed paired-band evidence")
            for endpoint in summary["bands"]["0"]:
                gate = endpoint["minimum_integer_score"]
                lost = np.count_nonzero((scores[0] >= gate) & (scores[64] < gate))
                if summary["band_lost_gate_hits"][str(endpoint["threshold"])] != lost:
                    raise ValueError("Changed paired lost-hit count")
    return plotted


def load(scores_path, scores_sha, receipt_path, receipt_sha):
    receipt = read_frozen(receipt_path, receipt_sha)
    raw = scores_path.read_bytes()
    if len(raw) > MAX_BYTES or hashlib.sha256(raw).hexdigest() != scores_sha:
        raise ValueError("Changed compressed observations")
    pin = next(p for p in receipt["files"] if p["path"].endswith("/frozen_null_score_observations_20261002.json.gz"))
    if pin["sha256"] != scores_sha or pin["bytes"] != len(raw):
        raise ValueError("Observations do not belong to the pinned receipt")
    with gzip.GzipFile(fileobj=io.BytesIO(raw)) as stream:
        data = stream.read(MAX_BYTES + 1)
    audit = receipt["audit"]
    if (len(data) > MAX_BYTES or len(data) != audit["result_bytes"]
            or hashlib.sha256(data).hexdigest() != audit["result_sha256"]
            or audit["all_endpoint_counts_and_intervals_recomputed"] is not True
            or audit["reference_python_scores_rechecked"] != 180
            or audit["all_native_scores_recomputed"] is not False
            or audit["calibration_established"] is not False):
        raise ValueError("Changed decoded result or audit scope")
    report = json.loads(data)
    return report, validate(report)


def plot(rows):
    figure, axes = plt.subplots(3, 3, figsize=(12, 10), sharex=True, sharey=True)
    figure.subplots_adjust(left=.1, right=.98, top=.82, bottom=.16, hspace=.28, wspace=.14)
    for i, length in enumerate(LENGTHS):
        for j, regime in enumerate(REGIMES):
            ax = axes[i, j]
            reference = sorted((r for r in rows if r["regime"] == regime and r["length"] == length
                                and r["band"] == 0), key=lambda r: r["threshold"])
            ax.plot([r["threshold"] for r in reference],
                    [r["poisson_model_tail_reference"] for r in reference],
                    "--", color="#555555", linewidth=1.1)
            # Offset band glyphs horizontally for visibility; exact cutoffs are in the TSV.
            for band, color, marker, shift in zip(BANDS, COLORS, ("o", "s"), (.91, 1.10)):
                selected = [r for r in rows if r["regime"] == regime and r["length"] == length and r["band"] == band]
                for row in selected:
                    x = row["threshold"] * shift
                    if row["zero_hit_upper_limit"]:
                        ax.plot(x, row["adjusted_high"], marker="v", markersize=5, color=color)
                    else:
                        ax.vlines(x, row["adjusted_low"], row["adjusted_high"], color=color, linewidth=1.2)
                        ax.plot(x, row["fraction"], marker=marker, markersize=4, color=color,
                                markeredgecolor="white", markeredgewidth=.4)
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlim(7e-5, 1.5)
            ax.set_ylim(1e-8, 1.8)
            ax.set_xticks((1e-4, 1e-2, 1))
            ax.set_yticks((1e-6, 1e-4, 1e-2, 1))
            ax.minorticks_off()
            ax.grid(alpha=.18)
            ax.tick_params(labelsize=9)
            ax.text(.95, .055, f"{length} residues", transform=ax.transAxes, ha="right", fontsize=10)
            if i == 0:
                ax.set_title(LABELS[j], fontsize=11, pad=10)
            if i == 2:
                ax.set_xlabel("Approximate E cutoff", fontsize=10)
            if j == 0:
                ax.set_ylabel("Forced-pair tail fraction", fontsize=10)
            for spine in ("top", "right"):
                ax.spines[spine].set_visible(False)
    figure.suptitle("Synthetic null-score tails", x=.06, ha="left", y=.97, fontsize=18)
    figure.text(.06, .925, "Frozen recurrence and significance formula; no parameter fitting", fontsize=11)
    legend = [Line2D([], [], color=color, marker=marker, linestyle="none", label=label)
              for color, marker, label in zip(COLORS, ("o", "s"), ("Full matrix", "Width 64"))]
    legend += [Line2D([], [], color="#555555", linestyle="--", label="Poisson-model reference"),
               Line2D([], [], color="#555555", marker="v", linestyle="none", label="Zero-hit upper bound")]
    figure.legend(handles=legend, loc="upper left", bbox_to_anchor=(.05, .9), frameon=False, ncol=4, fontsize=10)
    figure.text(.06, .087, "10,000 independent pairs per panel; whiskers: exact intervals adjusted over 90 planned endpoints.", fontsize=10)
    figure.text(.06, .056, "Forced scoring only, with one target per query. No prefilter, grouping or orthology inference.", fontsize=10)
    figure.text(.06, .025, "Zero-hit triangles mark interval upper bounds, not positive observed fractions; band glyphs are offset slightly.", fontsize=9)
    return figure


def export(scores_path, scores_sha, receipt_path, receipt_sha, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    inputs = [record(p.resolve()) for p in (scores_path, receipt_path, Path(__file__),
              Path(__file__).with_name("probe_frozen_null_scores.py"))]
    report, rows = load(scores_path, scores_sha, receipt_path, receipt_sha)
    output.mkdir(parents=True, exist_ok=False)
    with (output / "plotted_values.tsv").open("x") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    figure = plot(rows)
    try:
        figure.savefig(output / "frozen_null_scores.png", dpi=180)
        figure.savefig(output / "frozen_null_scores.pdf", metadata={"CreationDate": None, "ModDate": None})
    finally:
        plt.close(figure)
    for pin in inputs:
        check(pin)
    result = dict(status="frozen_null_score_figure_exported", checked_inputs=inputs,
        outputs=[record(p.resolve()) for p in sorted(output.iterdir())],
        revision=report["revision"], endpoints=len(rows), independent_pairs=90000,
        native_evaluations=180000, source_score_counts_recomputed=True,
        native_scoring_rerun=False, calibration_established=False,
        benchmark_scores_or_defaults_changed=False, publication_ready=False,
        figure_visually_reviewed=False,
        limitations=["Exact binomial intervals describe this synthetic generator, not orthology uncertainty.",
            "Zero-count glyphs are upper confidence limits, not positive observed fractions.",
            "Band glyphs are horizontally offset; exact cutoffs and all endpoints remain in the TSV.",
            "Export is not visual inspection, real-data false-positive measurement or calibration."])
    with (output / "manifest.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("scores", "receipt", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("scores-sha", "receipt-sha"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    export(args.scores, args.scores_sha, args.receipt, args.receipt_sha, args.output)
