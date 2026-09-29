"""Plot all 21 retained OrthoBench contrasts without recalculating intervals."""

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

BASELINE = "orthofinder_3_1_5_full"
METHODS = {
    "orthohmm_high_sensitivity": "OrthoHMM sensitive",
    "orthohmm_phylogeny_satellite_v2": "OrthoHMM phylogenetic",
    "orthofinder_3_1_5_sequence_only": "OrthoFinder sequence-only",
    "sonicparanoid_2_0_9": "SonicParanoid",
    "proteinortho_6_3_6": "ProteinOrtho",
    "fastoma_0_3_5": "FastOMA",
    "orthomcl_1_4": "OrthoMCL",
}
METRICS = ("f_score", "precision", "recall")


def endpoints(result):
    if (result["status"] != "complete_ob_exploratory_paired_uncertainty"
            or result["baseline"] != BASELINE or result["replicates"] != 100000
            or result["seed"] != 20260928 or result["alpha"] != .05
            or result["multiplicity"] != "Bonferroni tail adjustment over 21 planned endpoints"
            or result["families"] != [f"RefOG{i:03d}.txt" for i in range(1, 71)]
            or set(result["comparisons"]) != set(METHODS)
            or set(result["point_estimates_percent"]) != {BASELINE, *METHODS}):
        raise ValueError("Require the complete retained analysis design")
    rows = []
    for method in METHODS:
        contrast = result["comparisons"][method]
        if contrast["versus"] != BASELINE or set(contrast["metrics"]) != set(METRICS):
            raise ValueError("Wrong contrast inventory")
        for metric in METRICS:
            values = contrast["metrics"][metric]
            difference = values["difference_percentage_points"]
            nominal, adjusted = values["paired_percentile_ci"], values["bonferroni_percentile_ci"]
            if len(nominal) != 2 or len(adjusted) != 2:
                raise ValueError("Wrong interval shape")
            points = [result["point_estimates_percent"][m][metric] for m in (method, BASELINE)]
            if any(type(v) not in (float, int) or not math.isfinite(v)
                   for v in [difference, *nominal, *adjusted, *points]):
                raise ValueError("Nonfinite or nonnumeric endpoint")
            if not -100 <= adjusted[0] <= nominal[0] <= nominal[1] <= adjusted[1] <= 100:
                raise ValueError("Invalid or nonnested intervals")
            if any(not 0 <= p <= 100 for p in points) or abs(points[0] - points[1] - difference) > 1e-8:
                raise ValueError("Difference does not match retained point estimates")
            rows.append(dict(method=method, reference=BASELINE, metric=metric,
                             difference_pp=difference, nominal_low_pp=nominal[0], nominal_high_pp=nominal[1],
                             adjusted_low_pp=adjusted[0], adjusted_high_pp=adjusted[1]))
    return rows


def run(path, digest, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    source = record(path)
    if source["sha256"] != digest:
        raise ValueError("Result checksum differs")
    result = json.loads(path.read_text())
    rows = endpoints(result)
    refs = [source, record(Path(__file__)), record(Path(__file__).with_name("prepare_ob_candidate_neighborhood.py"))]
    output.mkdir(parents=True, exist_ok=False)
    with (output / "endpoints.tsv").open("x") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 10, "svg.hashsalt": "ob-complete-21"}):
        fig, axes = plt.subplots(1, 3, figsize=(15, 6), sharey=True)
        for axis, metric, title, color in zip(axes, METRICS, ("F1", "Precision", "Recall"),
                                              ("#167d8d", "#bd4148", "#9a710c")):
            axis.axvline(0, color="#777777", linestyle="--", linewidth=.8)
            for y, row in enumerate(r for r in rows if r["metric"] == metric):
                axis.plot([row["adjusted_low_pp"], row["adjusted_high_pp"]], [y, y], color=color, linewidth=1.4)
                axis.plot([row["nominal_low_pp"], row["nominal_high_pp"]], [y, y], color=color, linewidth=4)
                axis.plot(row["difference_pp"], y, "o", color=color, markersize=5)
            axis.set(title=title, xlim=(-80, 60), xlabel="Difference (percentage points)")
            axis.set_xticks([-80, -40, 0, 40])
            axis.spines[["top", "right", "left"]].set_visible(False)
            axis.grid(axis="y", color="#eeeeee", linewidth=.7)
            axis.tick_params(axis="y", length=0)
        axes[0].set_yticks(range(7), METHODS.values())
        axes[0].set_ylim(6.6, -.6)
        fig.suptitle("OrthoBench: seven methods minus full OrthoFinder 3.1.5", x=.025, ha="left", fontsize=16)
        fig.text(.025, .91, "70 RefOGs; 100,000 paired draws. Thick: nominal 95% interval. Thin: 21-endpoint adjusted interval.")
        fig.text(.025, .045, "Development-exposed, conditional on family exchangeability. Intervals are approximate; inclusion of zero does not establish equivalence.", fontsize=9)
        fig.text(.025, .018, "Sequence-only OrthoFinder is a pre-phylogeny diagnostic. Fixed method order; no new inference, resampling or significance tests.", fontsize=9)
        fig.subplots_adjust(left=.225, right=.985, top=.83, bottom=.18, wspace=.18)
        for suffix in ("png", "pdf", "svg"):
            fig.savefig(output / f"ob_complete_uncertainty.{suffix}", dpi=180)
        plt.close(fig)
    for ref in refs:
        check(ref)
    manifest = dict(status="complete_ob_uncertainty_plotted", inputs=refs, endpoints=len(rows),
                    units="percentage_points", outputs=[record(p) for p in sorted(output.iterdir())],
                    matplotlib_version=matplotlib.__version__, new_statistics=False,
                    publication_ready=False, visual_review_complete=False)
    with (output / "manifest.json").open("x") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.results, args.sha256, args.output)
