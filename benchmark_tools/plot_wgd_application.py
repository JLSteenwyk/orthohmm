"""Plot audited WGD diagnostics and retain every prespecified example."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from benchmark_tools.assemble_wgd_application import unchanged
from benchmark_tools.run_wgd_application import pinned
from benchmark_tools.snapshot_orthohmm_input_order import record

METHODS = ("orthohmm_high_sensitivity", "orthohmm_satellite_v2", "orthofinder_full", "sonicparanoid")
NAMES = ("OrthoHMM high sensitivity", "OrthoHMM phylogeny", "OrthoFinder full", "SonicParanoid")
COLORS = ("#00858a", "#a33b45", "#303030", "#98700c")
ENDPOINTS = ("separation_rate", "supported_separation_rate", "mean_non_scer_coverage")
TITLES = ("Distinct anchor groups", "Homolog-supported separation", "Mean homolog coverage")


def plot(report):
    if report["cohort_pairs"] != 240 or len(report["uncertainty"]["comparisons"]) != 12:
        raise ValueError("Unexpected frozen population or contrasts")
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    fig.subplots_adjust(left=.19, right=.97, top=.82, bottom=.19, wspace=.23, hspace=.58)
    fig.suptitle("Whole-genome-duplicate application", x=.025, ha="left", y=.975, fontsize=18)
    fig.text(.025, .927, "Development-exposed Saccharomyces data | all panels use the same 231 shared-pillar pairs", fontsize=11)
    fig.text(.025, .887, "Homolog support is not copy-specific orthology truth. Coverage alone can reward merged paralogs.", fontsize=10)
    contrast_names = ("OH phylogeny - OH high", "OH phylogeny - OrthoFinder", "OH phylogeny - Sonic", "OH high - OrthoFinder")
    for column, endpoint in enumerate(ENDPOINTS):
        top, bottom = axes[:, column]
        for y, (method, color) in enumerate(zip(METHODS, COLORS)):
            summary = report["methods"][method]["summary"]["endpoints"][endpoint]
            if summary["pairs"] != 231:
                raise ValueError("Method-dependent plotted population")
            value = summary["mean"] * 100
            top.plot(value, y, "o", color=color, markersize=7)
            top.annotate(f"{value:.1f}", (value, y), xytext=(-8, 9), textcoords="offset points", ha="right", fontsize=9)
        top.set(xlim=(0, 105), ylim=(3.6, -.6), xlabel="Observed percentage")
        top.set_title(f"{'ABC'[column]}  {TITLES[column]}", loc="left", fontsize=11, pad=15)
        top.set_yticks(range(4), NAMES if column == 0 else [""] * 4, fontsize=10)
        rows = [r for r in report["uncertainty"]["comparisons"] if r["endpoint"] == endpoint]
        bottom.axvline(0, color="#999999", linewidth=1)
        for y, row in enumerate(rows):
            adjusted, nominal, point = row["bonferroni12_pp"], row["nominal95_pp"], row["difference_pp"]
            if not np.isfinite([*adjusted, *nominal, point]).all():
                raise ValueError("Unavailable interval cannot be plotted as a score")
            bottom.plot(adjusted, [y, y], color="#303030", linewidth=1)
            bottom.plot(nominal, [y, y], color="#303030", linewidth=4, solid_capstyle="butt")
            bottom.plot(point, y, "o", color="#a33b45", markersize=5)
        bottom.set(xlim=(-90, 90), ylim=(3.6, -.6), xlabel="Difference (percentage points)")
        bottom.set_xticks([-80, -40, 0, 40, 80])
        bottom.set_title(f"{'DEF'[column]}  Paired differences", loc="left", fontsize=11, pad=12)
        bottom.set_yticks(range(4), contrast_names if column == 0 else [""] * 4, fontsize=9)
    for ax in axes.flat:
        ax.grid(axis="x", alpha=.16)
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
    fig.text(.025, .11, "Thick: nominal 95% interval. Thin: Bonferroni interval across 12 contrasts. 20,000 paired pillar resamples; seed 20260920.", fontsize=10)
    fig.text(.025, .065, "The OrthoHMM configurations differ in more than reconciliation; their contrast is not a pure phylogeny ablation.", fontsize=10)
    fig.text(.025, .025, "All 240 experimental pairs, exclusions, the diagnostic MCL checkpoint and individual memberships remain in the supplement.", fontsize=10)
    return fig


def cases(report):
    lines = ["# All Prespecified WGD Examples", "",
             "The six examples were selected by hash before native inference. None was replaced after seeing outcomes.",
             "Support counts describe homologs in each anchor group, not validated ancestral-copy assignments.",
             "High/Low/Sparse are source experimental classes, not orthology confidence labels.", ""]
    names = dict(zip(METHODS, NAMES))
    names["orthofinder_mcl_checkpoint"] = "OrthoFinder MCL checkpoint (diagnostic)"
    for example in report["prespecified_examples"]:
        lines += ["## " + " / ".join(example["orf_pair"]) + " (" + example["experimental_class"] + ")", "",
                  "| Method | State | Support per anchor | Homolog coverage | Pillar groups | Foreign-pillar members |",
                  "| --- | --- | --- | --- | ---: | ---: |"]
        for method, name in names.items():
            row = example["methods"][method]
            support = ", ".join(map(str, row["homolog_support_by_anchor"])) if row["reference_eligible"] else "not evaluated"
            coverage = f"{row['coverage_numerator']}/{row['coverage_denominator']}" if row["reference_eligible"] else "not evaluated"
            groups = str(row["pillar_native_group_count"]) if row["reference_eligible"] else "not evaluated"
            foreign = str(len(row["foreign_pillar_members"])) if row["reference_eligible"] else "not evaluated"
            lines.append(f"| {name} | {row['assignment_state']} | {support} | {coverage} | {groups} | {foreign} |")
        row = example["methods"][METHODS[0]]
        if not row["reference_eligible"]:
            lines += ["", "Reference exclusion retained: " + json.dumps(row["reference_reasons"], sort_keys=True) + "."]
        lines.append("")
    lines += ["## Interpretation Limits", "",
              "These tables establish final membership differences, not which search, tree or reconciliation decision caused them.",
              "A pillar distributed over additional groups can reflect fragmentation, but the reference does not resolve correct ancestral copies.",
              "See the complete machine-readable report for group IDs, sizes, unmapped members, missing assignments and all 240 pairs."]
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", required=True, type=Path)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    audit = pinned({"path": str(args.audit), "sha256": args.sha256})
    if audit["status"] != "native_membership_and_interval_audit_passed":
        raise ValueError("Unadmitted results")
    unchanged(audit["auditor"])
    report = pinned(audit["report"])
    if args.output.exists():
        raise FileExistsError(args.output)
    fig = plot(report)
    args.output.mkdir(parents=True)
    outputs = []
    for extension in ("png", "pdf", "svg"):
        path = args.output / ("wgd_application." + extension)
        fig.savefig(path, dpi=180)
        outputs.append(record(path))
    plt.close(fig)
    path = args.output / "PRESPECIFIED_EXAMPLES.md"
    path.write_text(cases(report))
    outputs.append(record(path))
    manifest = {"audit": record(args.audit), "source_results": audit["report"], "plotter": record(__file__),
                "matplotlib": matplotlib.__version__, "outputs": outputs}
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
