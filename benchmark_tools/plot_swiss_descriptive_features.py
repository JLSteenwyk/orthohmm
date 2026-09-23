"""Plot pinned descriptive identity/fragment tables without adding subgroup inference."""

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

METHODS = {
    "orthohmm_high_sensitivity": "OrthoHMM high sensitivity",
    "orthohmm_phylogeny_satellite_v2": "OrthoHMM phylogenetic",
    "orthofinder_3_1_5_full": "OrthoFinder 3.1.5 full",
    "orthofinder_3_1_5_sequence_only": "OrthoFinder sequence-only*",
    "sonicparanoid_2_0_9": "SonicParanoid",
    "proteinortho_6_3_6": "ProteinOrtho",
    "fastoma_0_3_5": "FastOMA (supplied tree)",
    "orthomcl_1_4": "OrthoMCL 1.4 (unavailable)",
}
REFERENCE = "orthofinder_3_1_5_full"
METRICS = ("F1", "PPV", "TPR")
IDENTITY_BINS = dict(all=18, higher_identity=9, lower_identity=9, missing_identity=0)
FRAGMENT_BINS = dict(all=18, historical_annotation_positive=5, historical_all_matched_unflagged=13,
    historical_missing_without_positive=0, baseline_only_annotation_positive=5,
    baseline_only_all_matched_unflagged=11, baseline_only_missing_without_positive=2)
SOURCES = {
    "identity": ("swiss_identity_strata_20260923/scores.tsv", "be8ddb6540b62582ebe53e05c33ffdc5893ed10651bfa57eb374d9791b3af664", IDENTITY_BINS),
    "fragment": ("swiss_fragment_strata_20260923/scores.tsv", "4ba79b1e6c9d0f9fabdd0136eddf43fa6ab2c2f9edd2660defe8781ee2aa432b", FRAGMENT_BINS),
}
PANELS = {
    "identity": [("lower_identity", "Lower identity (9 families)"), ("higher_identity", "Higher identity (9 families)")],
    "fragment": [("historical_annotation_positive", "Annotation-positive (5)"),
                 ("historical_all_matched_unflagged", "Historical unflagged (13)"),
                 ("baseline_only_all_matched_unflagged", "Baseline-only unflagged (11)"),
                 ("baseline_only_missing_without_positive", "Baseline-only missing (2)")],
}


def validate(rows, bins):
    expected = {(m, b) for m in METHODS for b in bins}
    indexed = {(r["method"], r["stratum"]): r for r in rows}
    if len(indexed) != len(rows) or set(indexed) != expected:
        raise ValueError("Incomplete or duplicated method/stratum table")
    for (method, stratum), row in indexed.items():
        unavailable = method == "orthomcl_1_4" or bins[stratum] == 0
        status = "method_not_admitted" if method == "orthomcl_1_4" else "empty_bin" if bins[stratum] == 0 else "descriptive"
        if row["status"] != status or int(row["families"]) != bins[stratum]:
            raise ValueError("Changed availability or family count")
        for metric in METRICS:
            for key in (metric, "delta_" + metric):
                value = row[key]
                if unavailable:
                    if value != "NA":
                        raise ValueError("Unavailable value is not missing")
                elif not math.isfinite(float(value)) or not (-1 <= float(value) <= 1):
                    raise ValueError("Invalid numeric endpoint")
            if not unavailable:
                value, reference = float(row[metric]), float(indexed[REFERENCE, stratum][metric])
                if value < 0 or not math.isclose(float(row["delta_" + metric]), value-reference, abs_tol=1e-12, rel_tol=0):
                    raise ValueError("Invalid score or reference difference")
        if not unavailable:
            precision, recall = float(row["PPV"]), float(row["TPR"])
            f1 = 2*precision*recall/(precision+recall) if precision+recall else 0
            if not math.isclose(f1, float(row["F1"]), abs_tol=1e-12, rel_tol=0):
                raise ValueError("F1 is not harmonic macro precision/recall")
    return indexed


def figure(indexed, kind):
    strata = PANELS[kind]
    colors, markers = ("#007f87", "#b1394b", "#555555", "#ad7700"), ("o", "s", "^", "D")
    fig, axes = plt.subplots(1, 3, figsize=(15, 9))
    fig.subplots_adjust(left=.23, right=.98, top=.78, bottom=.22, wspace=.13)
    title = "Sequence identity" if kind == "identity" else "Historical fragment annotations"
    fig.suptitle("QfO SwissTrees: " + title, x=.025, y=.96, ha="left", fontsize=18)
    fig.text(.025, .91, "Descriptive differences versus full OrthoFinder 3.1.5; no subgroup confidence intervals or significance tests.", fontsize=11)
    handles = [Line2D([], [], color=colors[i], marker=markers[i], linestyle="None", label=label)
               for i, (_, label) in enumerate(strata)]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(.02, .88), ncol=2, frameon=False, fontsize=10)
    endpoints = []
    for column, (ax, metric) in enumerate(zip(axes, METRICS)):
        values = [100*float(indexed[m, s]["delta_" + metric]) for m in METHODS for s, _ in strata
                  if indexed[m, s]["delta_" + metric] != "NA"]
        lower, upper = math.floor((min(values)-3)/10)*10, math.ceil((max(values)+3)/10)*10
        ax.axvline(0, color="#888888", linestyle="--", linewidth=1)
        for y, method in enumerate(METHODS):
            for i, (stratum, _) in enumerate(strata):
                row = indexed[method, stratum]
                value = None if row["delta_" + metric] == "NA" else float(row["delta_" + metric])
                endpoints.append(dict(method=method, stratum=stratum, metric=metric, difference=value, status=row["status"]))
                if value is not None:
                    offset = (i-(len(strata)-1)/2)*.15
                    ax.plot(value*100, y+offset, marker=markers[i], color=colors[i], markersize=5, linestyle="None")
            if method == "orthomcl_1_4":
                ax.text(.5, y, "Not yet admitted", transform=ax.get_yaxis_transform(), ha="center", va="center", fontsize=9, color="#666666")
        ax.set(xlim=(lower, upper), ylim=(7.6, -.6), xlabel="Difference (percentage points)")
        ax.set_yticks(range(8), list(METHODS.values()) if column == 0 else [""]*8, fontsize=10)
        ax.set_title(f"{'ABC'[column]}  " + dict(F1="F1", PPV="Precision", TPR="Recall")[metric], loc="left", fontsize=12)
        ax.grid(axis="x", alpha=.15)
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
    notes = ["Development-exposed associations, not causal effects. F1 is harmonic macro precision/recall; panels use separate x-axis ranges.",
             "*Sequence-only OrthoFinder uses group-clique pairs. FastOMA uses a supplied tree. Configuration contrasts are not pure ablations."]
    if kind == "identity":
        notes.append("Family mean identity median: 0.429122; 9 families per bin. Alignment-dependent identity is not calibrated evolutionary distance.")
        notes.append("Missing-identity bin: 0 families. All 18-family scores and empty-bin rows remain in the source table.")
    else:
        notes.append("Unflagged is not proven complete. Baseline-only treats 14 later-version records as missing; two families leave the unflagged bin.")
        notes.append("The positive bin is unchanged in the baseline-only view. Historical missing bin: 0 families. No points represent fabricated zero scores.")
    for y, note in zip((.155, .115, .075, .035), notes):
        fig.text(.025, y, note, fontsize=9)
    return fig, endpoints


def render(root, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    tables, inputs = {}, [record(__file__)]
    for kind, (relative, digest, bins) in SOURCES.items():
        path = root / "benchmark_tools/results" / relative
        item = record(path)
        if item["sha256"] != digest:
            raise ValueError("Changed descriptive source table")
        with path.open() as stream:
            tables[kind] = validate(list(csv.DictReader(stream, delimiter="\t")), bins)
        inputs.append(item)
    output.mkdir(parents=True)
    rows = []
    for kind, table in tables.items():
        fig, endpoints = figure(table, kind)
        try:
            for extension in ("png", "pdf", "svg"):
                fig.savefig(output / f"swiss_{kind}_descriptive.{extension}", dpi=180)
        finally:
            plt.close(fig)
        rows.extend(dict(panel=kind, **r) for r in endpoints)
    with (output / "endpoints.tsv").open("x") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({k: "NA" if v is None else v for k, v in row.items()} for row in rows)
    for item in inputs:
        check(item)
    manifest = dict(inputs=inputs, matplotlib_version=matplotlib.__version__, endpoints=len(rows),
        outputs=[record(p) for p in sorted(output.iterdir())], new_inferential_claims=False,
        publication_ready=False, display_conversion="TSV raw differences; figures multiply by 100")
    with (output / "manifest.json").open("x") as stream:
        json.dump(manifest, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    render(args.root.resolve(), args.output.absolute())
