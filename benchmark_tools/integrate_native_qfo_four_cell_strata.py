"""Integrate every checked profile-bin contrast without altering frozen drafts."""

import argparse
import json
import math
from pathlib import Path

from benchmark_tools.export_native_qfo_three_cell_strata import record, require, write_tsv


PINS = {
    "result": ("native_qfo_four_cell_strata_20261009_v2/report.json",
               "ef3b9de78197a946399543b5aec0b35f0cbe9c89cce643da1a7024cae48fda73"),
    "reader": ("native_qfo_four_cell_strata_readback_20261009_v2.json",
               "1260c1c87ec4d3cbaa606a52225e94fc1fa084a689ba8af9e721f2f14966d1c2"),
    "parent": ("native_qfo_comparator_manuscript_20261009_v1.md",
               "cb6741dcc53f5b67bb20deccdefe5f165c20f001385b227f697e0f19571d2673"),
}
SUITES = ("sequence", "domain", "duplication", "model_distance")
METRICS = ("F1", "PPV", "TPR")
FIELDS = ("suite", "stratum", "families", "status", *METRICS)
METHOD_ANCHOR = "### Fixed-Stratum SwissTrees Error Analyses\n"
RESULT_ANCHOR = "### Fixed-Stratum Errors And Estimated Sequence Divergence\n"


def profile_rows(result, reader):
    require(result["schema"] == "native_qfo_four_cell_strata_v2"
            and reader["schema"] == "native_qfo_four_cell_strata_rational_readback_v2"
            and reader["report_schema_checked"] == result["schema"], "Wrong integrated evidence schema")
    checks = ("family_rows_checked", "score_rows_checked", "differences_checked",
              "inherited_score_rows_reproduced", "inherited_difference_rows_reproduced",
              "human_rows_checked", "compatibility_rows_checked")
    require(tuple(reader[k] for k in checks) == (72, 92, 69, 69, 46, 23, 3)
            and [len(result[k]) for k in ("family_rows", "rows", "differences")] == [72, 92, 69]
            and result["prior_score_rows_reproduced"] == 69 and result["prior_difference_rows_reproduced"] == 46,
            "Incomplete integrated result/readback")
    for doc in (result, reader):
        require(type(doc["new_bootstrap_draws"]) is int and doc["new_bootstrap_draws"] == 0
                and all(doc[k] is False for k in ("publication_ready", "independent_confirmation", "new_uncertainty",
                    "new_accuracy_or_resource_admission", "scientific_timings_admitted")), "Inflated integrated scope")
    selected = [r for r in result["differences"] if r["contrast"] == "P_at_C0_R1"]
    expected = [(s, n) for s in SUITES for n in sorted(result["bins"][s])]
    require(len(expected) == 23 and [(r["suite"], r["stratum"]) for r in selected] == expected,
            "Incomplete or reordered profile bins")
    rows = []
    for row in selected:
        group = result["bins"][row["suite"]][row["stratum"]]
        require(row["candidate"] == "p1_c0_r1" and row["reference"] == "p0_c0_r1"
                and row["family_members"] == group and row["families"] == len(group)
                and row["status"] == ("descriptive" if group else "empty_bin"), "Changed profile contrast")
        require(all(row[m] is None if not group else type(row[m]) in (int, float)
                    and math.isfinite(row[m]) and -1 <= row[m] <= 1 for m in METRICS), "Invalid profile value")
        rows.append({k:row[k] for k in FIELDS})
    return rows


def count_changes(result):
    counts = {(r["cell"], r["family"]):r["counts_without_prior"] for r in result["family_rows"]}
    changes = []
    for family in sorted(result["memberships"]):
        delta = {k:counts["p1_c0_r1", family][k] - counts["p0_c0_r1", family][k] for k in ("TP", "FP", "FN", "TN")}
        if any(delta.values()):
            changes.append(dict(family=family, **delta))
    return changes


def sections(result, reader):
    rows = profile_rows(result, reader)
    changes = count_changes(result)
    index = {(r["suite"], r["stratum"]):r for r in rows}
    require(changes == [dict(family="CASP", TP=0, FP=-3, FN=0, TN=3),
                        dict(family="GH14", TP=-1, FP=0, FN=1, TN=0)], "Changed localized family pattern")
    require(all(index[s, n]["families"] == 9 for s, n in
                (("sequence", "lower_entropy"), ("sequence", "higher_entropy"),
                 ("model_distance", "higher_than_median"), ("model_distance", "lower_or_equal_median")))
            and sum(not r["families"] for r in rows) == 5, "Changed descriptive bin scope")
    for suite, positive, negative in (("sequence", "lower_entropy", "higher_entropy"),
                                      ("model_distance", "higher_than_median", "lower_or_equal_median")):
        require("CASP" in result["bins"][suite][positive] and "GH14" not in result["bins"][suite][positive]
                and "GH14" in result["bins"][suite][negative] and "CASP" not in result["bins"][suite][negative],
                "Changed CASP/GH14 bin localization")
    for suite, changed, unchanged in (("sequence", "not_short_relative", "short_relative"),
            ("domain", "median_pfam_types_below_two", "median_pfam_types_at_least_two"),
            ("domain", "repeated_type_fraction_below_quarter", "repeated_type_fraction_at_least_quarter"),
            ("duplication", "upper_duplication_fraction", "lower_duplication_fraction")):
        require({"CASP", "GH14"} <= set(result["bins"][suite][changed])
                and not {"CASP", "GH14"} & set(result["bins"][suite][unchanged])
                and all(index[suite, unchanged][m] == 0 for m in METRICS), "Changed complementary bin pattern")
    def points(suite, name):
        return f"{100*index[suite, name]['F1']:+.6f}"
    methods = """### Four-Cell Profile Fixed-Stratum Extension

The retrospective extension projects admitted P1/C0/R1 counts into the same
23 fixed bins used for the three profile-off cells: 11 sequence, five domain,
four reference-duplication and three model-distance bins. It adds only the
conditional P1/C0/R1-minus-P0/C0/R1 contrast. Original memberships, cutoffs,
integer counts, halving/unit-prior conversion and macro-family harmonic F1
are unchanged. Empty bins remain NULL/NA. All 69 inherited scores and 46
inherited differences reproduce within absolute 1e-12; the independent
rational readback checks 72 family, 92 score and 69 contrast records.

The original export failed on two retained vocabularies for native pairs.
A separately tested compatibility adapter normalizes only a copied metadata
view and records all three original labels; retained evidence is not edited.
No inference, scoring, feature extraction, bootstrap, subgroup interval,
cutoff tuning, default selection or new admission occurs. Initial HMM search
stays on. The older three-cell descriptions and figure below retain their
original scope; this extension supplies the previously absent profile contrast.
[Frozen protocol](NATIVE_QFO_FOUR_CELL_STRATA_PROTOCOL_20261009.md),
[compatibility amendment](NATIVE_QFO_FOUR_CELL_COMPATIBILITY_AMENDMENT_20261009.md),
[executed independent reader](native_qfo_four_cell_strata_readback_20261009_v2.json).

"""
    lines = ["### Profile Refinement Across All Fixed Bins", "",
        "All 23 conditional profile-bin differences are reported below in percentage",
        "points. Overlapping and repeated all-family bins are audit rows, not independent",
        "findings. NA denotes an empty bin, not zero. No subgroup intervals were computed.", "",
        "| Suite | Bin | Families | F1 Difference | PPV Difference | TPR Difference |",
        "| --- | --- | ---: | ---: | ---: | ---: |"]
    for row in rows:
        values = ["NA" if row[m] is None else f"{100*row[m]:+.3f}" for m in METRICS]
        lines.append("| " + " | ".join([row["suite"], row["stratum"], str(row["families"]), *values]) + " |")
    lines.extend(["", f"Overall F1 changes by {points('sequence', 'all')} percentage points.",
        f"Lower-entropy families have {points('sequence', 'lower_entropy')} points, while higher-entropy",
        f"families have {points('sequence', 'higher_entropy')}. Above-median model-distance families",
        f"have {points('model_distance', 'higher_than_median')}, versus {points('model_distance', 'lower_or_equal_median')}",
        "at or below the fixed median. Each entropy/distance bin has nine families.", "",
        "These patterns reflect the same retained changed families rather than many",
        "independent successes. Profile-minus-reference integer count changes are:", "",
        "| Family | TP | FP | FN | TN |", "| --- | ---: | ---: | ---: | ---: |"])
    for change in changes:
        lines.append("| " + " | ".join([change["family"], *[f"{change[k]:+d}" for k in ("TP", "FP", "FN", "TN")]]) + " |")
    lines.extend(["", f"The other {len(result['memberships'])-len(changes)} families have unchanged counts. CASP contributes",
        "the precision gain in the lower-entropy/higher-distance bins; GH14 contributes",
        "the recall loss in the complementary bins. Both occur in the not-short-relative,",
        "lower-Pfam-type, lower-repeat-fraction and upper-reference-duplication bins.",
        "Their complementary bins are unchanged. Five empty bins support no effect.",
        "The [existing complete pair localization](NATIVE_PROFILE_LOCALIZATION_RESULT_20261007.md)",
        "places all four removed calls at pre-reconciliation candidate separation.",
        "It does not identify a causal HMM edge or establish correct gene/species trees.", "",
        "These are development-exposed descriptive associations, not a general entropy",
        "or divergence benefit, an independent validation, or an initial-HMM effect.",
        "Existing aggregate intervals are not subgroup intervals. Fragment flags, Pfam",
        "descriptors, reference duplication fractions and model distances retain their",
        "proxy limitations. Missing factorial cells and failed/ineligible timing remain.",
        "[Complete four-cell table](native_qfo_four_cell_strata_20261009_v2/TABLE.md),",
        "[full-precision result](native_qfo_four_cell_strata_20261009_v2/report.json),",
        "[actual commands/readback](native_qfo_four_cell_strata_execution_20261009_v2.json).", "", ""])
    return methods, "\n".join(lines), rows, changes


def manuscript(parent, result, reader):
    require(parent.count(METHOD_ANCHOR) == parent.count(RESULT_ANCHOR) == 1, "Ambiguous manuscript anchors")
    methods, results, rows, changes = sections(result, reader)
    require(methods not in parent and results not in parent, "Already integrated")
    revised = parent.replace(METHOD_ANCHOR, methods + METHOD_ANCHOR).replace(RESULT_ANCHOR, results + RESULT_ANCHOR)
    require(revised.replace(methods, "", 1).replace(results, "", 1) == parent, "Changed frozen parent body")
    return revised, methods, results, rows, changes


def run(root, output, revised_path):
    directory = Path(root).resolve() / "benchmark_tools/results"
    output, revised_path = Path(output).absolute(), Path(revised_path).absolute()
    for path in (output, revised_path):
        if path.exists() or path.is_symlink():
            raise FileExistsError(path)
    require(output.parent.resolve() == revised_path.parent.resolve() == directory, "Preserve manuscript-relative links")
    docs, refs = {}, {}
    for key, (name, digest) in PINS.items():
        path = directory / name
        refs[key] = record(path)
        require(refs[key]["sha256"] == digest, "Changed integrated input: " + key)
        docs[key] = path.read_text() if key == "parent" else json.loads(path.read_text())
    require(docs["reader"]["report"] == refs["result"], "Mixed executed reader/report binding")
    revised, methods, results, rows, changes = manuscript(docs["parent"], docs["result"], docs["reader"])
    source = record(__file__)
    for ref in [*refs.values(), source]:
        require(record(ref["path"]) == ref, "Changed input before writing")
    output.mkdir(parents=True)
    (output / "methods_section.md").write_text(methods)
    (output / "results_section.md").write_text(results)
    write_tsv(output / "profile_bins.tsv", rows, FIELDS)
    revised_path.write_text(revised)
    for ref in [*refs.values(), source]:
        require(record(ref["path"]) == ref, "Changed input while writing")
    manifest = dict(schema="native_qfo_four_cell_strata_integration_v1", inputs=refs, source=source,
                    outputs=[record(p) for p in sorted(output.iterdir())], manuscript=record(revised_path),
                    profile_bin_rows=23, profile_value_cells=69, changed_family_counts=changes,
                    parent_body_unchanged_except_insertions=True, table_scale=100, tsv_scale=1,
                    new_bootstrap_draws=0, new_inference_or_scoring=False, publication_ready=False,
                    manuscript_rendered=False, visual_reviewed=False)
    with (output / "manifest.json").open("x") as stream:
        json.dump(manifest, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "output", "manuscript"):
        parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.root, args.output, args.manuscript), allow_nan=False))
