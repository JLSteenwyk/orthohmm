"""Integrate admitted four-cell native QfO evidence without rewriting prior snapshots."""

import argparse
import hashlib
import json
import math
from pathlib import Path

PINS = {
    "parent": ("PUBLICATION_MAIN_TEXT_20261007_v3.md", "41eab259cadea1d0339bf262cb1ad50a488aec4f29031f6aecce3e8d038c4951"),
    "snapshot": ("native_qfo_scientific_scores_20261007_v3/report.json", "7aab1cbd31fb650df42e6e80a14ee0167b35a810f08202768c89c514755211a2"),
    "binding": ("native_qfo_four_cell_swiss_uncertainty_20261007_v1.json", "e8f554a178841da24a1892a7e42f8680f9fe8e4827dec6d6bfe4643cdb890bac"),
    "rational": ("native_qfo_profile_swiss_readback_20261007_v1.json", "3f62b487014cb9ba99ea36b5adc15535f62c7e39c9a99f3cb4eb8c537523005d"),
    "figure": ("native_qfo_four_cell_figure_20261007_v1/manifest.json", "4732f525575b03d548bdf18a4047625b565ac1f28df005513d9457ac86f36d75"),
    "figure_reader": ("native_qfo_four_cell_figure_readback_20261007_v1.json", "71aa0807750223a22f7e6080db4de01a3a9f7ccae69c9d2ff5170418c514b614"),
    "localization": ("native_qfo_profile_localization_20261007_v1.json", "0d09567de8828b46312aa280f0bdabe6bbfbfaa01f591cbbca9190241bdb5ae0"),
    "tree_reader": ("native_qfo_profile_newick_readback_20261007_v1.json", "289ac094115f493b8607073b1395ca2d97493ba4a00d5f1a083ceb3aa030507f"),
    "visual": ("NATIVE_FOUR_CELL_FIGURE_VISUAL_REVIEW_20261007.md", "e61b3dc1b87f016bcb37a1396a1f8f06962fa9b3abe284dec94ce7aade470b3a"),
}
CELLS = ("p0_c0_r0", "p0_c0_r1", "p0_c1_r0", "p1_c0_r1")
ENDPOINTS = ("VGNC", "SwissTrees", "TreeFam-A", "GO", "EC", "FAS")
CONTRASTS = ("C_at_P0_R0", "R_at_P0_C0", "P_at_C0_R1")
START = "### Fresh Native QfO Trade-Off And Pair Composition\n\n"
END = "Candidate expansion recovers 162 asserted VGNC TPs and adds 2,133 scored\n"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def record(path):
    path = Path(path).resolve(strict=True)
    data = path.read_bytes()
    return dict(path=str(path), bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def inputs(directory):
    docs, refs = {}, {}
    for key, (name, digest) in PINS.items():
        path = directory / name
        ref = record(path)
        require(ref["sha256"] == digest, "Changed manuscript evidence: " + key)
        refs[key] = ref
        docs[key] = path.read_text() if key in ("parent", "visual") else json.loads(path.read_text())
    return docs, refs


def validate(docs, refs):
    expected = dict(snapshot="allocated_native_qfo_scientific_reporting_snapshot_v1",
        binding="allocated_native_qfo_retained_swiss_uncertainty_binding_v1",
        rational="allocated_native_qfo_profile_swiss_rational_readback_v1",
        figure="native_qfo_four_cell_figure_v1", figure_reader="native_qfo_four_cell_figure_readback_v1",
        localization="allocated_native_qfo_profile_pair_localization_v1",
        tree_reader="allocated_native_qfo_profile_newick_readback_v1")
    for key, schema in expected.items():
        require(docs[key]["schema"] == schema and docs[key]["publication_ready"] is False, "Changed input scope: " + key)
    snapshot, binding, figure = (docs[key] for key in ("snapshot", "binding", "figure"))
    rows = {row["cell"]: row for row in snapshot["rows"]}
    require(len(snapshot["rows"]) == len(rows) == 7
            and {c for c, r in rows.items() if r["accuracy_admitted"] is True} == set(CELLS)
            and all(all(v is None for v in r["scores"].values()) for c, r in rows.items() if c not in CELLS),
            "Changed admitted or unavailable cells")
    require(all(set(rows[c]["scores"]) == set(ENDPOINTS) for c in CELLS), "Changed endpoint inventory")
    for cell in CELLS:
        row = rows[cell]
        require(all(type(v) in (int, float) and math.isfinite(v) and 0 <= v <= 1 for v in row["scores"].values())
                and math.isclose(sum(row["scores"].values()) / 6, row["secondary_mean"], rel_tol=0, abs_tol=1e-12),
                "Invalid score or secondary mean")
        require(type(row["input_accessions"]) is int and row["input_accessions"] == 984137
                and type(row["relation_accessions"]) is int and 0 <= row["relation_accessions"] <= row["input_accessions"]
                and type(row["submitted_pairs"]) is int and row["submitted_pairs"] >= 0
                and row["relation_coverage"] == row["relation_accessions"] / row["input_accessions"],
                "Invalid native coverage/count scope")
    require(rows[CELLS[1]]["resources"] is None and rows[CELLS[1]]["timing_eligible"] is False
            and rows[CELLS[1]]["timing_admitted"] is False
            and rows[CELLS[3]]["scientific_timings_admitted"] is False, "Manuscript repairs failed or shared timing")
    require(binding["snapshot"] == figure["snapshot"] == refs["snapshot"]
            and figure["swiss_binding"] == docs["rational"]["binding"] == refs["binding"]
            and figure["profile_readback"] == refs["rational"]
            and docs["figure_reader"]["manifest"] == refs["figure"]
            and docs["tree_reader"]["report"] == refs["localization"], "Mixed manuscript evidence bindings")
    require(binding["replicates_reused"] == 100000 and binding["seed_reused"] == 20260922
            and binding["multiplicity_endpoints"] == 42 and type(binding["new_bootstrap_draws"]) is int
            and binding["new_bootstrap_draws"] == 0
            and len(binding["families"]) == 18, "Changed interval protocol")
    effects = {row["name"]: row for row in binding["contrasts"]}
    require(len(effects) == len(binding["contrasts"]) == 14
            and {c for c, r in effects.items() if r["status"] == "native_records_matched"} == set(CONTRASTS)
            and all(r["metrics"] is None for c, r in effects.items() if c not in CONTRASTS)
            and docs["rational"]["contrast"] == effects["P_at_C0_R1"], "Changed supported contrast inventory")
    require((docs["figure_reader"]["score_endpoints_checked"], docs["figure_reader"]["coverage_endpoints_checked"],
             docs["figure_reader"]["contrast_endpoints_checked"]) == (24, 4, 9)
            and docs["figure_reader"]["exact_table_readback"] is True, "Incomplete figure readback")
    require(docs["localization"]["changed_pairs_traced"] == docs["tree_reader"]["changed_pairs_checked"] == 4
            and docs["localization"]["summary"] == docs["tree_reader"]["summary"] == {"candidate_separation": 4}
            and docs["tree_reader"]["source_families_checked"] == 6
            and docs["tree_reader"]["tree_leaves_checked"] == 70
            and docs["tree_reader"]["distinct_lcas_checked"] == 2
            and docs["localization"]["species_tree_bytes_identical"] is False
            and docs["tree_reader"]["species_tree_bytes_identical"] is False, "Changed localization/readback scope")
    return rows, effects


def section(docs, refs):
    rows, effects = validate(docs, refs)
    lines = [START.rstrip(), "", "Four of seven fresh corrected-QfO cells have admitted accuracy; three remain",
        "unavailable in this dated snapshot. Initial sensitive HMM search stays on",
        "in every cell. P denotes downstream profile refinement, C candidate expansion",
        "and R phylogenetic reconciliation. This is not the selected-default all-tool",
        "comparison above or a total-HMM-off control. The recovered P0/C0/R1 result",
        "retains failed timing, null resources and timing ineligibility. Accuracy",
        "admission does not repair that failure.", "",
        "| Endpoint | Statistic | P0/C0/R0 | P0/C0/R1 | P0/C1/R0 | P1/C0/R1 |",
        "| --- | --- | ---: | ---: | ---: | ---: |"]
    for endpoint in ENDPOINTS:
        statistic = "F1" if endpoint in ENDPOINTS[:3] else "Sample mean" if endpoint == "FAS" else "Similarity"
        lines.append("| " + " | ".join([endpoint, statistic, *[f"{rows[c]['scores'][endpoint]:.8f}" for c in CELLS]]) + " |")
    lines.extend(["", "[Admitted native endpoint snapshot](native_qfo_scientific_scores_20261007_v3/report.json).",
        "The project-defined secondary means are " + ", ".join(f"{rows[c]['secondary_mean']:.8f}" for c in CELLS)
        + ", respectively, not official QfO F1 or a single accuracy ranking.",
        "Prediction coverage uses all 984,137 input accessions, not the number of",
        "assessed reference pairs. Mapping loss is zero in these conversions.", "",
        "| Cell | Submitted pairs | Input-accession coverage | Output semantics |",
        "| --- | ---: | ---: | --- |"])
    for cell in CELLS:
        label = cell.upper().replace("_", "/")
        semantics = "Inferred ortholog pairs" if cell.endswith("r1") else "Group-clique pairs"
        lines.append(f"| {label} | {rows[cell]['submitted_pairs']:,} | {100*rows[cell]['relation_coverage']:.4f}% | {semantics} |")
    lines.extend(["", "Three of 14 planned contrasts have exactly matched native SwissTrees family",
        "records; the other 11 remain unavailable. The original 100,000 paired-family",
        "draws (seed 20260922) and all 42 adjusted endpoints are reused, with no new",
        "draws. Effects and adjusted intervals below are percentage points; wins/ties/",
        "losses count the 18 families, not independent gene pairs.", "",
        "| Contrast | Endpoint | Difference (pp) | Adjusted interval (pp) | Wins/ties/losses |",
        "| --- | --- | ---: | --- | --- |"])
    for name, label in zip(CONTRASTS, ("C at P0 R0", "R at P0 C0", "P at C0 R1")):
        for metric, endpoint in (("F1", "F1"), ("PPV", "Precision"), ("TPR", "Recall")):
            values = effects[name]["metrics"][metric]
            interval = "[" + ", ".join(f"{100*x:.4f}" for x in values["bonferroni_percentile_ci"]) + "]"
            signs = "/".join(str(values[k]) for k in ("family_wins", "family_ties", "family_losses"))
            lines.append(f"| {label} | {endpoint} | {100*values['difference']:+.4f} | {interval} | {signs} |")
    lines.extend(["", "All three adjusted F1 intervals include zero. Candidate expansion raises recall",
        "while lowering precision; reconciliation shifts that trade-off in the opposite",
        "direction. Profile refinement at C0/R1 has a slightly lower F1 point estimate,",
        "with one family improving, 16 tying and one worsening. These are conditional",
        "pipeline contrasts, not identified main effects, interactions, equivalence or",
        "independent confirmation. Count-based arithmetic differs slightly from native",
        "serialized endpoints and remains separately reported. Development exposure,",
        "family exchangeability and approximate percentile coverage limit the intervals;",
        "they do not apply to other QfO endpoints or the secondary mean.",
        "[Guarded interval binding](native_qfo_four_cell_swiss_uncertainty_20261007_v1.json),",
        "[independent arithmetic readback](native_qfo_profile_swiss_readback_20261007_v1.json).", "",
        "**Native QfO Figure.** The [reviewed six-panel figure](native_qfo_four_cell_figure_20261007_v1/native_qfo_four_cell.pdf)",
        "separates F1 from functional similarities and coverage, and shows the three",
        "condition-specific SwissTrees contrasts. It is a four-cell ablation, not a",
        "complete fresh factorial or a timing comparison. All 24 endpoint scores, four",
        "coverage rows and nine intervals were separately read back; both actual PNG",
        "and decoded PDF preview were visually inspected.",
        "[Figure evidence](NATIVE_FOUR_CELL_FIGURE_RESULT_20261007.md).", "",
        "Every changed assessed SwissTrees profile pair was traced, avoiding favorable",
        "example selection. Three CASP false positives disappeared and one GH14 true",
        "positive was lost. All four pairs shared a candidate family with a speciation",
        "pair-event before refinement, then spanned different candidate families",
        "afterward. Thus these absences already follow from candidate separation before",
        "reconciliation, not a new duplication exclusion at a shared LCA. Saved Newick",
        "and independent root-membership checks cover six families, 70 leaves and the",
        "two before-state LCAs. Species-tree byte hashes also differ. These checks",
        "explain software inclusion, not biological correctness, individual HMM-edge",
        "causality or an isolated species-tree effect. Equal family weighting and",
        "different precision/recall changes explain why removing three false positives",
        "while losing one true positive need not improve macro F1.",
        "[Complete localization and limits](NATIVE_PROFILE_LOCALIZATION_RESULT_20261007.md),",
        "[independent saved-tree check](native_qfo_profile_newick_readback_20261007_v1.json).", "",
        "The [earlier three-cell table](native_qfo_scientific_scores_20261006_v2/report.json),",
        "[binding](native_qfo_candidate_swiss_uncertainty_20261006_v1.json) and",
        "[figure](native_qfo_three_cell_figure_20261006_v1/native_qfo_three_cell.pdf)",
        "remain historical and unchanged. Their fixed-bin diagnostics below still",
        "describe P0 cells and do not estimate the new profile-refinement effect.", ""])
    return "\n".join(lines) + "\n"


def manuscript(docs, refs):
    parent = docs["parent"]
    require(parent.count("## Abstract\n") == parent.count(START) == parent.count(END) == 1,
            "Missing or duplicate manuscript anchors")
    body = parent[parent.index("## Abstract\n"):]
    a, b = body.index(START), body.index(END)
    body = body[:a] + section(docs, refs) + body[b:]
    replacements = (
        ("Every R-on GO scored pair (78,607) occurs in the R-off set (145,142), and every\n",
         "For P0/C0, every R-on GO scored pair (78,607) occurs in the R-off set (145,142), and every\n"),
        ("They do not contain this later three-cell/complete-path\nsource snapshot.",
         "They do not contain this later four-cell native QfO source snapshot."),
        ("This source is [mechanically extended](../prepare_error_strata_main_text.py)\n",
         "This source is [mechanically extended](../prepare_four_cell_main_text.py)\n"),
    )
    for old, new in replacements:
        require(body.count(old) == 1, "Changed ancillary manuscript anchor")
        body = body.replace(old, new, 1)
    return ("# OrthoHMM: HMM-Centered Group Inference With Phylogenetic Refinement\n\n"
        "Native QfO evidence revision of 7 October 2026, v4. Not submission-ready.\n"
        "This version preserves the [reviewed v3 source](PUBLICATION_MAIN_TEXT_20261007_v3.md),\n"
        "its 21-page PDF, earlier figures and archives. The native QfO section now\n"
        "integrates four admitted cells, three supported conditional contrasts and\n"
        "the independently checked profile-pair localization. All-tool comparisons,\n"
        "frozen methods, prior P0 strata and citations retain their previous scope.\n"
        "No new scoring, inference, bootstrap or default selection is performed.\n"
        "A new source does not itself establish render review or archive inclusion.\n\n" + body)


def generate(directory, output, receipt):
    require(output.parent.resolve() == directory.resolve(), "Keep manuscript beside relative evidence links")
    require(output != receipt and all(not p.exists() and not p.is_symlink() for p in (output, receipt)),
            "Require distinct fresh manuscript and receipt")
    docs, refs = inputs(directory)
    text = manuscript(docs, refs)
    for ref in refs.values():
        require(record(ref["path"]) == ref, "Evidence changed during generation")
    with output.open("x") as stream:
        stream.write(text)
    result = dict(schema="four_cell_native_main_generation_v1", source=record(__file__),
        inputs=refs, output=record(output), native_cells=list(CELLS), endpoint_rows=6, coverage_rows=4,
        conditional_interval_rows=9, missing_score_cells=3, missing_contrasts=11,
        historical_source_modified=False, scientific_evidence_replayed=False,
        new_scoring_or_admission=False, new_bootstrap_draws=0, independent_confirmation=False,
        publication_ready=False, render_review_established=False)
    with receipt.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("evidence-directory", "output", "receipt"):
        parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()
    result = generate(args.evidence_directory, args.output, args.receipt)
    print(json.dumps({k: result[k] for k in ("native_cells", "endpoint_rows", "coverage_rows", "conditional_interval_rows")}))


if __name__ == "__main__":
    main()
