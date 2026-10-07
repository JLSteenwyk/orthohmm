"""Generate a new provisional manuscript from bound, already-reviewed summaries."""

import argparse
import hashlib
import json
import math
from pathlib import Path

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


ROOT = Path(__file__).resolve().parents[1]
BASE = Path("benchmark_tools/results")
PINS = {
    "original": ("PUBLICATION_MAIN_TEXT_20261006.md", "748cf8ea335ef7bc0fedd10b583498e7dc29d4981e065f62dec482f55e524681"),
    "scores": ("native_qfo_scientific_scores_20261006_v2/report.json", "7916d3e23808edbb92b016b5c40ac56b9417a863d1af01a7b93a8ef5dbb53a63"),
    "swiss": ("native_qfo_candidate_swiss_uncertainty_20261006_v1.json", "76e6f4d33a13a12b7f0a396f6b795b63b654666a038485fb8af2538259f553c2"),
    "swiss_readback": ("native_qfo_candidate_swiss_readback_20261006_v1.json", "07941ba9c55ba7dc14c5ba53df37d28cb3b2f60f599b9d177b10c9afa99c1329"),
    "vgnc": ("native_qfo_candidate_vgnc_20261006_v1/report.json", "0c83ab917417ac4e97570b8a68b4bbb0bdfd72eeb301b930bf24446e1a623aea"),
    "vgnc_readback": ("native_qfo_candidate_vgnc_readback_20261006_v1.json", "f75cc5775eb63f5a4aa9c72d844782f285d39b00481a380dca6450669c109c57"),
    "aliases": ("native_qfo_candidate_alias_group_20261006_v1/report.json", "3b3ab8f1132e70362e1450ac990caf7a9dad536e528bc54d6af39db16183cd53"),
    "aliases_readback": ("native_qfo_candidate_alias_group_readback_20261006_v1.json", "4afc4aa6853a69027b9dc710602309f44923315b1e68dbae4cd1e61866f4b92d"),
    "ledger": ("native_qfo_candidate_alias_group_20261006_v1/changed_pair_groups.tsv", "3f0598f7e00cfc0942146604ab1863ad9f03f2d5017efd5244c7e76a29d43642"),
    "figure": ("native_qfo_three_cell_figure_20261006_v1/manifest.json", "bb55259208cf4c45949a64923bd601e1e150331d4542b4e3464bd34f27d742d3"),
    "figure_pdf": ("native_qfo_three_cell_figure_20261006_v1/native_qfo_three_cell.pdf", "c3e6ca725b1a0e48653179956f87d41c8f957093517d9610324843387f8613f1"),
}
CELLS = ("p0_c0_r0", "p0_c0_r1", "p0_c1_r0")
ENDPOINTS = ("VGNC", "SwissTrees", "TreeFam-A", "GO", "EC", "FAS")
NATIVE_START = "### Fresh Native QfO Trade-Off And Pair Composition\n"
NATIVE_END = "### Synthetic Null Tails Depend On Composition\n"
FUNCTIONAL_START = "Every R-on GO scored pair (78,607)"
AVAILABILITY_START = "The second revision's 34-page review/rc3 archive"


def bounded(text, start, end):
    if text.count(start) != 1 or text.count(end) != 1:
        raise ValueError("Manuscript markers missing or ambiguous")
    left, rest = text.split(start, 1)
    middle, right = rest.split(end, 1)
    return left, middle, right


def close(left, right):
    return math.isclose(left, right, rel_tol=0, abs_tol=1e-12)


def validate(scores, swiss, vgnc, aliases, reader):
    if scores["publication_ready"] is not False or scores["new_scoring_or_admission"] is not False:
        raise ValueError("Only existing provisional scores may be presented")
    rows = scores["rows"]
    admitted = [r for r in rows if r["accuracy_admitted"]]
    if ([r["index"] for r in rows] != list(range(6, 13))
            or tuple(r["cell"] for r in admitted) != CELLS
            or len({r["cell"] for r in rows}) != len(rows)):
        raise ValueError("Unexpected native admission cohort")
    for row in rows:
        if set(row["scores"]) != set(ENDPOINTS):
            raise ValueError("Unexpected endpoint selection")
        if not row["accuracy_admitted"]:
            if any(v is not None for v in row["scores"].values()) or row["secondary_mean"] is not None:
                raise ValueError("Unadmitted native score must remain unavailable")
            continue
        for value in row["scores"].values():
            if type(value) not in (float, int) or not math.isfinite(value) or not 0 <= value <= 1:
                raise ValueError("Invalid admitted score")
        if not close(row["secondary_mean"], sum(row["scores"].values()) / len(ENDPOINTS)):
            raise ValueError("Secondary mean differs from retained endpoints")
        if not close(row["relation_coverage"], row["relation_accessions"] / row["input_accessions"]):
            raise ValueError("Coverage denominator differs")
    recovered = admitted[1]
    if (recovered["resources"] is not None or recovered["timing_eligible"] is not False
            or recovered["timing_admitted"] is not False):
        raise ValueError("Recovered scientific score must not repair failed timing")
    matched = [c for c in swiss["contrasts"] if c["metrics"] is not None]
    if (set(swiss["bound_cells"]) != set(CELLS) or len(swiss["contrasts"]) != 14
            or {c["name"] for c in matched} != {"C_at_P0_R0", "R_at_P0_C0"}
            or swiss["replicates_reused"] != 100000 or swiss["seed_reused"] != 20260922
            or swiss["multiplicity_endpoints"] != 42 or swiss["new_bootstrap_draws"] != 0
            or len(swiss["families"]) != 18 or swiss["publication_ready"] is not False):
        raise ValueError("Unexpected retained SwissTrees contrast scope")
    for contrast in matched:
        if set(contrast["metrics"]) != {"F1", "PPV", "TPR"}:
            raise ValueError("Unexpected SwissTrees metric selection")
        for metric in contrast["metrics"].values():
            low, high = metric["bonferroni_percentile_ci"]
            if not math.isfinite(low) or not math.isfinite(high) or low > high:
                raise ValueError("Invalid retained interval")
            if sum(metric[k] for k in ("family_wins", "family_ties", "family_losses")) != 18:
                raise ValueError("Incomplete family comparison")
        low, high = contrast["metrics"]["F1"]["bonferroni_percentile_ci"]
        if not low <= 0 <= high:
            raise ValueError("Frozen F1 zero-crossing statement no longer valid")
    if (aliases["localized_summary"] != reader["localized_summary"]
            or aliases["bridges"] != reader["bridges"]
            or any(aliases[k] is not True or reader[k] is not True for k in (
                "complete_pair_localization", "whole_candidate_partition_reconstructed",
                "original_scored_identifiers_preserved"))
            or any(aliases[k] is not False or reader[k] is not False for k in (
                "new_scoring_or_admission", "accuracy_rescored", "uncertainty_admitted", "publication_ready"))):
        raise ValueError("Missing independently matched complete path evidence")
    paths = {}
    for row in aliases["localized_summary"]:
        key = (row["iteration"], row["connection_path"], row["candidate_state"])
        if (key in paths or key[0] not in (0, 1)
                or key[1] not in ("direct_cross_endpoint", "transitive_union")
                or key[2] not in ("TP", "FP") or type(row["pairs"]) is not int or row["pairs"] < 0):
            raise ValueError("Unexpected grouping path summary")
        paths[key] = row["pairs"]
    if len(paths) != 8 or sum(paths.values()) != aliases["changed_pairs"]:
        raise ValueError("Incomplete changed-pair localization")
    transitions = {(r["baseline"], r["candidate"]): r["pairs"] for r in vgnc["transition_counts"]}
    if len(transitions) != len(vgnc["transition_counts"]) or set(transitions) != {
            ("FN", "FN"), ("FN", "TP"), ("FP", "FP"), ("TP", "TP"), ("not_scored", "FP")}:
        raise ValueError("Unexpected reference-state transitions")
    for state, transition in (("TP", ("FN", "TP")), ("FP", ("not_scored", "FP"))):
        if sum(n for (_, _, s), n in paths.items() if s == state) != transitions[transition]:
            raise ValueError("Group paths differ from scored transition totals")
    for row, which in ((admitted[0], "baseline"), (admitted[2], "candidate")):
        if not close(row["scores"]["VGNC"], vgnc[which]["metrics"]["f1"]):
            raise ValueError("VGNC decomposition differs from admission")
    if aliases["baseline_groups"] - aliases["candidate_groups"] != aliases["accepted_merges"]:
        raise ValueError("Whole-group reconstruction differs")
    return admitted, matched, paths, transitions


def label(cell):
    return cell.upper().replace("_", "/")


def native_section(admitted, matched, paths, transitions, vgnc, aliases):
    text = [NATIVE_START, "Three of seven fresh corrected-QfO cells have admitted accuracy; four remain\n"
        "unavailable in this dated snapshot. All shown cells have downstream P0 and\n"
        "initial sensitive HMM search on. They are not the selected-default all-tool\n"
        "comparison above or an initial-HMM-off control. The recovered R1 scientific\n"
        "result retains original failed 22437 timing, null resources and timing\n"
        "ineligibility. Scientific admission does not repair that failed measurement.\n",
        "| Endpoint | Statistic | " + " | ".join(label(r["cell"]) for r in admitted) + " |",
        "| --- | --- | ---: | ---: | ---: |"]
    for endpoint in ENDPOINTS:
        statistic = "F1" if endpoint in ENDPOINTS[:3] else ("Sample mean" if endpoint == "FAS" else "Similarity")
        text.append(f"| {endpoint} | {statistic} | " + " | ".join(f'{r["scores"][endpoint]:.8f}' for r in admitted) + " |")
    text.extend(["", "[Admitted native endpoint snapshot](native_qfo_scientific_scores_20261006_v2/report.json).\n"
        "The project-defined secondary means are " + ", ".join(f'{r["secondary_mean"]:.8f}' for r in admitted) +
        ", respectively, not official QfO F1 or a single accuracy ranking.\n"
        "Prediction coverage uses input accessions, not the number of assessed\n"
        "reference pairs; mapping loss is zero in these conversions.", "",
        "| Cell | Submitted pairs | Input-accession coverage | Output semantics |",
        "| --- | ---: | ---: | --- |"])
    for row in admitted:
        semantics = "Inferred ortholog pairs" if row["cell"].endswith("r1") else "Group-clique pairs"
        text.append(f'| {label(row["cell"])} | {row["submitted_pairs"]:,} | '
                    f'{100 * row["relation_coverage"]:.4f}% | {semantics} |')
    text.extend(["", "Two of 14 planned contrasts have exactly matched native SwissTrees family\n"
        "records; the other 12 remain unavailable. The original 100,000 draws\n"
        "(seed 20260922) and all 42 adjusted endpoints are reused, with no new draws.\n"
        "Differences and adjusted intervals below are percentage points; wins/ties/\n"
        "losses count the 18 retained families, not independent gene pairs.", "",
        "| Contrast | Endpoint | Difference (pp) | Adjusted interval (pp) | Wins/ties/losses |",
        "| --- | --- | ---: | --- | --- |"])
    for contrast in matched:
        for key, name in (("F1", "F1"), ("PPV", "Precision"), ("TPR", "Recall")):
            metric = contrast["metrics"][key]
            low, high = metric["bonferroni_percentile_ci"]
            counts = "/".join(str(metric[k]) for k in ("family_wins", "family_ties", "family_losses"))
            text.append(f'| {contrast["name"].replace("_", " ")} | {name} | '
                        f'{100 * metric["difference"]:+.4f} | [{100 * low:.4f}, {100 * high:.4f}] | {counts} |')
    text.extend(["", "Both adjusted F1 intervals include zero. Candidate expansion raises recall\n"
        "while lowering precision; reconciliation shifts that trade-off in the\n"
        "opposite direction. Count-based arithmetic differs slightly from native\n"
        "serialized endpoints and remains separately reported, not substituted.\n"
        "Development exposure, family exchangeability and percentile-coverage\n"
        "assumptions limit these conditional intervals; they establish neither\n"
        "equivalence nor independent confirmation or other-endpoint uncertainty.\n"
        "[Guarded interval binding](native_qfo_candidate_swiss_uncertainty_20261006_v1.json),\n"
        "[independent arithmetic readback](native_qfo_candidate_swiss_readback_20261006_v1.json).", "",
        "**Native QfO Figure.** The [reviewed six-panel figure](native_qfo_three_cell_figure_20261006_v1/native_qfo_three_cell.pdf)\n"
        "separates F1 endpoints from GO/EC/FAS, shows prediction coverage and\n"
        "SwissTrees precision/recall, and retains both adjusted contrast intervals.\n"
        "This is a three-cell ablation, not a complete fresh factorial or a timing\n"
        "figure. Original-Python numerical replay and visual review retain their\n"
        "own scopes. [Figure evidence](NATIVE_QFO_THREE_CELL_FIGURE_RESULT_20261006.md).", "",
        "The earlier [two-cell snapshot](native_qfo_scientific_scores_20261006_v1/report.json),\n"
        "[two-cell binding](native_qfo_swiss_uncertainty_binding_22449_20261006.json),\n"
        "[its arithmetic readback](recovered_native_qfo_swiss_readback_22449_20261006.json),\n"
        "[four-panel figure](native_qfo_p0c0_figure_20261006_v1/native_qfo_p0c0.pdf) and\n"
        "[figure review](NATIVE_QFO_FIGURE_RESULT_20261006.md) remain historical;\n"
        "they are not relabeled as this three-cell snapshot.", ""])
    tp, fp = transitions[("FN", "TP")], transitions[("not_scored", "FP")]
    text.append(f"Candidate expansion recovers {tp:,} asserted VGNC TPs and adds {fp:,} scored\n"
        f"FPs while retaining baseline TPs/FPs. Across {vgnc['union_scored_pairs']:,} union rows,\n"
        f"precision changes {100 * vgnc['differences']['precision']:+.6f} percentage points,\n"
        f"recall {100 * vgnc['differences']['recall']:+.6f} and F1 {100 * vgnc['differences']['f1']:+.6f}.\n"
        "All added FPs cross fixed reference overlap blocks. Those blocks are not\n"
        "validated independent uncertainty units; an unscored baseline pair is\n"
        "not a true negative. This explains the score arithmetic, not biological\n"
        "orthology or causal false-homology mechanisms.\n"
        "[VGNC decomposition](native_qfo_candidate_vgnc_20261006_v1/report.json),\n"
        "[independent readback](native_qfo_candidate_vgnc_readback_20261006_v1.json).\n")
    text.append(f"The complete accepted-union reconstruction preserves {aliases['genes']:,} genes\n"
        f"and transforms {aliases['baseline_groups']:,} baseline groups into {aliases['candidate_groups']:,}\n"
        f"candidate groups through {aliases['accepted_merges']:,} unions. The initial export\n"
        "failure is retained; a separately versioned follow-up proves the two\n"
        "scored/native accession bridges from the original protein map and DB\n"
        "rows, without suffix guessing, then independently localizes every\n"
        f"one of the {aliases['changed_pairs']:,} changed scored pairs.\n")
    text.extend(["| Candidate round | Group-connection path | Recovered TPs | Added FPs |",
                 "| ---: | --- | ---: | ---: |"])
    for iteration in (0, 1):
        for path, name in (("direct_cross_endpoint", "Direct group endpoints"), ("transitive_union", "Transitive union")):
            text.append(f"| {iteration} | {name} | {paths[(iteration, path, 'TP')]:,} | {paths[(iteration, path, 'FP')]:,} |")
    text.extend([f"| Total | All localized paths | {tp:,} | {fp:,} |", "",
        "Most added FPs have direct group-endpoint connections, ruling out an\n"
        "all-transitive software-path explanation. Directness here does not prove\n"
        "a direct protein-pair HMM hit, causal erroneous evidence or a biological\n"
        "mechanism; transitivity does not prove absence of supporting evidence.\n"
        "Original scored identifiers and categories remain unchanged.\n"
        "[Complete localization](native_qfo_candidate_alias_group_20261006_v1/report.json),\n"
        "[full pair ledger](native_qfo_candidate_alias_group_20261006_v1/changed_pair_groups.tsv),\n"
        "[independent graph readback](native_qfo_candidate_alias_group_readback_20261006_v1.json),\n"
        "[alias proof and retained failure](NATIVE_QFO_CANDIDATE_ALIAS_GROUP_RESULT_20261006.md).", ""])
    return "\n".join(text) + "\n"


def compose(original, scores, swiss, vgnc, aliases, reader):
    admitted, matched, paths, transitions = validate(scores, swiss, vgnc, aliases, reader)
    before, native, after = bounded(original, NATIVE_START, NATIVE_END)
    if native.count(FUNCTIONAL_START) != 1:
        raise ValueError("Functional diagnostic boundary missing")
    functional = FUNCTIONAL_START + native.split(FUNCTIONAL_START, 1)[1]
    current = before + native_section(admitted, matched, paths, transitions, vgnc, aliases) + functional + NATIVE_END + after
    heading = "This evidence-integrated working revision preserves the\n"
    if current.count(heading) != 1:
        raise ValueError("Chronology boundary ambiguous")
    left, rest = current.split(heading, 1)
    title = left.split("\n\n", 1)[0]
    current = title + "\n\nCondensed scientific draft, second native-evidence revision of 6 October 2026. Not submission-ready.\n"
    current += "This provisional snapshot preserves the [first 6 October text](PUBLICATION_MAIN_TEXT_20261006.md)\n"
    current += "and the [third 4 October text](PUBLICATION_MAIN_TEXT_20261004_v3.md),\n"
    current += "with all earlier rendered PDFs and archives. Its three-cell native QfO\n"
    current += "figure and complete candidate-pair paths are later additions, not retroactive\n"
    current += "members of those artifacts. Native inference and scientific/reproducibility\n"
    current += "requirements remain incomplete; dated receipts govern render/package scope.\n"
    current += "The prior revision's chronology remains recorded below.\n\n" + heading + rest
    methods = "### Shared-Host Resource Measurement\n"
    if current.count(methods) != 1:
        raise ValueError("Methods boundary ambiguous")
    new_methods = "A separately reviewed candidate diagnostic replays all accepted group unions\n"
    new_methods += "and joins changed scored pairs through the original protein map and scoped\n"
    new_methods += "original database rows. An independent standard-library graph traversal\n"
    new_methods += "checks whole partitions and all pair paths without importing the primary\n"
    new_methods += "union implementation. Group directness and transitivity describe software\n"
    new_methods += "connectivity, not calibrated homology confidence or causal biology. The\n"
    new_methods += "initial failed join remains retained, not rewritten as successful.\n"
    new_methods += "[Identity/path protocol](NATIVE_QFO_CANDIDATE_ALIAS_GROUP_PROTOCOL_20261006.md).\n\n"
    current = current.replace(methods, new_methods + methods, 1)
    before, _, after = bounded(current, AVAILABILITY_START, "## References\n")
    availability = "The second revision's 34-page review/rc3 archive, third revision's\n"
    availability += "35-page review/rc4 archive and first 6 October direct-review component\n"
    availability += "remain historical. They do not contain this later three-cell/complete-path\n"
    availability += "source snapshot. Fresh native tables, bounded interval bindings, existing\n"
    availability += "figures and pair diagnostics retain their source identities and failures.\n"
    availability += "This source is [mechanically generated](../prepare_native_main_text_v2.py)\n"
    availability += "from already-reviewed summaries, without new scoring, bootstrap draws,\n"
    availability += "admission, inference, default selection or independent-validation claims.\n"
    availability += "Any new local render/review or direct-link package includes only the\n"
    availability += "content and checks in its own dated receipt. Such reporting portability\n"
    availability += "does not supply transitive raw/dependency closure, a hermetic full-study\n"
    availability += "reproduction, redistribution clearance, public deposition or an archival DOI.\n"
    availability += "Existing replay receipts remain valid only within their recorded scopes.\n\n"
    return before + availability + "## References\n" + after


def build(output, receipt, root=ROOT):
    root = Path(root).resolve()
    output, receipt = Path(output).absolute(), Path(receipt).absolute()
    for path in (output, receipt):
        if path.exists() or path.is_symlink():
            raise FileExistsError(path)
        if path.parent.resolve() != root / BASE:
            raise ValueError("New manuscript and receipt must stay beside original local links")
    if output == receipt:
        raise ValueError("Distinct manuscript and receipt required")
    inputs, values = {}, {}
    for key, (name, expected) in PINS.items():
        path = root / BASE / name
        data = path.read_bytes()
        if hashlib.sha256(data).hexdigest() != expected:
            raise ValueError(f"Bound reporting input changed: {key}")
        inputs[key] = record(path)
        values[key] = data.decode("utf-8") if key == "original" else (
            json.loads(data) if name.endswith(".json") else None)
    source = record(__file__)
    text = compose(values["original"], values["scores"], values["swiss"], values["vgnc"],
                   values["aliases"], values["aliases_readback"])
    for ref in [*inputs.values(), source]:
        check(ref)
    with output.open("x", encoding="utf-8") as stream:
        stream.write(text)
    result = dict(schema="native_main_text_summary_generation_v1", status="provisional_source_generated",
        source=source, inputs=inputs, output=record(output), publication_ready=False,
        new_scoring_or_admission=False, new_bootstrap_draws=0, native_inference_reexecuted=False,
        original_evidence_modified=False, admitted_cells=list(CELLS), unavailable_cells=4,
        matched_swiss_contrasts=2, unavailable_swiss_contrasts=12,
        limitations=["Bound-summary presentation, not another scientific or raw-evidence audit.",
            "No render, visual review, archive restoration, transitive closure or rights implied.",
            "Fresh native factorial, independent-family validation and other-endpoint uncertainty remain incomplete."])
    with receipt.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--receipt", required=True, type=Path)
    args = parser.parse_args()
    result = build(args.output, args.receipt)
    print(json.dumps({k: result[k] for k in ("status", "output", "publication_ready")}))
