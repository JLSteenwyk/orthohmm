"""Generate a bounded native QfO supplement from existing verified reports."""

import argparse
import csv
import json
import math
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check

CELLS = ("p0_c0_r0", "p0_c0_r1")
ENDPOINTS = ("VGNC", "SwissTrees", "TreeFam-A", "GO", "EC", "FAS")
SCHEMAS = dict(
    snapshot="native_qfo_scientific_reporting_snapshot_v1",
    intervals="native_qfo_retained_swiss_uncertainty_binding_v1",
    functional="native_qfo_functional_pair_composition_v1",
    functional_reader="native_qfo_functional_pair_sql_readback_v1",
    transitions="native_qfo_swiss_pair_transitions_v1",
    transitions_reader="native_qfo_swiss_transition_sql_readback_v1",
    reconciliation="native_qfo_swiss_reconciliation_localization_v1",
    reconciliation_reader="native_qfo_swiss_reconciliation_newick_readback_v1",
    search="native_qfo_swiss_direct_search_support_v1",
    search_reader="native_qfo_swiss_search_support_code_readback_v1",
    graph="native_qfo_swiss_graph_support_v1",
    graph_reader="native_qfo_swiss_graph_support_floyd_readback_v1",
    strata="native_qfo_swiss_sequence_strata_v1",
    strata_reader="native_qfo_swiss_sequence_strata_readback_v1",
    vgnc="native_qfo_vgnc_blocks_v1", vgnc_reader="native_qfo_vgnc_blocks_readback_v1",
    score_figure="native_qfo_p0c0_figure_v1", score_figure_reader="native_qfo_figure_content_readback_v1",
    mechanism_figure="native_qfo_swiss_mechanism_figure_v1",
    mechanism_figure_reader="native_qfo_swiss_mechanism_figure_readback_v1")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def validate(data, refs):
    require(set(data) == set(refs) == set(SCHEMAS), "Incomplete selected evidence")
    for key, schema in SCHEMAS.items():
        require(data[key].get("schema") == schema and data[key].get("publication_ready") is False,
            "Changed evidence schema/scope: " + key)
    for primary in ("functional", "transitions", "reconciliation", "search", "graph", "strata", "vgnc",
                    "score_figure", "mechanism_figure"):
        field = "composition" if primary == "functional" else "manifest" if primary.endswith("figure") else "report"
        require(data[primary + "_reader"][field] == refs[primary], "Readback binding differs: " + primary)
    for key in ("intervals", "functional", "transitions", "strata", "vgnc", "score_figure"):
        require(data[key]["snapshot"] == refs["snapshot"], "Mixed scientific snapshots")
    require(data["reconciliation"]["transition"] == refs["transitions"]
        and data["search"]["localization"] == refs["reconciliation"]
        and data["graph"]["search"] == refs["search"]
        and data["graph"]["search_readback"] == refs["search_reader"]
        and data["strata"]["binding"] == refs["intervals"]
        and data["score_figure"]["swiss_binding"] == refs["intervals"], "Mixed diagnostic chain")
    for primary in ("reconciliation", "search", "strata"):
        require(data["mechanism_figure"]["inputs"][primary] == [refs[primary], refs[primary + "_reader"]],
            "Mixed mechanism figure inputs")
    snapshot = data["snapshot"]
    require(snapshot.get("new_scoring_or_admission") is False
        and snapshot.get("recovered_inference_resources_admitted") is False, "Changed snapshot admission scope")
    rows = snapshot["rows"]
    require(len(rows) == 7 and [r["index"] for r in rows] == list(range(6, 13)), "Changed planned identities")
    admitted = [r for r in rows if r.get("accuracy_admitted") is True]
    require([r["cell"] for r in admitted] == list(CELLS), "Require two completed native cells only")
    require(admitted[1].get("timing_eligible") is False and admitted[1].get("timing_admitted") is False
        and admitted[1].get("resources", "missing") is None, "Recovered timing relabeled")
    for row in admitted:
        require(set(row["scores"]) == set(ENDPOINTS)
            and all(type(v) in (int, float) and math.isfinite(v) and 0 <= v <= 1 for v in row["scores"].values()),
            "Invalid score inventory")
        require(math.isclose(sum(row["scores"].values()) / 6, row["secondary_mean"], abs_tol=1e-12, rel_tol=0),
            "Secondary mean arithmetic differs")
    for row in rows[2:]:
        require(row.get("accuracy_admitted") is False and all(v is None for v in row["scores"].values()),
            "Missing identity has a supplied score")
    intervals = data["intervals"]
    require(type(intervals["new_bootstrap_draws"]) is int and intervals["new_bootstrap_draws"] == 0
        and intervals["replicates_reused"] == 100000
        and intervals["seed_reused"] == 20260922 and intervals["multiplicity_endpoints"] == 42
        and len(intervals["families"]) == 18, "Changed retained uncertainty protocol")
    matches = [c for c in intervals["contrasts"] if c["name"] == "R_at_P0_C0"]
    require(len(matches) == 1 and matches[0]["status"] == "native_records_matched", "Missing matched contrast")
    for report in (data["vgnc"], data["graph"], data["functional"], data["transitions"]):
        require(report["uncertainty_admitted"] is False and report["new_scoring_or_admission"] is False,
            "Diagnostic promoted to scientific admission")
    require(data["graph"]["graph_original_admission_established"] is False
        and data["reconciliation"]["node_annotations_previously_inventoried"] is False,
        "Observed artifacts promoted to original admission")
    require([m["cell"] for m in data["vgnc"]["methods"]] == list(CELLS), "Changed VGNC cohort")
    for method, row in zip(data["vgnc"]["methods"], admitted):
        require(math.isclose(method["metrics"]["f1"], row["scores"]["VGNC"], abs_tol=1e-12, rel_tol=0),
            "VGNC and snapshot differ")
    changed = data["transitions"]["comparison"]["changed_relations"]
    require(data["graph"]["changed_pairs"] == data["search"]["changed_pairs"]
        == data["reconciliation"]["changed_pairs_traced"] == changed, "Changed diagnostic cohort")
    for cell in CELLS:
        require(sum(r["pairs"] for r in data["graph"]["summary"] if r["cell"] == cell) == changed,
            "Graph summary lost pairs")
    require(not any(r["pairs"] for r in data["graph"]["summary"] if r["graph_support"] == "disconnected"),
        "Observed connectivity statement no longer holds")
    return admitted, matches[0]


def tables(data, admitted, contrast):
    scores = [[endpoint, "F1" if endpoint in ENDPOINTS[:3] else "Similarity" if endpoint in ("GO", "EC") else "FAS",
        *[r["scores"][endpoint] for r in admitted]] for endpoint in ENDPOINTS]
    scores.append(["Six-metric mean", "Secondary, not F1", *[r["secondary_mean"] for r in admitted]])
    intervals = [[metric, 100 * item["difference"], *[100 * v for v in item["bonferroni_percentile_ci"]],
        item["family_wins"], item["family_ties"], item["family_losses"]]
        for metric, item in contrast["metrics"].items()]
    vgnc = [[m["cell"], *[m["counts"][c] for c in ("TP", "FP", "FN")],
        *[m["metrics"][c] for c in ("precision", "recall", "f1")]] for m in data["vgnc"]["methods"]]
    vg_changes = [[r["r0"], r["r1"], r["pairs"]] for r in data["vgnc"]["transition_counts"]]
    swiss_changes = [[*key.split("->"), value] for key, value in data["transitions"]["comparison"]["transitions"].items()]
    graph = [[cell, before, *[sum(r["pairs"] for r in data["graph"]["summary"]
        if r["cell"] == cell and r["before"] == before and r["graph_support"] == kind)
        for kind in ("direct_edge", "indirect_path", "disconnected")]] for cell in CELLS for before in ("TP", "FP")]
    strata = [[r["stratum"], r["families"], *[None if r[k] is None else 100 * r[k] for k in ("F1", "PPV", "TPR")]]
        for r in data["strata"]["differences"]]
    functional = []
    for r in data["functional"]["comparisons"]:
        item, metric = r["result"], r["metric"]
        suffix = "sample_pairs" if metric == "FAS" else "pairs"
        functional.append([metric, item["right_" + suffix], item["left_" + suffix], item["shared_" + suffix],
            item["shared_pairs_with_different_serialized_scores"]])
    return dict(scores=(["Endpoint", "Statistic", "R0", "R1"], scores),
        swiss_intervals=(["Metric", "Difference_pp", "Adjusted_lower_pp", "Adjusted_upper_pp", "Wins", "Ties", "Losses"], intervals),
        vgnc_counts=(["Cell", "TP", "FP", "FN", "Precision", "Recall", "F1"], vgnc),
        vgnc_transitions=(["R0", "R1", "Pairs"], vg_changes),
        swiss_transitions=(["R0", "R1", "Pairs"], swiss_changes),
        graph_support=(["Cell", "Removed", "Direct_edge", "Indirect_path", "Disconnected"], graph),
        strata_differences=(["Stratum", "Families", "F1_pp", "Precision_pp", "Recall_pp"], strata),
        functional_overlap=(["Endpoint", "R0_scored_pairs", "R1_scored_pairs", "Shared", "Changed_shared_scores"], functional))


def markdown_table(header, rows):
    def value(v):
        return "Unavailable" if v is None else f"{v:.6f}" if isinstance(v, float) else str(v)

    return "\n".join(["| " + " | ".join(header) + " |", "| " + " | ".join(["---"] * len(header)) + " |",
        *["| " + " | ".join(value(v) for v in row) + " |" for row in rows]])


def document(data, refs, output, admitted, contrast, derived, figure_assets):
    def link(key, title):
        return f"[{title}]({Path(os.path.relpath(refs[key]['path'], output)).as_posix()})"

    def display(name, header=None, rows=None):
        original_header, original_rows = derived[name]
        return markdown_table(header or original_header, original_rows if rows is None else rows)

    missing = [r["cell"] for r in data["snapshot"]["rows"] if r.get("accuracy_admitted") is False]
    text = ["# OrthoHMM Native QfO Supplement", "",
        "Working manuscript supplement, 6 October 2026. Not submission-ready. Generated from existing verified reports; "
        "no new inference, scoring, accuracy admission, bootstrap or timing repair. Historical main text and archive bytes "
        "are unchanged. This document supplements, not supersedes, their frozen methods and all-tool results.", "",
        "## Methods And Evidence Scope", "",
        "P0/C0/R0 and P0/C0/R1 retain initial HMM search; downstream profile refinement and candidate expansion are off. "
        "R changes group-derived cross-species clique predictions to native phylogenetically inferred pairs. This is not "
        "a total-HMM ablation, the selected-default high-sensitivity/satellite_v2 comparison or a sensitivity-matched "
        "non-HMM control. QfO and OrthoBench remain development-exposed primary benchmarks; Three Kingdoms is supplementary.", "",
        f"Only {len(admitted)} of {len(data['snapshot']['rows'])} native QfO identities have supplied accuracy admissions. "
        "Missing from this snapshot: " + ", ".join(f"`{cell}`" for cell in missing) + ". "
        "Unavailable cells are not imputed and their interactions cannot be evaluated here. R1 accuracy was recovered "
        "after a native measurement failure; its inference resource values remain null and timing ineligible. "
        "Current job states belong in the progress ledger, not this fixed score snapshot.", "",
        "Report and code identities are directly checked, with each diagnostic's original independent readback bound. "
        "Existing scientific receipts are reused; their raw scans and transitive admissions are not repeated for writing. "
        "Present-day checksum agreement is not evidence of uninterrupted integrity. " + link("snapshot", "Score snapshot") + ".", "",
        "## Native Endpoint Results", "", display("scores"), "",
        "VGNC, SwissTrees and TreeFam-A are F1 endpoints. GO/EC similarity and FAS are not F1. "
        "The six-metric arithmetic mean is project-defined and secondary, not official QfO F1 or a superiority ranking. "
        "Native serialized scores are retained rather than replaced by unrounded diagnostic arithmetic.", "",
        "## SwissTrees Conditional Uncertainty", "",
        "Paired intervals are attached only after every complete native family record matches the retained count manifest. "
        "The analysis reuses 100,000 family draws, seed 20260922, and correction over all 42 originally planned endpoints. "
        "It does not shrink multiplicity to the single available contrast or create new independent confirmation. "
        "The statistic is the harmonic mean of macro precision and recall, not mean family F1 or pooled pair F1.", "",
        display("swiss_intervals", ["Metric", "Change (pp)", "Adjusted Lower", "Adjusted Upper", "Wins", "Ties", "Losses"]), "",
        "The adjusted F1 interval includes zero: a conclusive improvement is not established. The intervals are conditional "
        "on 18 development-exposed reference families, exchangeability and approximate percentile coverage. "
        "They do not apply to VGNC, TreeFam-A, GO, EC, FAS or the secondary mean. " + link("intervals", "Retained interval binding") + ".", "",
        "## VGNC Precision-Recall Trade-Off", "", display("vgnc_counts"), "", display("vgnc_transitions"), "",
        "The actual native scored-pair union retains all category changes, not historical selected-default counts. "
        "`not_scored` is not a TN or biological non-orthology label. Shared-protein reference labels form overlap blocks; "
        "neither these blocks nor method-dependent prediction components establish independent biological units. "
        "Rare-error/shared-clade uncertainty-validation failures remain unresolved; no arbitrary block bootstrap is used. "
        "Full prediction-database hashes were checked before and after reference-mapping reads, but prediction edges were "
        "not requeried and omitted-FP completeness was not independently rescored. " + link("vgnc", "Native decomposition") + ".", "",
        "## SwissTrees Pair And Stage Localization", "",
        display("swiss_transitions", rows=[r for r in derived["swiss_transitions"][1] if r[2]]), "",
        "These pair counts describe all scored reference relations; they are not the macro F1 denominator. "
        "Every changed pair was traced, avoiding selection of representative successes. " + link("transitions", "Pair transitions") + ".", "",
        display("graph_support"), "",
        "A direct graph edge or within-family path is observed homology support, not true orthology, a causal recruitment "
        "test or evolutionary distance. Full candidate-family members, including same-species and non-reference genes, "
        "were eligible intermediates. All changed pairs without a direct significant hit had indirect graph paths; "
        "direct-hit absence alone does not isolate a prefilter failure. The combined final graph does not isolate "
        "reciprocal-best-normalized-hit and singleton-attachment decisions. " + link("search", "Search evidence") + "; "
        + link("graph", "Graph evidence") + ".", "",
        "All changed pairs share their original candidate family and have an exclusion LCA annotated as a duplication "
        "under the saved positive-paralogy rule. Some excluded true pairs lie within the same final root HOG: inspecting "
        "only split groups misses those exclusions. Saved Newicks independently reproduce LCAs, but inferred duplication "
        "annotations and node confidence are not validated biological histories. Graphs and node annotations were not "
        "inventoried by the original validators; later observed bytes cannot be promoted into original admission. "
        + link("reconciliation", "Reconciliation localization") + ".", "",
        "## Prespecified Sequence Strata", "",
        display("strata_differences", ["Input-Only Bin", "Families", "F1 Change (pp)", "Precision Change", "Recall Change"],
            [r for r in derived["strata_differences"][1] if r[0] in ("higher_entropy", "lower_entropy", "short_relative", "not_short_relative")]), "",
        "Unchanged September input-only bins describe associations, not validated subgroup effects or mechanisms. "
        "Global residue entropy is not divergence or local low complexity; relative shortness does not prove fragmentation. "
        "Empty bins remain unavailable in the complete generated table. No subgroup CI, tuning or default promotion follows. "
        + link("strata", "Complete native projections") + ".", "",
        "## Functional Scores And Pair Selection", "", display("functional_overlap"), "",
        "Every R1 GO/EC scored pair is present in R0 with the same six-decimal serialized score. Their higher native means "
        "therefore reflect changed scored-pair membership and denominators at retained precision, not improved values on "
        "shared pairs. This does not prove orthology correctness or full-precision equality. FAS rows are realized samples, "
        "not all eligible predictions; their small shared subset, unseeded method-specific mixtures and missing-score "
        "attrition do not identify inclusion probabilities or support a paired population CI. " + link("functional", "Functional composition") + ".", "",
        "## Verified Companion Figures", ""]
    captions = dict(score_figure="Native P0/C0 endpoint scores and the conditional SwissTrees reconciliation contrast. "
        "F1 uncertainty remains inconclusive; other endpoints are point estimates without validated difference intervals.",
        mechanism_figure="Complete SwissTrees exclusions, direct significant-search evidence and unchanged input-only strata. "
        "This earlier companion does not depict the subsequent graph/VGNC diagnostic; those are reported above.")
    for number, key in enumerate(("score_figure", "mechanism_figure"), 1):
        assets = figure_assets[key]
        relative = {kind: Path(os.path.relpath(ref["path"], output)).as_posix() for kind, ref in assets.items()}
        # Inline caption text prevents Pandoc's implicit figure caption duplication.
        text += [f"![Figure S{number}]({relative['png']})\n[Vector PDF]({relative['pdf']}). " + captions[key], ""]
    text += ["## Claim Boundaries And Remaining Work", "",
        "Supported: the two admitted native cells and their explicit precision-recall/selection differences; conditional "
        "SwissTrees interval reuse; complete observed error-stage localization. Not established: conclusive native F1 "
        "superiority, verified evolutionary causes, total-HMM advantage, isolated efficiency or genome-wide generalization.", "",
        "Timing measurements were collected on a shared Threadripper while other analyses were running. Competition for "
        "CPU, memory bandwidth and I/O may have affected elapsed times, with an unknown and potentially tool-dependent "
        "impact. These are observed shared-host timings, not estimates of isolated performance. Matching resource limits "
        "does not isolate speed; failed R1 timing remains ineligible. This supplement does not create a new timing panel.", "",
        "Remaining work includes the missing native identities and interactions, sensitivity-matched search, valid wider "
        "uncertainty, independent-generalization limitations, remaining error strata, original TreeFam files and all-tool "
        "provenance/resource gaps. The full manuscript, executable archive/release and external-deposition scope remain "
        "unfinished. No dedicated host or quiet window is required; safe capacity and valid accounting still apply.", "",
        "## Direct Evidence Index", "",
        *["- " + link(key, key.replace("_", " ")) for key in SCHEMAS], ""]
    return "\n".join(text)


def assemble(root, selection, selection_sha, output):
    root, output = Path(root).resolve(), Path(output).absolute()
    require(output.resolve().is_relative_to(root), "Output must remain inside repository")
    require(not output.exists() and not output.is_symlink(), "Output already exists")
    selection_ref = record(selection)
    require(selection_ref["sha256"] == selection_sha, "Changed selection")
    chosen = json.loads(Path(selection).read_text())
    require(chosen.get("schema") == "native_qfo_supplement_selection_v1"
        and chosen.get("publication_ready") is False and set(chosen["inputs"]) == set(SCHEMAS), "Changed selection scope")
    refs, data, checked = {}, {}, [selection_ref, record(__file__)]
    for key, ref in chosen["inputs"].items():
        relative = Path(ref["path"])
        require(not relative.is_absolute() and ".." not in relative.parts, "Unsafe input path")
        path = (root / relative).resolve()
        require(path.is_relative_to(root), "Input escapes repository")
        actual = record(path)
        require((actual["bytes"], actual["sha256"]) == (ref["bytes"], ref["sha256"]), "Changed selected input: " + key)
        refs[key], data[key] = actual, json.loads(path.read_text())
        checked.extend([actual, data[key]["source"]])
        check(data[key]["source"])
    admitted, contrast = validate(data, refs)
    derived = tables(data, admitted, contrast)
    figure_assets = {}
    for key in ("score_figure", "mechanism_figure"):
        figure_assets[key] = {}
        for suffix in ("png", "pdf"):
            matches = [r for r in data[key]["outputs"] if r["path"].endswith("." + suffix)]
            require(len(matches) == 1, "Ambiguous figure asset")
            asset = matches[0]
            require(Path(asset["path"]).resolve().is_relative_to(root), "Figure escapes repository")
            check(asset)
            figure_assets[key][suffix] = asset
            checked.append(asset)
    output.mkdir(parents=True)
    outputs = []
    for name, (header, rows) in derived.items():
        path = output / (name + ".tsv")
        with path.open("x", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
            writer.writerow(header)
            writer.writerows(["Unavailable" if v is None else v for v in row] for row in rows)
        outputs.append(record(path))
    manuscript = output / "supplement.md"
    manuscript.write_text(document(data, refs, output, admitted, contrast, derived, figure_assets))
    outputs.append(record(manuscript))
    for ref in checked:
        check(ref)
    report = dict(schema="native_qfo_supplement_assembly_v1", status="bounded_native_supplement_generated",
        source=record(__file__), selection=selection_ref, checked_records=checked, selected_inputs=refs,
        outputs=outputs, table_rows={k: len(v[1]) for k, v in derived.items()},
        admitted_cells=list(CELLS), missing_cells=[r["cell"] for r in data["snapshot"]["rows"] if not r["accuracy_admitted"]],
        scientific_evidence_replayed=False, new_scoring_or_admission=False, new_bootstrap_draws=False,
        original_manuscript_or_archive_changed=False, recovered_timing_admitted=False, publication_ready=False,
        limitations=["Direct report/source/readback bindings, not transitive raw scientific re-admission.",
            "Presentation of development-exposed evidence; wider uncertainty and missing native identities remain open.",
            "No rendered PDF or visual review is implied by Markdown generation.",
            "Local absolute evidence paths and original runtime assets remain requirements, not a portable full-study archive."])
    with (output / "assembly.json").open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--selection-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = assemble(args.root, args.selection, args.selection_sha256, args.output)
    print(json.dumps(dict(status=result["status"], table_rows=result["table_rows"])))
