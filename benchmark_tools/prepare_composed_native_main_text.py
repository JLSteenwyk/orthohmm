"""Integrate actual composed-native evidence into a separate manuscript revision."""

import argparse
import csv
import json
import math
import os
from pathlib import Path

from benchmark_tools import plot_composed_native_qfo_scores as figure
from benchmark_tools.export_native_factorial_progress import finite, load, require, same_bytes
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

PARENT_SHA = "0b7012e3dd46bd3b028857bdb3cbe3a9e9048936b90f3df9e4833cb342d5705a"
FAILURE_SHA = "4f1a5c536af8ce992d46866374b3ab44ee3158077316d69ba2e9f5fb34cee203"
START = "### Fresh Native QfO Trade-Off And Pair Composition\n\n"
MECHANISM = "Every changed assessed SwissTrees profile pair was traced, avoiding favorable\n"
END = "Candidate expansion recovers 162 asserted VGNC TPs and adds 2,133 scored\n"
ABSTRACT_OLD = ("Fresh native ablations retain neutral\n"
    "downstream-profile effects and a SwissTrees precision-recall trade-off whose\n"
    "adjusted F1 interval includes zero.")


def link(title, ref, directory):
    relative = Path(os.path.relpath(ref["path"], directory)).as_posix()
    return f"[{title}]({relative})"


def inputs(specifications):
    require(set(specifications) == {"parent", "snapshot", "intervals", "figure", "reader", "failure"},
        "Require all six manuscript inputs")
    docs, refs, evidence = {}, {}, []
    for key, (path, digest) in specifications.items():
        if key == "parent":
            ref = record(path)
            require(ref["sha256"] == digest == PARENT_SHA, "Changed frozen v4 manuscript")
            docs[key], refs[key] = Path(path).read_text(), ref
            evidence.append(ref)
        else:
            docs[key], refs[key] = load(path, digest, evidence)
    require(refs["failure"]["sha256"] == FAILURE_SHA, "Changed retained scoring-failure addendum")
    return docs, refs, evidence


def validate(docs, refs):
    rows, scores, coverage, statuses, intervals = figure.figure_data(docs["snapshot"], docs["intervals"])
    sources = dict(snapshot="export_composed_native_qfo_scientific_scores.py",
        intervals="bootstrap_composed_native_qfo_swiss.py", figure="plot_composed_native_qfo_scores.py",
        reader="review_composed_native_qfo_figure.py", failure="export_composed_qfo_failure_addendum.py")
    for key, name in sources.items():
        require(docs[key]["source"] == record(Path(__file__).with_name(name)), "Changed input source: " + key)
    snapshot, uncertainty, plotted, reader, failure = (docs[key] for key in (
        "snapshot", "intervals", "figure", "reader", "failure"))
    require(snapshot["failure_addendum"] == refs["failure"]
        and uncertainty["snapshot"] == plotted["snapshot"] == reader["snapshot"] == refs["snapshot"]
        and plotted["native_intervals"] == reader["native_intervals"] == refs["intervals"]
        and reader["manifest"] == refs["figure"], "Cross-artifact manuscript binding")
    require(plotted["schema"] == "composed_native_qfo_figure_v1" and plotted["plotted_cells"] == list(figure.CELLS)
        and plotted["plotted_score_endpoints"] == 30 and plotted["plotted_coverage_endpoints"] == 5
        and plotted["plotted_swiss_contrast_endpoints"] == len(intervals)
        and plotted["new_bootstrap_draws"] == 0 and plotted["consumed_native_bootstrap_draws"] == 100000
        and plotted["publication_ready"] is False and plotted["new_scoring_or_admission"] is False,
        "Changed figure scope")
    require(reader["schema"] == "composed_native_qfo_figure_readback_v1"
        and reader["score_rows"] == 42 and reader["coverage_rows"] == reader["status_rows"] == 7
        and reader["interval_rows"] == len(intervals) and reader["scientific_tables_matched"] is True
        and reader["assets_decoded"] is True and reader["publication_ready"] is False
        and reader["new_scoring_or_admission"] is False, "Incomplete independent figure readback")
    require(failure["schema"] == "composed_native_qfo_failure_reporting_v1"
        and failure["publication_ready"] is False and failure["new_scientific_admission"] is False
        and same_bytes(failure["references"]["manuscript"], refs["parent"]), "Changed retained failure scope")
    failed, row = failure["row"], rows["p1_c1_r0"]
    require(failed["assessment_job_id"] == 24038 and failed["native_index"] == 11
        and failed["admitted_endpoint_score_count"] == 0 and failed["accuracy_admitted"] is False
        and failed["secondary_six_metric_mean"] is None and failed["scoring_status"] == "OUT_OF_MEMORY"
        and failed["failed_endpoint"] == "FAS" and failed["failed_endpoint_exit_code"] == 137
        and failed["scoring_peak_memory_bytes"] is None and failed["scoring_cpu_slots"] == 8
        and failed["scoring_memory_limit_bytes"] == 32*1024**3
        and all(row[key] == failed[key] for key in ("input_accessions", "relation_accessions",
            "submitted_pairs", "relation_coverage", "prediction_semantics")), "Failed scoring contributes no score")
    for cell in figure.CELLS:
        row = rows[cell]
        mean = finite(row["secondary_mean"], "secondary mean", high=1)
        require(math.isclose(sum(row["scores"].values())/6, mean, rel_tol=0, abs_tol=1e-12),
            "Secondary mean arithmetic differs")
    return rows, [row for row in uncertainty["comparisons"] if row["metrics"] is not None]


def interpretation(effects):
    result = dict(positive=[], negative=[], includes_zero=[])
    for effect in effects:
        lower, upper = effect["metrics"]["F1"]["bonferroni_percentile_ci"]
        key = "positive" if lower > 0 else "negative" if upper < 0 else "includes_zero"
        result[key].append(effect["name"])
    return result


def claims(docs, refs, directory):
    rows, effects = validate(docs, refs)
    return [dict(claim="Native endpoint comparison", scope="Five admitted development-exposed cells; two scores missing",
        evidence=link("Native snapshot", refs["snapshot"], directory), excluded="General superiority; independent validation"),
        dict(claim="Conditional SwissTrees uncertainty", scope=f"{len(effects)}/14 planned contrasts; 42 adjusted endpoints",
            evidence=link("Actual native-count intervals", refs["intervals"], directory),
            excluded="Missing-cell imputation; other-challenge or secondary-mean intervals"),
        dict(claim="Scoring failure retention", scope="Native11 conversion coverage only; FAS137/OUT_OF_MEMORY; zero admitted endpoints",
            evidence=link("Retained failure", refs["failure"], directory), excluded="Partial score or zero-filled mean; automatic retry"),
        dict(claim="Figure arithmetic and assets", scope="42 score rows, seven status/coverage rows; independent table and asset readback",
            evidence=link("Figure readback", refs["reader"], directory), excluded="New scoring; manuscript render review; full-study archive"),
        dict(claim="Resource comparison", scope="Shared-host observations; inference and scoring resources remain separate",
            evidence=link("Native provenance", refs["snapshot"], directory), excluded="Isolated speed ranking; repaired failed timing")]


def section(docs, refs, directory, claims_path):
    rows, effects = validate(docs, refs)
    lines = [START.rstrip(), "", "Five of seven fresh corrected-QfO cells have admitted accuracy; two remain unavailable.",
        "P denotes downstream profile refinement, C candidate expansion and R phylogenetic",
        "pair inference. Initial sensitive HMM search remains on in every cell. This is",
        "a development-exposed component analysis, not the selected-default all-tool",
        "comparison above, a total-HMM-off control or a complete fresh factorial.", "",
        "| Endpoint | Statistic | P0/C0/R0 | P0/C0/R1 | P0/C1/R0 | P1/C0/R1 | P1/C1/R1 |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: |"]
    for endpoint in figure.ENDPOINTS:
        statistic = "F1" if endpoint in figure.ENDPOINTS[:3] else "Sample mean" if endpoint == "FAS" else "Similarity"
        lines.append("| " + " | ".join([endpoint, statistic,
            *[f"{rows[cell]['scores'][endpoint]:.8f}" for cell in figure.CELLS]]) + " |")
    lines.extend(["", link("Admitted endpoint snapshot", refs["snapshot"], directory)+".",
        "The project-defined secondary means are "+", ".join(f"{rows[cell]['secondary_mean']:.8f}" for cell in figure.CELLS)
        +" respectively; these are not official QfO F1 or a single accuracy ranking.", "",
        "| Cell | Accuracy status | Submitted pairs | All-input relation coverage | Output semantics |",
        "| --- | --- | ---: | ---: | --- |"])
    for row in docs["snapshot"]["rows"]:
        count = "Unavailable" if row["submitted_pairs"] is None else f"{row['submitted_pairs']:,}"
        coverage = "Unavailable" if row["relation_coverage"] is None else f"{100*row['relation_coverage']:.4f}%"
        lines.append("| "+" | ".join([row["cell"], row["status"], count, coverage, row["prediction_semantics"]])+" |")
    lines.extend(["", "Coverage is not accuracy and uses all input accessions, not reference pairs.",
        "Native11 completed inference and group-clique conversion, but its six-endpoint",
        "QfO assessment failed: FAS exited137 and the scheduler reported OUT_OF_MEMORY",
        "with 8 CPU slots/32 GiB. Five completed endpoint tasks supply zero admitted scores",
        "and no secondary mean. Its scoring peak memory is unknown, not measured zero.",
        link("Retained scoring failure", refs["failure"], directory)+".",
        "Final native12 scoring uses a separately prespecified 8 CPU/128 GiB envelope.",
        "That scoring allocation is not part of the matched-resource inference timing",
        "panel or evidence of isolated efficiency. Recovered P0/C0/R1 accuracy retains",
        "failed timing, null resource measurements and timing ineligibility.", "",
        f"{len(effects)} of 14 planned SwissTrees contrasts are estimable from actual admitted",
        "native family counts. This analysis uses 100,000 new shared paired-family draws",
        "(seed 20260922), not cached interval substitution. All 42 planned metric endpoints",
        "remain in the adjustment; unavailable contrasts retain null metrics. F1 is",
        "the harmonic mean of family-mean precision and recall, recomputed within each",
        "replicate, not mean family F1. Effects and intervals below are percentage points.", "",
        "| Contrast | Endpoint | Difference (pp) | Adjusted interval (pp) | Wins/ties/losses |",
        "| --- | --- | ---: | --- | --- |"])
    for effect in effects:
        for metric in ("F1", "PPV", "TPR"):
            value = effect["metrics"][metric]
            interval = "["+", ".join(f"{100*x:.4f}" for x in value["bonferroni_percentile_ci"])+"]"
            signs = "/".join(str(value[key]) for key in ("family_wins", "family_ties", "family_losses"))
            lines.append(f"| {effect['name']} | {metric} | {100*value['difference']:+.4f} | {interval} | {signs} |")
    summary = interpretation(effects)
    lines.extend(["", "Adjusted F1 intervals: "+", ".join(f"{len(summary[key])} {label}" for key, label in (
        ("positive", "strictly above zero"), ("negative", "strictly below zero"), ("includes_zero", "include zero")))+".",
        "These are conditional pipeline differences, not identified causal HMM effects,",
        "equivalence, family-disjoint confirmation or superiority over OrthoFinder.",
        "The 18 development-exposed families retain exchangeability, selection and",
        "approximate percentile-coverage limitations. No intervals follow for other",
        "QfO challenges or the secondary mean.", link("Actual native intervals", refs["intervals"], directory)+".", "",
        "**Native QfO Figure.** "+link("Endpoint, coverage and conditional-interval figure",
            next(ref for ref in docs["figure"]["outputs"] if ref["path"].endswith(".pdf")), directory),
        "separates orthology F1 from GO/EC/FAS similarities and coverage. Its 42 score",
        f"rows, seven coverage/status rows and {3*len(effects)} interval rows were independently",
        "read against the bound source records, and PDF/PNG assets decoded. A machine",
        "asset check does not establish manual visual review or a new manuscript",
        "render/archive receipt. "+link("Independent figure readback", refs["reader"], directory)+".",
        "The earlier four-cell figure, intervals and pair-localization reports remain",
        "historical and unchanged. The following localization concerns only the original",
        "P0/C0/R1 versus P1/C0/R1 contrast, not a new mechanistic diagnosis of native12.",
        link("Native claim-to-evidence scope", {"path": str(claims_path)}, directory)+".", ""])
    return "\n".join(lines)+"\n"


def manuscript(docs, refs, directory, claims_path):
    parent = docs["parent"]
    require(all(parent.count(anchor) == 1 for anchor in ("## Abstract\n", START, MECHANISM, END, ABSTRACT_OLD)),
        "Missing or duplicate manuscript anchors")
    body = parent[parent.index("## Abstract\n"):]
    a, b = body.index(START), body.index(MECHANISM)
    body = body[:a]+section(docs, refs, directory, claims_path)+body[b:]
    abstract = ("Fresh native QfO ablations admit five of seven cells; two scores remain\n"
        "unavailable. Conditional SwissTrees intervals use actual native-family counts\n"
        "with 42-endpoint adjustment; missing cells do not supply imputed contrasts.")
    body = body.replace(ABSTRACT_OLD, abstract, 1)
    replacements = (
        ("They do not contain this later four-cell native QfO source snapshot.",
         "They do not contain this later composed-native five-cell QfO source snapshot."),
        ("This source is [mechanically extended](../prepare_four_cell_main_text.py)\n",
         "This source is [mechanically extended](../prepare_composed_native_main_text.py)\n"),
    )
    for old, new in replacements:
        require(body.count(old) == 1, "Changed availability anchor")
        body = body.replace(old, new, 1)
    old = "A read-only functional diagnostic binds six native raw GO/EC/FAS tables through\n"
    require(body.count(old) == 1, "Changed native-protocol anchor")
    body = body.replace(old, "A separate prespecified native-count analysis recomputes supported QfO\n"
        "contrasts from actual admitted family counts using 100,000 shared draws, the\n"
        "same seed 20260922 and all 42 adjusted endpoints. Missing native cells are not\n"
        "substituted with cached predictions. This is conditional uncertainty, not\n"
        "independent confirmation.\n\n"+old, 1)
    header = ("# OrthoHMM: HMM-Centered Group Inference With Phylogenetic Refinement\n\n"
        "Composed-native QfO evidence revision, v5. Not submission-ready.\n"
        "Preserves "+link("the frozen reviewed v4 source", refs["parent"], directory)+" and all prior artifacts.\n"
        "Integrates five admitted cells, explicit missing outcomes, actual native-count\n"
        "intervals and independent figure readback. All-tool comparisons, frozen\n"
        "methods, simulations, transfer tests, prior diagnostics and citations retain\n"
        "their original scope. This generator performs no scoring, bootstrap draws,\n"
        "inference, default selection or independent confirmation. The consumed\n"
        "native-count analysis has 100,000 new draws, distinguished from this reporting.\n"
        "New manuscript render review and archive inclusion require separate evidence.\n\n")
    return header+body


def generate(directory, specifications, output, receipt, claims_path):
    paths = (output, receipt, claims_path)
    require(len({path.resolve() for path in paths}) == 3
        and all(path.parent.resolve() == directory.resolve() for path in paths)
        and all(not path.exists() and not path.is_symlink() for path in paths), "Require distinct fresh outputs beside evidence")
    docs, refs, evidence = inputs(specifications)
    text, claim_rows = manuscript(docs, refs, directory, claims_path), claims(docs, refs, directory)
    for key in ("snapshot", "intervals", "figure"):
        evidence.extend(docs[key].get("evidence", []))
        evidence.extend(docs[key].get("outputs", []))
        evidence.extend(docs[key].get("helpers", []))
        evidence.append(docs[key]["source"])
    evidence.extend([*docs["reader"]["checked_inputs"], docs["reader"]["pdf_preview"],
        docs["reader"]["source"], *docs["failure"]["outputs"], docs["failure"]["source"]])
    for ref in evidence: check(ref)
    with output.open("x") as stream: stream.write(text)
    with claims_path.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(claim_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(claim_rows)
    for ref in evidence: check(ref)
    _, effects = validate(docs, refs)
    result = dict(schema="composed_native_main_generation_v1", source=record(__file__), inputs=refs,
        evidence=evidence, output=record(output), claims=record(claims_path), native_cells=list(figure.CELLS),
        admitted_endpoint_count=30, unavailable_score_cells=2, conditional_interval_rows=3*len(effects),
        adjusted_f1_interpretation=interpretation(effects), historical_source_modified=False,
        new_scoring_or_admission=False, new_bootstrap_draws=0, consumed_native_bootstrap_draws=100000,
        independent_confirmation=False, publication_ready=False, render_review_established=False)
    with receipt.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("evidence-directory", "output", "receipt", "claims-output"):
        parser.add_argument("--"+name, type=Path, required=True)
    for name in ("parent", "snapshot", "intervals", "figure", "reader", "failure"):
        parser.add_argument("--"+name, nargs=2, metavar=("PATH", "SHA256"), required=True)
    args = parser.parse_args()
    specifications = {name: getattr(args, name) for name in ("parent", "snapshot", "intervals", "figure", "reader", "failure")}
    result = generate(args.evidence_directory, specifications, args.output, args.receipt, args.claims_output)
    print(json.dumps(dict(admitted_endpoints=30, conditional_interval_rows=result["conditional_interval_rows"])))


if __name__ == "__main__":
    main()
