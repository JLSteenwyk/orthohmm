"""Render the finite-panel mechanism addendum from checked counts."""

import argparse
import hashlib
import json
from pathlib import Path


def render(report):
    if (report["status"] != "independent_count_and_candidate_decomposition_verified"
            or len(report["cells"]) != 70 or len(report["summary"]) != 7
            or any(r["scored_cells"] != 10 for r in report["summary"])):
        raise ValueError("Incomplete checked diagnostic")
    lines = ["# Fixed-Candidate Gene-Tree Oracle Results", "",
        "## Methods", "",
        "This development-exposed mechanism diagnostic retains all70 generating-",
        "species-tree simulation cells, seven conditions and ten seeds each. It",
        "holds candidates, species trees, frozen reconciliation rules and constraint",
        "policy fixed. Only actually inferred single-ancestral-family candidates",
        "receive generating gene trees, either with their generating root or with",
        "frozen minimum-duplication/loss rerooting. All bypasses/ineligible candidates",
        "retain native predictions. The [prospective protocol](SIMULATION_GENE_TREE_ORACLE_PROTOCOL_20261004.md)",
        "and tested runner were committed/pushed before new scoring.", "",
        "## Accuracy", "",
        "F1 values are percentages, arithmetic means of ten dataset scores; changes",
        "are percentage points. These are descriptive finite-panel effects, not",
        "population confidence intervals, significance tests or independent validation.", "",
        "| Condition | Inferred | Generating root | Generating, rerooted | Generating root - inferred |",
        "| --- | ---: | ---: | ---: | ---: |"]
    for row in report["summary"]:
        means = row["mean_metrics"]
        values = [100 * means[a]["f1"] for a in ("inferred", "generating_root", "generating_rerooted")]
        lines.append(f"| {row['condition']} | {values[0]:.3f} | {values[1]:.3f} | {values[2]:.3f} | {values[1]-values[0]:+.3f} |")
    candidate_count = sum(sum(r["candidate_counts"].values()) for r in report["summary"])
    topo = report["topology"]
    lines.extend(["", f"All70 baseline pair sets reproduce native output; independent readback also",
        f"matches every original retained generating-species-tree score. There are",
        f"{candidate_count:,} candidates and {topo['eligible']:,} eligible gene-tree controls.",
        f"Among eligible controls, {topo['unrooted_disagreement']} have unrooted topology disagreement;",
        f"{topo['root_only_disagreement']} more differ only in root ({topo['rooted_disagreement']} rooted disagreements total).",
        "This does not establish that every disagreement causes an accuracy error.", "",
        "## Residual Errors", "",
        "Counts below sum ten seeds per condition; they are not the denominator of",
        "the macro-mean F1 above. The cross-candidate count is independently derived",
        "from truth and the native candidate partition, not inferred from a score.", "",
        "| Condition | Generating-root FN | True pairs across candidates | Across-candidate share | Eligible within-candidate FN |",
        "| --- | ---: | ---: | ---: | ---: |"])
    for row in report["summary"]:
        fn, across = row["oracle_total_fn"], row["true_pairs_across_candidates"]
        share = f"{100 * across / fn:.3f}%" if fn else "undefined"
        within = row["residual_by_arm"]["generating_root"]["oracle_eligible"]["fn"]
        lines.append(f"| {row['condition']} | {fn:,} | {across:,} | {share} | {within:,} |")
    lines.extend(["", "In both divergent conditions more than99.7% of residual generating-root",
        "false negatives are fixed upstream by different candidate membership.",
        "Replacing gene trees alone cannot recover those pairs. This localizes this",
        "deficit upstream of gene-tree reconciliation, but does not distinguish HMM",
        "search, grouping and candidate expansion. Generating trees are not uniformly",
        "perfect: retain the within-candidate errors, bypass false positives and",
        "losses as well as gains. No defaults or benchmark scores are changed.", "",
        "## Limits And Reproduction", "",
        "True gene trees/roots are unavailable in ordinary inference. The oracle",
        "also changes lengths/support and is not topology-only causal evidence or an",
        "accuracy upper bound. Native pair scoring does not validate ancestral-copy",
        "root-HOG membership. This narrow mechanism result does not close all error-",
        "analysis, independent-generalization or uncertainty requirements.", "",
        "The [compact readback](simulation_gene_tree_oracle_readback_20261004.json)",
        "retains all70 score/decomposition rows; the detailed local report records",
        "every candidate and all source/input identities. Readback independently",
        "checks counts/partitions and original scores, not the new oracle pair sets.",
        f"It rehashes {report['input_identities_rechecked']:,} input/source identities.", "",
        f"Detailed local report: `{report['detailed_report']['path']}`.",
        f"SHA256: `{report['detailed_report']['sha256']}`.", "",
        "The [execution receipt](simulation_gene_tree_oracle_execution_20261004.json)",
        "records successful execution, tests and the failed first readback. The first",
        "readback9e450025 rejected before output because it expected `scored` rather",
        "than the retained `complete` status. Corrected9893f6a0 explicitly validates",
        "all70 baseline completions; no seed or outcome was dropped. The44 diagnostic/",
        "readback tests and two renderer tests pass. No completed inference or timing",
        "measurement was repeated.", "",
        "This is a small read-only diagnostic on the shared Threadripper, not a",
        "matched timing run. No unrelated analysis was stopped or modified. Existing",
        "rc2/PDF artifacts remain unchanged; this is a later scientific addendum.", ""])
    return "\n".join(lines)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    raw = args.report.read_bytes()
    if hashlib.sha256(raw).hexdigest() != args.report_sha256:
        raise ValueError("Changed checked report")
    result = render(json.loads(raw))
    with args.output.open("x") as handle:
        handle.write(result)
