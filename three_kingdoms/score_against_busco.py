#!/usr/bin/env python3
"""Score a tool's predicted orthogroups against the BUSCO reference set.

Metric (BUSCO reference-gene pairs, micro-averaged):

For each reference orthogroup R (a set of N reference genes), we map
each ref-gene to the predicted OG it landed in. The prediction is
*correct* on a gene-pair (g_i, g_j) in R when both genes are in the
same predicted OG. Aggregating across all gene-pairs gives micro
precision/recall:

    TP = pairs of ref-genes co-predicted (in same predicted OG)
    FP = pairs of predicted-genes (where each is also in some ref-OG)
         that are co-predicted but NOT co-referenced
    FN = pairs of ref-genes NOT co-predicted

We restrict the universe to genes that appear in *some* reference OG —
that's the only set where ground-truth co-membership is defined.

Outputs precision, recall, F-score, and per-tool diagnostics.

Usage:
    score_against_busco.py --predictions <orthogroups.txt> \\
                           --reference  busco/reference_orthogroups.txt
"""
import argparse
from collections import defaultdict
from itertools import combinations
from pathlib import Path


def read_ogs(path):
    """Return list of frozensets of gene IDs."""
    out = []
    observed = {}
    with open(path) as f:
        for line_number, line in enumerate(f, start=1):
            parts = line.split()
            if not parts:
                continue
            # Strip trailing colons used by some tools (OG0:, OG1:, ...)
            if parts[0].endswith(":"):
                parts = parts[1:]
            if parts:
                group = frozenset(parts)
                if len(group) != len(parts):
                    raise ValueError(f"Duplicate gene within {path}:{line_number}")
                for gene in group:
                    if gene in observed:
                        raise ValueError(
                            f"Gene {gene!r} occurs in multiple groups in {path}: "
                            f"lines {observed[gene]} and {line_number}"
                        )
                    observed[gene] = line_number
                out.append(group)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--predictions", required=True, type=Path)
    ap.add_argument("--reference",   required=True, type=Path)
    ap.add_argument("--label",       default="tool")
    args = ap.parse_args()

    ref_ogs  = read_ogs(args.reference)
    pred_ogs = read_ogs(args.predictions)

    # Build gene -> predicted-OG-index map, and the universe of ref-genes
    pred_by_gene = {}
    for i, og in enumerate(pred_ogs):
        for g in og:
            pred_by_gene[g] = i

    ref_genes = set()
    for og in ref_ogs:
        ref_genes |= og

    n_ref_in_pred = sum(1 for g in ref_genes if g in pred_by_gene)
    n_ref_total   = len(ref_genes)

    # Gene-pair TP / FN over reference: each ref-OG generates C(|R|,2) gene-pairs;
    # TP if both genes share a predicted OG.
    tp = fn = 0
    for og in ref_ogs:
        members = list(og)
        for a, b in combinations(members, 2):
            pa = pred_by_gene.get(a); pb = pred_by_gene.get(b)
            if pa is not None and pb is not None and pa == pb:
                tp += 1
            else:
                fn += 1

    # FP: pairs of REF genes that share a predicted OG but are in
    # *different* ref OGs. We only count over genes that appear in some
    # ref OG, since that's the universe where "should they be together"
    # is well-defined (uses union over ref OGs).
    ref_gene_to_ref_og = {}
    for i, og in enumerate(ref_ogs):
        for g in og:
            ref_gene_to_ref_og[g] = i

    pred_to_ref_genes = defaultdict(list)
    for g, p in pred_by_gene.items():
        if g in ref_gene_to_ref_og:
            pred_to_ref_genes[p].append(g)

    fp = 0
    for p, genes in pred_to_ref_genes.items():
        for a, b in combinations(genes, 2):
            if ref_gene_to_ref_og[a] != ref_gene_to_ref_og[b]:
                fp += 1

    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall    = tp / (tp + fn) if (tp + fn) else 0.0
    f         = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0

    print(f"=== {args.label} ===")
    print(f"  reference OGs            : {len(ref_ogs):>10,}")
    print(f"  reference genes          : {n_ref_total:>10,}")
    print(f"  ref-genes in prediction  : {n_ref_in_pred:>10,} ({100*n_ref_in_pred/max(1,n_ref_total):.1f}%)")
    print(f"  predicted OGs            : {len(pred_ogs):>10,}")
    print(f"  TP gene pairs            : {tp:>10,}")
    print(f"  FP gene pairs            : {fp:>10,}")
    print(f"  FN gene pairs            : {fn:>10,}")
    print(f"  precision                : {precision:>10.4f}")
    print(f"  recall                   : {recall:>10.4f}")
    print(f"  F-score                  : {f:>10.4f}")

    return precision, recall, f


if __name__ == "__main__":
    main()
