"""Frozen label-independent simulation transformations and pair scoring."""

import hashlib
import itertools


def ranked_names(names, seed, condition):
    if not isinstance(seed, int) or seed <= 0:
        raise ValueError("Positive integer seed required")
    if len(set(names)) != len(names) or any(not isinstance(n, str) or not n for n in names):
        raise ValueError("Nonempty unique names required")
    return sorted(names, key=lambda n: (hashlib.sha256(f"{seed}:{condition}:{n}".encode("utf-8")).digest(), n))


def select_removals(gene_species, tree, seed, condition):
    """Selection has no access to reference labels, pairs, or sequences."""
    if not isinstance(seed, int) or seed <= 0:
        raise ValueError("Positive integer seed required")
    tips = [n.name for n in tree.get_terminals()]
    if not tips or len(set(tips)) != len(tips) or set(gene_species.values()) - set(tips):
        raise ValueError("Invalid species tree or unmatched input species")
    removed_species = []
    if condition == "missing20":
        removed_genes = ranked_names(list(gene_species), seed, condition)[:len(gene_species) // 5]
    elif condition in {"uneven_taxa", "taxon_count_control"}:
        candidates = []
        for clade in tree.get_nonterminals():
            if clade is tree.root:
                continue
            names = tuple(sorted(n.name for n in clade.get_terminals()))
            if 2 <= len(names) <= 4:
                candidates.append(names)
        if not candidates:
            return {"status": "inapplicable", "reason": "No eligible non-root clade", "condition": condition}
        chosen = min(candidates, key=lambda names: (-len(names), names))
        removed_species = list(chosen[1:])
        if condition == "taxon_count_control":
            removed_species = ranked_names(tips, seed, condition)[:len(removed_species)]
        removed_genes = [g for g, s in gene_species.items() if s in set(removed_species)]
    else:
        raise ValueError("Unknown frozen transformation")
    return {"status": "selected", "condition": condition, "seed": seed,
            "removed_genes": sorted(removed_genes), "removed_species": sorted(removed_species),
            "retained_species": sorted(set(tips) - set(removed_species))}


def canonical_pairs(pairs, gene_species):
    result = set()
    supplied = 0
    for pair in pairs:
        if isinstance(pair, str) or len(pair) != 2:
            raise ValueError("Pair must contain exactly two IDs")
        a, b = pair
        if a not in gene_species or b not in gene_species:
            raise ValueError("Unknown pair endpoint")
        if a == b or gene_species[a] == gene_species[b]:
            raise ValueError("Expected cross-species distinct gene pairs")
        supplied += 1
        result.add(tuple(sorted((a, b))))
    return result, supplied - len(result)


def project_truth(truth_pairs, gene_species, removals):
    if removals["status"] != "selected":
        raise ValueError("Cannot project an inapplicable transformation")
    pairs, duplicates = canonical_pairs(truth_pairs, gene_species)
    if duplicates:
        raise ValueError("Duplicate truth pairs")
    removed = set(removals["removed_genes"])
    if not removed <= gene_species.keys():
        raise ValueError("Unknown removed gene")
    retained = {g: s for g, s in gene_species.items() if g not in removed}
    projected = {p for p in pairs if set(p) <= retained.keys()}
    return retained, sorted(projected)


def group_pairs(groups, gene_species):
    seen = set()
    for group in groups:
        if isinstance(group, str) or not group:
            raise ValueError("Nonempty gene collection required")
        for gene in group:
            if gene not in gene_species or gene in seen:
                raise ValueError("Unknown or duplicate group membership")
            seen.add(gene)
        for a, b in itertools.combinations(group, 2):
            if gene_species[a] != gene_species[b]:
                yield tuple(sorted((a, b)))


def score_pairs(predictions, truth, gene_species):
    predicted, duplicates = canonical_pairs(predictions, gene_species)
    reference, truth_duplicates = canonical_pairs(truth, gene_species)
    if truth_duplicates:
        raise ValueError("Duplicate truth pairs")
    tp = len(predicted & reference)
    fp = len(predicted - reference)
    fn = len(reference - predicted)
    denominator = 2 * tp + fp + fn
    covered = {gene for pair in predicted for gene in pair}
    return {"tp": tp, "fp": fp, "fn": fn,
            "f1": 2 * tp / denominator if denominator else 0.0,
            "precision": tp / len(predicted) if predicted else 0.0,
            "recall": tp / len(reference) if reference else 0.0,
            "eligible_true_pairs": len(reference), "predicted_pairs": len(predicted),
            "duplicate_prediction_rows": duplicates, "input_genes": len(gene_species),
            "genes_in_predicted_pairs": len(covered),
            "pair_endpoint_coverage": len(covered) / len(gene_species) if gene_species else 0.0,
            "coverage_definition": "input genes occurring in at least one predicted cross-species pair",
            "undefined_ratios": [m for m, d in (("f1", denominator), ("precision", len(predicted)), ("recall", len(reference))) if not d]}
