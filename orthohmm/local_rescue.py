"""Label-independent, bounded nomination for experimental local tree rescue.

Nomination does not merge families or infer orthology. A nominee still
requires cross-profile validation and a joint-tree admission decision.
"""

from dataclasses import dataclass
import math
from numbers import Integral

import numpy as np


@dataclass(frozen=True)
class _DirectionalSupport:
    score: float
    queries: tuple[int, ...]
    targets: tuple[int, ...]


def _positive_integer(value, name, minimum=1):
    if (isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral)
            or value < minimum):
        raise ValueError(name + " must be a bounded positive integer")
    return int(value)


def _two_witnesses(old, first, second, sentinel):
    return tuple(sorted({*old, *[int(value) for value in (first, second) if value < sentinel]})[:2])


def _partition(clusters, names, species, max_joint_genes):
    count = len(names)
    if count > np.iinfo(np.int32).max or any(
        not isinstance(name, str) or not name for name in names
    ) or len(set(names)) != count:
        raise ValueError("Require unique nonempty global gene names within int32 bounds")
    if species.ndim != 1 or len(species) != count or species.dtype.kind not in "iu" or np.any(species < 0):
        raise ValueError("Require one nonnegative integer species label per gene")
    owner = np.full(count, -1, dtype=np.int32)
    eligible = np.zeros(len(clusters), dtype=bool)
    family_names, excluded = {}, []
    for family, cluster in enumerate(clusters):
        members = np.asarray(cluster)
        if (members.ndim != 1 or not len(members) or members.dtype.kind not in "iu"
                or np.any(members < 0) or np.any(members >= count)
                or len(np.unique(members)) != len(members)):
            raise ValueError("Candidate families must contain unique valid gene indices")
        if np.any(owner[members] >= 0):
            raise ValueError("Candidate families overlap")
        owner[members] = family
        species_count = len(np.unique(species[members]))
        reasons = []
        if len(members) < 3:
            reasons.append("too_few_genes")
        if len(members) > max_joint_genes - 3:
            reasons.append("family_exceeds_joint_budget")
        if species_count < 2:
            reasons.append("too_few_species")
        if reasons:
            excluded.append({"family_id": family, "genes": len(members), "species": species_count, "reasons": reasons})
        else:
            eligible[family] = True
            family_names[family] = tuple(sorted(names[int(index)] for index in members))
    if np.any(owner < 0):
        raise ValueError("Candidate partition does not cover the complete gene table")
    return owner, eligible, family_names, excluded


def _aggregate_chunk(stats, codes, queries, targets, scores, sentinel):
    if not len(codes):
        return
    order = np.argsort(codes, kind="stable")
    keys, q, t, values = codes[order], queries[order], targets[order], scores[order]
    starts = np.r_[0, np.flatnonzero(keys[1:] != keys[:-1]) + 1]
    sizes = np.diff(np.r_[starts, len(keys)])
    maxima = np.maximum.reduceat(values, starts)
    # Two distinct witnesses suffice for the gate; repeated hits add no support.
    first_q, first_t = np.minimum.reduceat(q, starts), np.minimum.reduceat(t, starts)
    second_q = np.minimum.reduceat(np.where(q != np.repeat(first_q, sizes), q, sentinel), starts)
    second_t = np.minimum.reduceat(np.where(t != np.repeat(first_t, sizes), t, sentinel), starts)
    for index, start in enumerate(starts):
        key = int(keys[start])
        previous = stats.get(key, _DirectionalSupport(0., (), ()))
        stats[key] = _DirectionalSupport(
            max(previous.score, float(maxima[index])),
            _two_witnesses(previous.queries, first_q[index], second_q[index], sentinel),
            _two_witnesses(previous.targets, first_t[index], second_t[index], sentinel),
        )


def nominate_local_contexts(
    clusters, gene_names, gene_to_species, hit_queries, hit_targets, hit_scores,
    *, max_contexts=128, max_joint_genes=200, chunk_size=500000,
):
    """Nominate disjoint unique mutual-best family pairs from existing HMM hits.

    Directional witnesses saturate at two distinct genes; they are not exact
    support counts or necessarily reciprocal individual gene-pair matches.
    No benchmark identity or reference annotations are accepted by this API.
    """
    max_contexts = _positive_integer(max_contexts, "Context budget")
    max_joint_genes = _positive_integer(max_joint_genes, "Joint gene budget", minimum=6)
    chunk_size = _positive_integer(chunk_size, "Hit chunk size")
    names = list(gene_names)
    species = np.asarray(gene_to_species)
    if not len(names) and species.shape == (0,):
        species = species.astype(np.int32)
    owner, eligible, family_names, excluded = _partition(clusters, names, species, max_joint_genes)
    queries, targets, scores = (np.asarray(array) for array in (hit_queries, hit_targets, hit_scores))
    if any(array.ndim != 1 for array in (queries, targets, scores)) or not len(queries) == len(targets) == len(scores):
        raise ValueError("HMM hit arrays must be aligned one-dimensional arrays")
    if not len(queries):
        queries, targets = queries.astype(np.int32), targets.astype(np.int32)
    if queries.dtype.kind not in "iu" or targets.dtype.kind not in "iu" or scores.dtype.kind not in "iuf":
        raise ValueError("HMM hit indices must be integers and scores real numbers")
    stats, used_hits = {}, 0
    families = len(clusters)
    for start in range(0, len(scores), chunk_size):
        q, t, s = queries[start:start + chunk_size], targets[start:start + chunk_size], scores[start:start + chunk_size]
        if (np.any(q < 0) or np.any(t < 0) or np.any(q >= len(names)) or np.any(t >= len(names))
                or not np.isfinite(s).all()):
            raise ValueError("Invalid HMM hit index or nonfinite score")
        left, right = owner[q], owner[t]
        keep = (left != right) & (species[q] != species[t]) & eligible[left] & eligible[right] & (s > 0)
        q, t, s, left, right = q[keep], t[keep], s[keep], left[keep], right[keep]
        used_hits += len(s)
        codes = left.astype(np.int64) * families + right
        _aggregate_chunk(stats, codes, q, t, s, len(names))

    unordered = sorted({tuple(sorted(divmod(code, families))) for code in stats})
    supported, pair_audit, neighbors = {}, [], {}
    for left, right in unordered:
        forward, reverse = stats.get(left * families + right), stats.get(right * families + left)
        if (forward is None or reverse is None or any(
            len(witnesses) < 2 for evidence in (forward, reverse) for witnesses in (evidence.queries, evidence.targets)
        )):
            pair_audit.append({"family_a": left, "family_b": right, "status": "insufficient_reciprocal_witnesses"})
            continue
        score = math.sqrt(forward.score) * math.sqrt(reverse.score)
        if not math.isfinite(score):
            raise ValueError("Nonfinite reciprocal family score")
        supported[(left, right)] = score
        neighbors.setdefault(left, []).append((score, right))
        neighbors.setdefault(right, []).append((score, left))
    ranking, best = {}, {}
    for family in sorted(family_names):
        ordered = sorted(neighbors.get(family, []), key=lambda item: (-item[0], item[1]))
        tied = [neighbor for score, neighbor in ordered if score == ordered[0][0]] if ordered else []
        selected = tied[0] if len(tied) == 1 else None
        best[family] = selected
        ranking[family] = {"family_id": family, "best_neighbor": selected,
            "best_score": ordered[0][0] if ordered else None, "best_tied_neighbors": sorted(tied),
            "runner_up_neighbor": ordered[1][1] if len(ordered) > 1 else None,
            "runner_up_score": ordered[1][0] if len(ordered) > 1 else None,
            "score_margin": ordered[0][0] - ordered[1][0] if len(ordered) > 1 else None}
    mutual = []
    for (left, right), score in sorted(supported.items()):
        row = {"family_a": left, "family_b": right, "score": score,
               "joint_genes": len(clusters[left]) + len(clusters[right])}
        if best[left] != right or best[right] != left:
            row["status"] = "not_unique_mutual_best"
        elif row["joint_genes"] > max_joint_genes:
            row["status"] = "joint_gene_budget"
        else:
            row["status"] = "eligible_before_context_budget"
            mutual.append(row)
        pair_audit.append(row)
    mutual.sort(key=lambda row: (-row["score"], *sorted((family_names[row["family_a"]], family_names[row["family_b"]]))))
    contexts = []
    for index, row in enumerate(mutual):
        if index >= max_contexts:
            row["status"] = "context_budget"
            continue
        row["status"] = "selected"
        left, right = row["family_a"], row["family_b"]
        contexts.append({**row, "context_id": f"LocalContext{index + 1:06d}",
            "genes_a": list(family_names[left]), "genes_b": list(family_names[right]),
            "a_to_b_score": stats[left * families + right].score,
            "b_to_a_score": stats[right * families + left].score,
            "family_a_ranking": ranking[left], "family_b_ranking": ranking[right]})
    selected_families = [family for row in contexts for family in (row["family_a"], row["family_b"])]
    if len(selected_families) != len(set(selected_families)):
        raise ValueError("Mutual-best nomination unexpectedly produced overlapping contexts")
    return {
        "schema_version": 1, "status": "nominated_not_profile_or_tree_validated",
        "parameters": {"max_contexts": max_contexts, "max_joint_genes": max_joint_genes,
                       "minimum_family_genes": 3, "minimum_family_species": 2,
                       "minimum_directional_witnesses": 2, "hit_chunk_size": chunk_size},
        "counts": {"genes": len(names), "families": families, "eligible_families": int(eligible.sum()),
                   "input_hits": len(scores), "eligible_positive_cross_hits": used_hits,
                   "directed_family_pairs": len(stats), "reciprocally_supported_pairs": len(supported),
                   "selected_contexts": len(contexts)},
        "excluded_families": excluded, "family_ranking": [ranking[key] for key in sorted(ranking)],
        "pair_audit": sorted(pair_audit, key=lambda row: (row["family_a"], row["family_b"])),
        "directional_support": [{"family_query": code // families, "family_target": code % families,
            "max_score": evidence.score, "query_witnesses": list(evidence.queries),
            "target_witnesses": list(evidence.targets)} for code, evidence in sorted(stats.items())],
        "contexts": contexts, "orthology_inferred": False, "families_merged": False,
    }


def assess_joint_outcome(outcome, genes_a, genes_b, baseline_high_pairs, gene_species):
    """Admit only a constrained reconciled context preserving baseline pairs.

    This gate does not write candidate groups or append native pairs. The
    final ordinary species-tree/reconciliation pass can still change its
    provisional relationships and must be evaluated separately.
    """
    genes_a, genes_b = tuple(genes_a), tuple(genes_b)
    if (not genes_a or not genes_b or any(not isinstance(gene, str) or not gene for gene in (*genes_a, *genes_b))
            or len(set(genes_a)) != len(genes_a) or len(set(genes_b)) != len(genes_b)
            or set(genes_a) & set(genes_b)):
        raise ValueError("Joint admission requires disjoint nonempty unique source families")
    left, right = set(genes_a), set(genes_b)
    universe = left | right
    if (len(outcome.genes) != len(universe) or set(outcome.genes) != universe
            or outcome.reconciliation is None):
        raise ValueError("Joint admission requires the complete reconciled two-family context")
    if any(gene not in gene_species or not isinstance(gene_species[gene], str) or not gene_species[gene]
           for gene in universe):
        raise ValueError("Joint context species ownership is missing or invalid")
    roots = {}
    for index, group in enumerate(outcome.root_groups):
        if not group:
            raise ValueError("Joint root groups must be nonempty")
        for gene in group:
            if gene not in universe or gene in roots:
                raise ValueError("Joint root groups have invalid or duplicate ownership")
            roots[gene] = index
    if set(roots) != universe:
        raise ValueError("Joint root groups do not cover the complete context")

    def pair_key(pair):
        if (len(pair) != 2 or any(not isinstance(gene, str) or gene not in universe for gene in pair) or pair[0] == pair[1]
                or gene_species[pair[0]] == gene_species[pair[1]]):
            raise ValueError("Invalid joint cross-species ortholog pair")
        return tuple(sorted(pair))

    pairs = {pair_key(pair) for pair in outcome.ortholog_pairs}
    if len(pairs) != len(outcome.ortholog_pairs) or any(roots[a] != roots[b] for a, b in pairs):
        raise ValueError("Joint pairs are duplicated or cross root-group boundaries")
    confidence = {}
    for row in outcome.ortholog_pair_confidence:
        if len(row) != 3 or row[2] not in ("high", "medium", "low"):
            raise ValueError("Invalid joint pair confidence")
        pair = pair_key(row[:2])
        if pair in confidence:
            raise ValueError("Duplicate joint pair confidence")
        confidence[pair] = row[2]
    if set(confidence) != pairs:
        raise ValueError("Joint pair/confidence universes differ")
    baseline_rows = tuple(baseline_high_pairs)
    baseline = {pair_key(pair) for pair in baseline_rows}
    if len(baseline) != len(baseline_rows):
        raise ValueError("Baseline high pairs are duplicated")
    if any((a in left) != (b in left) for a, b in baseline):
        raise ValueError("Baseline high pairs cross the original family boundaries")
    lost = baseline - pairs
    downgraded = {pair for pair in baseline & pairs if confidence[pair] != "high"}
    cross = {pair for pair in pairs if confidence[pair] == "high" and ((pair[0] in left) != (pair[1] in left))}
    supported_a = {gene for pair in cross for gene in pair if gene in left}
    supported_b = {gene for pair in cross for gene in pair if gene in right}
    species_a = {gene_species[gene] for gene in supported_a}
    species_b = {gene_species[gene] for gene in supported_b}
    reasons = []
    if lost:
        reasons.append("baseline_high_pair_loss")
    if min(len(supported_a), len(supported_b)) < 2:
        reasons.append("insufficient_novel_high_gene_support")
    if min(len(species_a), len(species_b)) < 2:
        reasons.append("insufficient_novel_high_species_support")
    reconciliation = outcome.reconciliation
    return {
        "accepted": not reasons, "reasons": reasons, "family_id": outcome.family_id,
        "genes": len(universe), "root_groups": len(outcome.root_groups),
        "baseline_high_pairs": len(baseline), "baseline_high_pairs_retained": len(baseline & pairs),
        "baseline_high_pairs_lost": [list(pair) for pair in sorted(lost)],
        "baseline_high_pairs_downgraded": [list(pair) for pair in sorted(downgraded)],
        "novel_cross_high_pairs": [list(pair) for pair in sorted(cross)],
        "supported_genes_a": sorted(supported_a), "supported_genes_b": sorted(supported_b),
        "supported_species_a": sorted(species_a), "supported_species_b": sorted(species_b),
        "duplications": reconciliation.duplication_count, "speciations": reconciliation.speciation_count,
        "uncertain_events": reconciliation.uncertain_count,
        "scope": "Provisional constrained joint-tree admission, not a final pair/group export or accuracy claim",
    }
