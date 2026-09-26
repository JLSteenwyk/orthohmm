"""Classify watched pairs from an unfiltered indexed search result.

This describes the observed search only. It cannot reconstruct rejected
candidates from a significant-hit cache or establish historical equivalence.
"""

import math


def watched_candidate_arrays(query_ids, target_ids, pairs):
    """Construct explicit candidates for a reference-conditioned diagnostic."""
    import numpy as np

    queries = {gene: i for i, gene in enumerate(query_ids)}
    targets = {gene: i for i, gene in enumerate(target_ids)}
    if len(queries) != len(query_ids) or len(targets) != len(target_ids):
        raise ValueError("Duplicate sequence identities")
    pairs = list(pairs)
    if len(set(pairs)) != len(pairs):
        raise ValueError("Duplicate watched pair")
    by_query = [[] for _ in query_ids]
    for q, t in pairs:
        if q not in queries or t not in targets:
            raise ValueError("Watched pair outside sequence universe")
        by_query[queries[q]].append(targets[t])
    offsets = np.zeros(len(query_ids) + 1, dtype=np.int64)
    offsets[1:] = np.cumsum([len(items) for items in by_query])
    candidates = np.array([i for items in by_query for i in sorted(items)], dtype=np.int32)
    return candidates, offsets


def classify_search_result(result, query_ids, target_ids, watched_pairs, threshold):
    """Return directed decisions, requiring complete unfiltered candidate output."""
    if not math.isfinite(threshold) or threshold <= 0:
        raise ValueError("Require a finite positive E-value threshold")
    if len(set(query_ids)) != len(query_ids) or len(set(target_ids)) != len(target_ids):
        raise ValueError("Duplicate sequence identities")
    sizes = [len(result.query_indices), len(result.target_indices),
             len(result.scores), len(result.evalues)]
    if len(set(sizes)) != 1 or sizes[0] != result.candidate_count:
        raise ValueError("Require complete unfiltered candidate results")
    observed = {}
    for q, t, score, evalue in zip(result.query_indices, result.target_indices,
                                  result.scores, result.evalues):
        if int(q) != q or int(t) != t or not 0 <= q < len(query_ids) or not 0 <= t < len(target_ids):
            raise ValueError("Candidate index outside sequence universe")
        pair = (query_ids[int(q)], target_ids[int(t)])
        if pair in observed:
            raise ValueError("Duplicate directed candidate")
        score, evalue = float(score), float(evalue)
        if not math.isfinite(score) or not math.isfinite(evalue) or evalue < 0:
            raise ValueError("Invalid numerical scoring result; do not label as rejection")
        observed[pair] = (score, evalue)
    query_set, target_set = set(query_ids), set(target_ids)
    watched = list(watched_pairs)
    if len(set(watched)) != len(watched):
        raise ValueError("Duplicate watched pair")
    rows = []
    for q, t in watched:
        if q not in query_set or t not in target_set:
            raise ValueError("Watched pair outside sequence universe")
        values = observed.get((q, t))
        if values is None:
            decision, score, evalue = "not_selected_by_prefilter", None, None
        else:
            score, evalue = values
            decision = "accepted" if evalue < threshold else "scored_not_significant"
        rows.append({"query": q, "target": t, "decision": decision,
                     "score": score, "evalue": evalue})
    return rows
