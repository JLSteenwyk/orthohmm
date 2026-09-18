"""Chunked, label-free diagnostics for canonically sorted numeric search hits."""

import numpy as np


SCORE_EDGES = np.array([0., .01, .03, .1, .3, 1., 3., 10., 30., 100., np.inf])


def blocks(queries, targets, scores, n, chunk_size=100000):
    if type(chunk_size) is not int or chunk_size < 1 or type(n) is not int or n < 1:
        raise ValueError("Positive integer universe and chunk size required")
    if n > int(np.sqrt(np.iinfo(np.int64).max)):
        raise ValueError("Gene universe cannot be encoded in int64")
    if (any(a.ndim != 1 for a in (queries, targets, scores))
            or not (len(queries) == len(targets) == len(scores))
            or queries.dtype.kind not in "iu" or targets.dtype.kind not in "iu"
            or scores.dtype.kind not in "fiu"):
        raise ValueError("Invalid hit arrays")
    previous = -1
    for start in range(0, len(scores), chunk_size):
        q, t, s = (a[start:start + chunk_size] for a in (queries, targets, scores))
        if (np.any(q < 0) or np.any(t < 0) or np.any(q >= n) or np.any(t >= n)
                or not np.isfinite(s).all() or np.any(s <= 0)):
            raise ValueError("Invalid hit indices or scores")
        codes = q.astype(np.int64) * n + t.astype(np.int64)
        if codes[0] <= previous or np.any(codes[1:] <= codes[:-1]):
            raise ValueError("Hits must have unique, increasing query-target keys")
        previous = int(codes[-1])
        yield codes, q, t, s


def summarize(queries, targets, scores, species, chunk_size=100000):
    species = np.asarray(species)
    if species.ndim != 1 or species.dtype.kind not in "iu" or np.any(species < 0):
        raise ValueError("Invalid species ownership")
    labels, owners = np.unique(species, return_inverse=True)
    n, k = len(species), len(labels)
    outgoing = np.zeros(n, dtype=np.int64)
    target_seen = np.zeros(n, dtype=bool)
    nonself_seen = np.zeros(n, dtype=bool)
    cross_seen = np.zeros(n, dtype=bool)
    # One bit-valued byte per gene/target species, independent of hit count.
    directions_seen = np.zeros((n, k), dtype=bool)
    directions = np.zeros(k * k, dtype=np.int64)
    histogram = np.zeros(len(SCORE_EDGES) - 1, dtype=np.int64)
    hits = self_hits = cross_hits = 0
    score_min = score_max = None
    for _, q, t, s in blocks(queries, targets, scores, n, chunk_size):
        hits += len(s)
        self_hits += int(np.count_nonzero(q == t))
        qs, ts = owners[q], owners[t]
        cross = qs != ts
        cross_hits += int(np.count_nonzero(cross))
        np.add.at(outgoing, q, 1)
        target_seen[t] = True
        nonself_seen[q[q != t]] = True
        cross_seen[q[cross]] = True
        directions_seen[q, ts] = True
        directions += np.bincount(qs * k + ts, minlength=k * k)
        histogram += np.histogram(s, bins=SCORE_EDGES)[0]
        low, high = float(s.min()), float(s.max())
        score_min = low if score_min is None else min(score_min, low)
        score_max = high if score_max is None else max(score_max, high)
    probabilities = [0, .25, .5, .75, .95, .99, 1]
    rows = []
    for source in range(k):
        mask = owners == source
        seen = directions_seen[mask].sum(axis=0)
        gene_count = int(mask.sum())
        for target in range(k):
            rows.append({"query_species": int(labels[source]), "target_species": int(labels[target]),
                "directed_hits": int(directions[source * k + target]),
                "queries_with_hits": int(seen[target]), "query_species_genes": gene_count})
    return {"genes": n, "directed_hits": hits, "self_hits": self_hits,
        "nonself_hits": hits - self_hits, "cross_species_hits": cross_hits,
        "queries_without_hits": int(np.count_nonzero(outgoing == 0)),
        "queries_without_nonself_hits": int(n - nonself_seen.sum()),
        "queries_without_cross_species_hits": int(n - cross_seen.sum()),
        "targets_without_hits": int(n - target_seen.sum()), "species_directions": rows,
        "quantile_probabilities": probabilities,
        "outgoing_hit_count_quantiles": np.quantile(outgoing, probabilities).tolist(),
        "normalized_score_min": score_min, "normalized_score_max": score_max,
        "normalized_score_histogram": {"edges": [float(x) for x in SCORE_EDGES[:-1]] + [None],
            "last_edge": "positive_infinity", "intervals": "left_closed_right_open",
            "counts": histogram.tolist()},
        "limitations": ["No biological reference labels or accuracy estimates are used.",
            "Score bins are descriptive, not calibrated across search engines.",
            "Reciprocity and exact score quantiles are not computed by this helper."]}


def overlap(first, second, n, chunk_size=100000):
    """Exact directed overlap using at most two encoded hit chunks at a time."""
    streams = [iter(blocks(*arrays, n, chunk_size)) for arrays in (first, second)]
    counts = [[0, 0], [0, 0]]

    def advance(index):
        item = next(streams[index], None)
        if item is None:
            return None
        codes, q, t, _ = item
        counts[index][0] += len(codes)
        counts[index][1] += int(np.count_nonzero(q != t))
        return codes

    a, b = advance(0), advance(1)
    common = [0, 0]
    while a is not None and b is not None:
        boundary = min(a[-1], b[-1])
        ai, bi = np.searchsorted(a, boundary, side="right"), np.searchsorted(b, boundary, side="right")
        shared = np.intersect1d(a[:ai], b[:bi], assume_unique=True)
        common[0] += len(shared)
        common[1] += int(np.count_nonzero(shared // n != shared % n))
        a = advance(0) if ai == len(a) else a[ai:]
        b = advance(1) if bi == len(b) else b[bi:]
    # Exhaust both streams even when one is empty: validate tails and count uniques.
    for index, current in enumerate((a, b)):
        while current is not None:
            current = advance(index)
    result = {}
    for column, label in enumerate(("all", "nonself")):
        left, right = counts[0][column], counts[1][column]
        shared = common[column]
        union = left + right - shared
        result[label] = {"intersection": shared, "first_only": left - shared,
            "second_only": right - shared, "union": union,
            "jaccard": shared / union if union else None,
            "fraction_of_first_recovered": shared / left if left else None,
            "fraction_of_second_recovered": shared / right if right else None}
    return result
