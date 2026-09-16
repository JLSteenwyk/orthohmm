"""Frozen YGOB co-membership statistic; not resolved pairwise orthology."""

from collections import Counter
import csv

import numpy as np


METRICS = ("f1", "precision", "recall")


def read_predictions(path, format):
    """Read explicit native group formats without silently deduplicating genes."""
    groups = {}
    with path.open() as handle:
        if format == "root_hogs":
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != ["root_hog", "source_family", "genes"]:
                raise ValueError("Unexpected root-HOG header")
            rows = []
            for row in reader:
                if None in row or any(value is None for value in row.values()):
                    raise ValueError("Malformed root-HOG row")
                rows.append((row["root_hog"], row["genes"].split(",")))
        elif format == "named_groups":
            rows = []
            for line in handle:
                if not line.strip():
                    continue
                name, separator, genes = line.partition(":")
                if not separator:
                    raise ValueError("Expected group name followed by colon")
                rows.append((name.strip(), genes.split()))
        else:
            raise ValueError("Specify root_hogs or named_groups")
    for name, genes in rows:
        if not name or name in groups:
            raise ValueError("Empty or duplicate group name")
        groups[name] = genes
    return groups


def membership(groups):
    index = {}
    for name, genes in groups.items():
        if not isinstance(name, str) or not name or isinstance(genes, str) or not genes:
            raise ValueError("Groups need names and nonempty gene collections")
        for gene in genes:
            if not isinstance(gene, str) or not gene or gene.strip() != gene:
                raise ValueError("Invalid gene identifier")
            if gene in index:
                raise ValueError(f"Duplicate gene membership: {gene}")
            index[gene] = name
    return index


def statistics(counts):
    """Return F1, precision, recall as fractions; undefined ratios are zero."""
    tp, fp, fn = np.moveaxis(np.asarray(counts, dtype=float), -1, 0)
    def ratio(n, d):
        return np.divide(n, d, out=np.zeros_like(n), where=d > 0)
    return np.stack((ratio(2 * tp, 2 * tp + fp + fn),
                     ratio(tp, tp + fp), ratio(tp, tp + fn)), axis=-1)


def score_groups(predictions, references, input_genes):
    reference_index = membership(references)
    predicted_index = membership(predictions)
    if not reference_index:
        raise ValueError("Empty reference")
    input_genes = set(input_genes)
    if not set(reference_index) <= input_genes:
        raise ValueError("Reference genes absent from inference inputs")
    if not set(predicted_index) <= input_genes:
        raise ValueError("Unknown prediction IDs: input mapping must be verified")
    records = {name: {"pillar": name, "genes": len(genes), "tp": 0, "fp": 0.0,
                      "fn": len(genes) * (len(genes) - 1) // 2,
                      "covered_genes": 0, "exact": False}
               for name, genes in sorted(references.items())}
    projected_groups = 0
    predicted_pairs = 0
    for genes in predictions.values():
        counts = Counter(reference_index[g] for g in genes if g in reference_index)
        size = sum(counts.values())
        projected_groups += size > 0
        predicted_pairs += size * (size - 1) // 2
        for name, count in counts.items():
            record = records[name]
            true_pairs = count * (count - 1) // 2
            record["tp"] += true_pairs
            record["fn"] -= true_pairs
            # Each cross-pillar pair contributes one half to each endpoint pillar.
            record["fp"] += count * (size - count) / 2
            record["covered_genes"] += count
            record["exact"] |= count == size == record["genes"]
    totals = [sum(r[key] for r in records.values()) for key in ("tp", "fp", "fn")]
    if totals[0] + totals[1] != predicted_pairs:
        raise AssertionError("Allocated counts do not recover predicted pairs")
    return {
        "counts": dict(zip(("tp", "fp", "fn"), totals)),
        "metrics": dict(zip(METRICS, statistics(totals).tolist())),
        "reference_genes": len(reference_index),
        "covered_reference_genes": len(set(predicted_index) & set(reference_index)),
        "reference_gene_coverage": len(set(predicted_index) & set(reference_index)) / len(reference_index),
        "reference_groups": len(records),
        "represented_reference_groups": sum(r["covered_genes"] > 0 for r in records.values()),
        "predicted_groups_before_projection": len(predictions),
        "predicted_groups_with_scored_genes": projected_groups,
        "predicted_group_coverage": projected_groups / len(predictions) if predictions else 0.0,
        "predicted_group_coverage_definition": "fraction of supplied groups retaining at least one scored gene",
        "excluded_predicted_genes": len(set(predicted_index) - set(reference_index)),
        "exact_reference_groups": sum(r["exact"] for r in records.values()),
        "records": list(records.values()),
    }


def paired_bootstrap(scores, baseline="orthofinder_full", replicates=20000,
                     seed=20260917, batch_size=128):
    """Paired pillar resampling with bounded memory, recomputing micro ratios."""
    if baseline not in scores or len(scores) < 2:
        raise ValueError("Need baseline and comparators")
    if replicates < 100 or batch_size < 1:
        raise ValueError("Invalid bootstrap size")
    arrays = {}
    expected = None
    for method, score in sorted(scores.items()):
        records = sorted(score["records"], key=lambda r: r["pillar"])
        signature = [(r["pillar"], r["genes"]) for r in records]
        if not records or len({r["pillar"] for r in records}) != len(records):
            raise ValueError("Nonempty unique pillars required")
        if expected is not None and signature != expected:
            raise ValueError("Methods must share reference pillars and sizes")
        expected = signature
        counts = np.array([[r[k] for k in ("tp", "fp", "fn")] for r in records], dtype=float)
        sizes = np.array([r["genes"] for r in records], dtype=float)
        if (not np.isfinite(counts).all() or (counts < 0).any()
                or not np.isfinite(sizes).all() or (sizes < 1).any()
                or (sizes != np.floor(sizes)).any()
                or not np.array_equal(counts[:, 0] + counts[:, 2], sizes * (sizes - 1) / 2)):
            raise ValueError("Invalid pillar sufficient statistics")
        arrays[method] = counts
    rng = np.random.Generator(np.random.PCG64(seed))
    draws = {m: np.empty((replicates, 3)) for m in arrays}
    n = len(expected)
    for start in range(0, replicates, batch_size):
        end = min(start + batch_size, replicates)
        weights = rng.multinomial(n, np.full(n, 1 / n), size=end - start)
        for method, counts in arrays.items():
            draws[method][start:end] = statistics(weights @ counts)
    observed = {m: statistics(a.sum(axis=0)) for m, a in arrays.items()}
    comparisons = {}
    contrasts = (len(arrays) - 1) * len(METRICS)
    for method in arrays:
        if method == baseline:
            continue
        delta = 100 * (draws[method] - draws[baseline])
        comparisons[method] = {metric: {
            "difference_percentage_points": float(100 * (observed[method][i] - observed[baseline][i])),
            "paired_95_percent_ci": np.quantile(delta[:, i], [0.025, 0.975]).tolist(),
            "bonferroni_ci": np.quantile(delta[:, i], [0.025 / contrasts, 1 - 0.025 / contrasts]).tolist(),
        } for i, metric in enumerate(METRICS)}
    return {"baseline": baseline, "replicates": replicates, "seed": seed,
            "rng": "PCG64 multinomial", "numpy_version": np.__version__,
            "multiplicity_count": contrasts, "comparisons": comparisons,
            "limitations": ["Approximate intervals assume exchangeable reference pillars.",
                            "Cross-pillar false positives are allocated half to each endpoint pillar.",
                            "Group recovery does not establish resolved pairwise orthology."]}
