"""Experimental candidate-only ordering policy; not a production default."""

import argparse
import itertools
from pathlib import Path

import numpy as np

POLICY = "canonical_directed_pair_v1"


def canonical_hit_order_v1(names, queries, targets, scores):
    """Preserve unique hit values, ordering by fixed lexical gene indices."""
    names = list(names)
    if (not names or any(not isinstance(n, str) or not n for n in names)
            or names != sorted(set(names))):
        raise ValueError("Require unique lexically sorted gene names")
    q, t, s = map(np.asarray, (queries, targets, scores))
    if (any(a.ndim != 1 for a in (q,t,s)) or not len(q) == len(t) == len(s)
            or q.dtype.kind not in "iu" or t.dtype.kind not in "iu" or s.dtype.kind not in "fiu"):
        raise ValueError("Invalid hit shape or dtype")
    if len(q) and (min(q.min(),t.min()) < 0 or max(q.max(),t.max()) >= len(names)):
        raise ValueError("Hit index outside gene universe")
    if not np.isfinite(s).all() or np.any(s <= 0):
        raise ValueError("Require finite positive normalized scores")
    order = np.lexsort((t,q))
    q, t, s = q[order], t[order], s[order]
    if len(q) > 1 and np.any((q[1:] == q[:-1]) & (t[1:] == t[:-1])):
        raise ValueError("Duplicate directed hits require an explicit upstream policy")
    return q, t, s


def fixture_inputs():
    rows = [(s,t,w) for s in (3,4) for t,w in zip((0,1,2),(.1,.2,.3))]
    rows += [(t,s,w) for s in (3,4) for t,w in zip((0,1,2),(.1,.2,.3))]
    return (np.array([r[0] for r in rows],dtype=np.int32),
            np.array([r[1] for r in rows],dtype=np.int32),
            np.array([r[2] for r in rows],dtype=np.float64))


def fixture_orders():
    for a in itertools.permutations(range(3)):
        for b in itertools.permutations(range(3)):
            yield np.array([*a, *(3+i for i in b), *(6+i for i in a), *(9+i for i in b)])


def fixture_run():
    from orthohmm import refinement
    from benchmark_tools.prepare_ob_candidate_neighborhood import record
    source = record(refinement.__file__)
    if source["sha256"] != "991f1eb6a5f73d0442529ed19095a34b7c6ba8bff8dfe43a1127e24ec73fb26d":
        raise ValueError("Fixture requires the retained frozen candidate implementation")
    inputs = fixture_inputs()
    names, clusters, species = list("abcde"), [[0,1,2],[3],[4]], np.arange(5)
    rows = []
    for order in fixture_orders():
        arrays = tuple(a[order] for a in inputs)
        groups = {}
        for policy, values in (("unchanged",arrays), (POLICY,canonical_hit_order_v1(names,*arrays))):
            result = refinement.merge_supported_satellite_candidate_clusters(
                clusters,*values,species,max_satellites_per_anchor=1,max_iterations=1)
            groups[policy] = sorted(sorted(g) for g in result[0])
        rows.append(dict(input_order=order.tolist(),groups=groups))
    return dict(status="synthetic_candidate_order_regression_exercised", source=record(__file__),
        scientific_source=source, numpy=np.__version__, policy=POLICY, genes=5, directed_hits=12,
        rows=rows, sequential_total=float(np.add.reduceat(np.array([.1,.2,.3]),[0])[0]),
        reversed_total=float(np.add.reduceat(np.array([.3,.2,.1]),[0])[0]),
        production_changed=False, accuracy_evaluated=False,
        limitations=["Synthetic floating-point/order mechanism fixture, not biological accuracy evidence.",
            "Canonicalization is experimental and not connected to production inference.",
            "Fixed lexical gene indexing, unique directed hits, fixed score values and runtime.",
            "Does not address dependency versions, cross-platform arithmetic, seed ordering or score perturbations.",
            "Full-data candidate and end-to-end validation remain required before any promotion."])


if __name__ == "__main__":
    from benchmark_tools.probe_installed_ob_clustering import write_json
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fixture-output", type=Path, required=True)
    a = p.parse_args()
    if a.fixture_output.exists():
        raise FileExistsError(a.fixture_output)
    write_json(a.fixture_output, fixture_run())
