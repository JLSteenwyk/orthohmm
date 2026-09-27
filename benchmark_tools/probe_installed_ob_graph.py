"""Factor retained score/order differences at the initial graph boundary."""

import argparse
import itertools
import json
from pathlib import Path
import pickle
import sys

import numpy as np

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.trace_ob_initial_edges import capture_edges, CORE_SHA

SEARCH_SHA = "247133956257e25cfbf5b2b79a19544ea382de07f9a7393f275d707ee8a658c5"


def align_scores(n, old_q, old_t, old_s, new_q, new_t, new_s):
    """Align two unique nonself hit sets while preserving each input order."""
    keys = []
    for q, t, s in ((old_q, old_t, old_s), (new_q, new_t, new_s)):
        if (any(a.ndim != 1 for a in (q, t, s)) or not len(q)
                or not len(q) == len(t) == len(s) or q.dtype.kind not in "iu"
                or t.dtype.kind not in "iu" or min(q.min(), t.min()) < 0
                or max(q.max(), t.max()) >= n or np.any(q == t)
                or not np.isfinite(s).all() or np.any(s <= 0)):
            raise ValueError("Invalid directed nonself hits")
        keys.append(q.astype(np.int64) * n + t)
    orders = [np.argsort(k, kind="stable") for k in keys]
    sorted_keys = [k[o] for k, o in zip(keys, orders)]
    if (not np.array_equal(*sorted_keys) or np.any(np.diff(sorted_keys[0]) == 0)):
        raise ValueError("Hit sets must be identical and unique")
    old_in_new = np.empty_like(old_s)
    new_in_old = np.empty_like(new_s)
    old_in_new[orders[1]] = old_s[orders[0]]
    new_in_old[orders[0]] = new_s[orders[1]]
    return old_in_new, new_in_old


def edge_comparison(a, b, n):
    keys, weights = [], []
    for graph in (a, b):
        k = graph.sources.astype(np.int64) * n + graph.targets
        order = np.argsort(k)
        if np.any(graph.sources >= graph.targets) or np.any(np.diff(k[order]) == 0):
            raise ValueError("Edges must be canonical and unique")
        keys.append(k[order])
        weights.append(graph.weights[order])
    common, i, j = np.intersect1d(*keys, assume_unique=True, return_indices=True)
    delta = np.abs(weights[0][i] - weights[1][j])
    return dict(left_edges=len(keys[0]), right_edges=len(keys[1]), shared_edges=len(common),
                left_only=len(keys[0])-len(common), right_only=len(keys[1])-len(common),
                shared_weights_changed=int(np.count_nonzero(delta)),
                max_shared_weight_difference=float(delta.max()) if len(delta) else 0.)


def run(repo, output):
    if output.exists():
        raise FileExistsError(output)
    path = repo / "benchmark_tools/results/installed_ob_search_comparison_20260926.json"
    if record(path)["sha256"] != SEARCH_SHA:
        raise ValueError("Changed retained search comparison")
    previous = json.loads(path.read_text())
    records = [record(path), *previous["checked_records"]]
    for item in records:
        check(item)
    from orthohmm import accuracy, helpers
    for module, expected in ((accuracy, CORE_SHA), (helpers,
            "ca0f4b572480bfb046399a4e87de605dd22ca7b211aa3e9f14da580567f40bc6")):
        source = record(module.__file__)
        if source["sha256"] != expected:
            raise ValueError("Changed graph implementation")
        records.append(source)
    checkpoint = next(Path(r["path"]).parent for r in records
                      if Path(r["path"]).name == "hit_queries.npy")
    cache = next(Path(r["path"]) for r in records if Path(r["path"]).suffix == ".pkl")
    # The parent receipt pins this previously admitted local pickle.
    with cache.open("rb") as stream:
        payload = pickle.load(stream)
    names = sorted(payload["all_gene_ids"])
    ids = {name: i for i, name in enumerate(names)}
    if len(names) != len(ids):
        raise ValueError("Duplicate historical gene names")
    species_labels = sorted(set(payload["gene_to_species"].values()))
    species_ids = {name: i for i, name in enumerate(species_labels)}
    species = np.array([species_ids[payload["gene_to_species"][g]] for g in names], dtype=np.int32)
    hits = payload["all_hits"]
    old_q = np.fromiter((ids[a] for a, b in hits), dtype=np.int32, count=len(hits))
    old_t = np.fromiter((ids[b] for a, b in hits), dtype=np.int32, count=len(hits))
    old_s = np.fromiter(hits.values(), dtype=np.float64, count=len(hits))
    new_names = (checkpoint / "gene_names.txt").read_text().splitlines()
    if set(new_names) != set(names) or len(new_names) != len(names):
        raise ValueError("Fresh gene universe mismatch")
    mapping = np.array([ids[g] for g in new_names], dtype=np.int32)
    fresh_species = np.load(checkpoint / "gene_to_species.npy", allow_pickle=False)
    if fresh_species.shape != mapping.shape:
        raise ValueError("Fresh species array shape differs")
    species_pairs = set(zip(species[mapping].tolist(), fresh_species.tolist()))
    if (len(species_pairs) != len(species_labels)
            or len({b for a, b in species_pairs}) != len(species_labels)):
        raise ValueError("Species memberships differ")
    q, t, s = [np.load(checkpoint / f"hit_{k}.npy", allow_pickle=False, mmap_mode="r")
               for k in ("queries", "targets", "scores")]
    keep = q != t
    new_q, new_t, new_s = mapping[q[keep]], mapping[t[keep]], np.asarray(s[keep])
    old_in_new, new_in_old = align_scores(len(names), old_q, old_t, old_s, new_q, new_t, new_s)
    del hits, payload
    output.mkdir(parents=True)
    arms = {"historical_order_historical_scores": (old_q, old_t, old_s),
            "historical_order_fresh_scores": (old_q, old_t, new_in_old),
            "fresh_order_historical_scores": (new_q, new_t, old_in_new),
            "fresh_order_fresh_scores": (new_q, new_t, new_s)}
    graphs, thresholds, artifacts = {}, {}, {}
    for name, arrays in arms.items():
        graph, threshold = capture_edges(accuracy.build_rbnh_edges, names, species, *arrays)
        graphs[name], thresholds[name] = graph, threshold
        target = output / f"{name}.npz"
        np.savez(target, sources=graph.sources, targets=graph.targets, weights=graph.weights,
                 thresholds=threshold)
        artifacts[name] = record(target)
    name_path = output / "gene_names.txt"
    name_path.write_text("\n".join(names) + "\n")
    contrasts = []
    for a, b in itertools.combinations(arms, 2):
        ta, tb = thresholds[a], thresholds[b]
        finite = np.isfinite(ta) & np.isfinite(tb)
        delta = np.abs(ta[finite] - tb[finite])
        contrasts.append(dict(left=a, right=b, **edge_comparison(graphs[a], graphs[b], len(names)),
            threshold_genes_changed=int(np.count_nonzero(ta != tb)),
            finite_threshold_presence_changed=int(np.count_nonzero(np.isfinite(ta) != np.isfinite(tb))),
            max_shared_finite_threshold_difference=float(delta.max()) if len(delta) else 0.))
    for item in records:
        check(item)
    result = dict(status="initial_graph_score_order_factorial_complete", source=record(__file__),
        checked_records=records, artifacts=artifacts, gene_names=record(name_path), contrasts=contrasts,
        numpy=np.__version__, python=sys.version, genes=len(names), hits=len(old_q),
        comparisons_use_fixed_gene_indexing=True, clustering_run=False, historical_scores_replaced=False,
        limitations=["Post-hoc diagnostic, not a prespecified accuracy comparison.",
            "One frozen graph implementation and audit runtime for all four arms, not historical runtime reconstruction.",
            "Input order and scores varied independently; self hits removed in all arms.",
            "Initial graph excludes singleton assignment, profile expansion, refinement and candidate merging.",
            "Does not quantify final-score causality or test clustering dependency and gene-index effects."])
    (output / "report.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.repo.resolve(), args.output.resolve())
