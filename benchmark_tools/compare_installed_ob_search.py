"""Compare retained normalized search evidence without rerunning inference."""

import argparse
from collections import Counter
import json
import math
from pathlib import Path
import pickle

import numpy as np

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.audit_installed_orthobench import compare_partitions

CACHE_SHA = "78a5af40ea2683a69e1baefefb0966c3549cffe71b24bda03b9ba4a6e29e65e1"
READBACK_SHA = "79fc6ce6a169d486235c811c7d6eec9019127958c07e56d8ff982675cd4149e0"


def compare_hits(historical, names, queries, targets, scores):
    """Consume a historical dictionary; collapse fresh duplicates by max score."""
    n = len(names)
    if not n or len(set(names)) != n or any(not name for name in names):
        raise ValueError("Invalid gene names")
    if any(a.ndim != 1 for a in (queries, targets, scores)) or not (
            len(queries) == len(targets) == len(scores)):
        raise ValueError("Hit array dimensions differ")
    if any(a.dtype.kind not in "iu" for a in (queries, targets)):
        raise ValueError("Noninteger hit indices")
    if len(queries) and (min(queries.min(), targets.min()) < 0
                        or max(queries.max(), targets.max()) >= n):
        raise ValueError("Out-of-universe hit index")
    if not np.all(np.isfinite(scores)) or np.any(scores <= 0):
        raise ValueError("Invalid fresh score")
    universe = set(names)
    for pair, score in historical.items():
        if (len(pair) != 2 or pair[0] not in universe or pair[1] not in universe
                or pair[0] == pair[1] or not math.isfinite(score) or score <= 0):
            raise ValueError("Invalid historical normalized hit")
    totals = Counter(historical_unique=len(historical), fresh_rows=len(queries))
    keys = queries.astype(np.int64) * n + targets
    order = np.argsort(keys, kind="stable")
    witnesses = {}
    max_abs = max_rel = 0.0

    def compare(key, value):
        nonlocal max_abs, max_rel
        q, t = divmod(int(key), n)
        pair = (names[q], names[t])
        old = historical.pop(pair, None)
        totals["fresh_unique_nonself"] += 1
        if old is None:
            category = "fresh_only"
        else:
            totals["shared"] += 1
            delta = abs(value - old)
            max_abs = max(max_abs, delta)
            max_rel = max(max_rel, delta / abs(old))
            category = "shared_exact" if value == old else "shared_changed"
            if not math.isclose(value, old, rel_tol=1e-12, abs_tol=1e-12):
                totals["shared_outside_1e12_tolerance"] += 1
        totals[category] += 1
        if category != "shared_exact" and category not in witnesses:
            witnesses[category] = dict(pair=pair, historical=old, fresh=value)

    last_key, best = None, None
    for i in order:
        if queries[i] == targets[i]:
            totals["fresh_self_rows_excluded"] += 1
            continue
        key, score = int(keys[i]), float(scores[i])
        if key == last_key:
            totals["fresh_duplicate_nonself_rows"] += 1
            best = max(best, score)
        else:
            if last_key is not None:
                compare(last_key, best)
            last_key, best = key, score
    if last_key is not None:
        compare(last_key, best)
    totals["historical_only"] = len(historical)
    if historical:
        pair = min(historical)
        witnesses["historical_only"] = dict(pair=pair, historical=historical[pair], fresh=None)
    return dict(counts=dict(totals), max_shared_absolute_difference=max_abs,
                max_shared_relative_difference=max_rel, witnesses=witnesses,
                directed_nonself_presence_equal=not (totals["fresh_only"] or totals["historical_only"]),
                normalized_scores_exact=not (totals["fresh_only"] or totals["historical_only"]
                                            or totals["shared_changed"]))


def partition(path, universe):
    groups, seen = [], set()
    for line in path.read_text().splitlines():
        genes = line.split()
        group = frozenset(genes)
        if not genes or len(group) != len(genes) or seen & group or not group <= universe:
            raise ValueError("Invalid candidate partition")
        seen.update(group)
        groups.append(group)
    if seen != universe:
        raise ValueError("Incomplete candidate partition")
    return groups


def audit(repo):
    readback = repo / "benchmark_tools/results/installed_orthobench_readback_20260926.json"
    if record(readback)["sha256"] != READBACK_SHA:
        raise ValueError("Changed installed readback")
    receipt = json.loads(readback.read_text())
    records = [record(readback), receipt["execution_record"]]
    for item in records:
        check(item)
    run = Path(receipt["execution_record"]["path"]).parent
    checkpoint = run / "inference/orthohmm_working_res/high_sensitivity_checkpoint"
    manifest = json.loads((checkpoint / "manifest.json").read_text())
    if manifest["complete"] is not True or manifest["schema_version"] != 1:
        raise ValueError("Incomplete checkpoint")
    records.append(record(checkpoint / "manifest.json"))
    if set(manifest["files"]) != {"gene_names.txt", "gene_to_species.npy", "hit_queries.npy",
                                 "hit_targets.npy", "hit_scores.npy"}:
        raise ValueError("Incomplete checkpoint inventory")
    for name, item in manifest["files"].items():
        if name not in {"gene_names.txt", "gene_to_species.npy", "hit_queries.npy",
                        "hit_targets.npy", "hit_scores.npy"}:
            raise ValueError("Unexpected checkpoint filename")
        expected = dict(item, path=str(checkpoint / name))
        check(expected)
        records.append(expected)
    names = (checkpoint / "gene_names.txt").read_text().splitlines()
    arrays = [np.load(checkpoint / f"hit_{n}.npy", mmap_mode="r", allow_pickle=False)
              for n in ("queries", "targets", "scores")]
    if len(names) != manifest["genes"] or len(arrays[0]) != manifest["hits"]:
        raise ValueError("Checkpoint count mismatch")
    cache = repo / "benchmarks/results/hits_BLOSUM62_mc100.pkl"
    cache_record = record(cache)
    if cache_record["sha256"] != CACHE_SHA:
        raise ValueError("Untrusted historical pickle")
    records.append(cache_record)
    # Only the previously admitted, checksum-pinned local cache is deserialized.
    with cache.open("rb") as stream:
        historical = pickle.load(stream)
    if set(historical["all_gene_ids"]) != set(names):
        raise ValueError("Historical gene universe differs")
    comparison = compare_hits(historical["all_hits"], names, *arrays)
    old = repo / "benchmarks/results/publication_ob_factorial_v1/candidates/p1_c1/orthohmm_working_res/phylogeny_candidate_superfamilies.txt"
    new = run / "inference/orthohmm_working_res/phylogeny_candidate_superfamilies.txt"
    records.extend([record(old), record(new)])
    if records[-2]["sha256"] != "44f201dbc7e4eccf06d9401e5ce1b20fdf6ad523e15ee6007fe1c58f3941842d":
        raise ValueError("Changed historical candidate partition")
    candidates = compare_partitions(partition(old, set(names)), partition(new, set(names)))
    for item in records:
        check(item)
    return dict(status="retained_search_and_candidate_comparison_complete", source=record(__file__),
                checked_records=records, search=comparison, candidate_partitions=candidates,
                native_inference_rerun=False, scores_replaced=False,
                limitations=["Post-hoc provenance diagnostic, not a new accuracy endpoint.",
                    "Compares accepted normalized directed nonself hits; absent raw candidates cannot identify a filter/scorer cause.",
                    "Historical cache and fresh checkpoint are different execution paths/environments.",
                    "Differences before clustering do not quantify their causal contribution to final score changes.",
                    "No raw-score, E-value or complete upstream runtime equivalence is established."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.repo.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
