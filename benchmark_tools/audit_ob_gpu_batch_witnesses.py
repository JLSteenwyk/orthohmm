"""Find accepted short-target witnesses for the all-long GPU routing condition."""

import argparse
import itertools
import json
import math
from pathlib import Path
import pickle
import sys

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.trace_ob_initial_edges import TRACE_SHA, CACHE_SHA


def witnesses(hits, owners, lengths):
    if set(owners) != set(lengths) or any(type(n) is not int or n <= 0 for n in lengths.values()):
        raise ValueError("Invalid gene length universe")
    labels = sorted(set(owners.values()))
    rows = {(a, b): {"query_species": a, "target_species": b, "retained_hits": 0,
                     "eligible_target_hits": 0, "witness": None}
            for a, b in itertools.product(labels, repeat=2)}
    for (query, target), score in hits.items():
        if query not in owners or target not in owners or not math.isfinite(score) or score <= 0:
            raise ValueError("Invalid retained hit")
        row = rows[owners[query], owners[target]]
        row["retained_hits"] += 1
        if lengths[target] <= 1998:
            row["eligible_target_hits"] += 1
            candidate = (lengths[target], target, query)
            previous = row["witness"]
            if previous is None or candidate < (previous["target_length"], previous["target"], previous["query"]):
                row["witness"] = {"query": query, "target": target,
                                  "target_length": lengths[target], "retained_score": float(score)}
    for row in rows.values():
        row["assessment"] = ("eligible_target_witness_present" if row["witness"] is not None
                             else "unresolved_no_eligible_retained_target")
    return list(rows.values())


def run(root, output):
    if output.exists():
        raise FileExistsError(output)
    trace_path = root / "benchmark_tools/results/ob_family_trace_verified_20260916.json"
    trace = read_frozen(trace_path, TRACE_SHA)
    cache = [r for r in trace["inputs"] if r["sha256"] == CACHE_SHA]
    fastas = [r for r in trace["inputs"] if Path(r["path"]).suffix == ".fa"]
    if len(cache) != 1 or len(fastas) != 12:
        raise ValueError("Unexpected admitted input inventory")
    checked = [record(trace_path), *cache, *fastas,
               record(root / "benchmarks/benchmark_orthobench.py")]
    for item in checked:
        check(item)
    owners, lengths = {}, {}
    for item in fastas:
        path = Path(item["path"])
        for seq in SeqIO.parse(path, "fasta"):
            if seq.id in owners:
                raise ValueError("Duplicate gene ID")
            owners[seq.id], lengths[seq.id] = path.name, len(seq.seq)
    # Only deserialize the admitted, checksum-pinned local cache.
    with Path(cache[0]["path"]).open("rb") as stream:
        data = pickle.load(stream)
    if (len(owners) != 251378 or owners != data["gene_to_species"]
            or lengths != data["gene_length_dict"]
            or len(data["all_gene_ids"]) != len(owners)
            or set(data["all_gene_ids"]) != set(owners)):
        raise ValueError("Cache universe/ownership/lengths differ from admitted FASTAs")
    rows = witnesses(data["all_hits"], owners, lengths)
    for item in checked:
        check(item)
    result = {"status": "retained_hit_target_length_witnesses_audited", "source": record(__file__),
              "checked_records": checked, "genes": len(owners), "species": len(fastas),
              "hits": len(data["all_hits"]), "directions": rows,
              "witnessed_directions": sum(r["witness"] is not None for r in rows),
              "unresolved_directions": sum(r["witness"] is None for r in rows),
              "accuracy_evaluated": False, "publication_ready": False,
              "limitations": ["Conditional on one complete species-pair candidate batch and authentic retained search outputs.",
                  "The retained driver uses this batching, but its current source is not proof of the exact historical executed revision.",
                  "An accepted eligible target rules out the all-long condition for that batch; no witness is unresolved, not proof of failure.",
                  "Does not establish GPU availability, rejected-candidate identity, score correctness or historical runtime integrity.",
                  "No inference rerun, reference-label accuracy evaluation or other-dataset non-impact claim."]}
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    run(args.root.resolve(), args.output.resolve())
