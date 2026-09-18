"""Bound one frozen RBNH named-array snapshot, not total peak RAM."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_accuracy_checkpoint import audit
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_simulation_methods import read_frozen

CORE_SHA = "1a35944ab7fea859143f599b1262aa111272787f6f735b9acc37d299427f2ad6"


def bounds(genes, species_slots, hits, self_hits):
    values = (genes, species_slots, hits, self_hits)
    if any(type(v) is not int for v in values) or genes < 1 or species_slots < 1 or hits < 0 or not 0 <= self_hits <= hits:
        raise ValueError("Invalid numeric checkpoint dimensions")
    eligible = hits - self_hits
    slots = genes * species_slots
    # At the line after best_targets = np.full(...), all these locals coexist.
    fixed = {"eligible": hits, "best_queries": 4 * eligible,
             "best_targets_input": 4 * eligible, "best_input_scores": 8 * eligible,
             "keys": 8 * eligible, "best_scores": 8 * slots, "best_targets": 4 * slots}
    active = eligible > 0
    base = sum(fixed.values()) if active else 0
    # W: tied winning rows; K: occupied query/target-species slots. 1 <= K <= W <= E.
    winner_min = 32 if active else 0
    winner_max = 8 * eligible + 24 * min(eligible, slots) if active else 0
    return {"genes": genes, "species_slot_extent": species_slots, "hits": hits,
            "self_hits": self_hits, "eligible_nonself_finite_hits": eligible,
            "gene_species_slots": slots, "snapshot_reached": active,
            "fixed_named_array_bytes": fixed if active else {},
            "winner_array_expression_bytes": "8*W + 24*K",
            "named_array_payload_lower_bytes": base + winner_min,
            "named_array_payload_upper_bytes_at_this_snapshot_only": base + winner_max,
            "input_array_logical_bytes_separate": 16 * hits + 4 * genes,
            "total_peak_ram_bound_available": False, "graph_feasibility_admitted": False}


def estimate(checkpoint, checkpoint_sha, core):
    if np.dtype(np.intp).itemsize != 8:
        raise ValueError("This allocation model requires 64-bit NumPy indices")
    core_record = record(core)
    if core_record["sha256"] != CORE_SHA:
        raise ValueError("Memory model requires the frozen accuracy.py bytes")
    before = audit(checkpoint, checkpoint_sha)
    manifest = read_frozen(checkpoint / "manifest.json", checkpoint_sha)
    inputs = [record(checkpoint / name) for name in sorted(manifest["files"])]
    for item in inputs:
        expected = manifest["files"][Path(item["path"]).name]
        if any(item[key] != expected[key] for key in ("bytes", "sha256")):
            raise ValueError("Checkpoint changed after numeric audit")
    species = np.load(checkpoint / "gene_to_species.npy", mmap_mode="r", allow_pickle=False)
    summary = before["summary"]
    result = bounds(summary["genes"], int(species.max()) + 1, summary["hits"], summary["self_hits"])
    for item in [core_record, before["manifest"], *inputs]:
        check(item)
    return {"status": "frozen_rbnh_named_array_payload_bound", "publication_ready": False,
            "checkpoint_audit": before, "checkpoint_files": inputs, "core": core_record,
            "source": record(__file__), "numpy": np.__version__, "estimate": result,
            "snapshot": "build_rbnh_edges: after best_targets allocation, before best_targets[winning_keys] assignment",
            "limitations": [
                "Array payload calculation, not measured RSS, cgroup memory, runtime or an allocation recommendation.",
                "The upper bound covers only the named arrays at this snapshot, not the complete algorithm.",
                "Excludes sorting temporaries, allocator overhead, gene strings, input page residency and shared libraries.",
                "Excludes later deduplication, graph-library construction, clustering, refinement and singleton assignment.",
                "Input logical bytes are reported separately; memory mapping does not guarantee resident memory savings.",
                "Numerical integrity is checked, not source-search equivalence or scientific admission.",
                "No hits are truncated and no inference settings are changed."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--core", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = estimate(args.checkpoint.resolve(), args.sha256, args.core.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps(result["estimate"], indent=2, sort_keys=True))
