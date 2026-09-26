"""Isolated installed-production graph inference from label-free numeric hits."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import time


def record(path):
    path = Path(path).resolve()
    data = path.read_bytes()
    return dict(path=str(path), bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def validate_numeric(data):
    fields = {"gene_names", "gene_to_species", "hit_queries", "hit_targets", "hit_scores"}
    if set(data) != fields:
        raise ValueError("Unexpected numeric input fields")
    names, species = data["gene_names"], data["gene_to_species"]
    if not names or any(not isinstance(g, str) or not g or any(c.isspace() for c in g) for g in names):
        raise ValueError("Invalid gene name")
    if names != sorted(set(names)) or len(species) != len(names):
        raise ValueError("Gene/species universe mismatch")
    if any(type(s) is not int or s < 0 for s in species) or set(species) != set(range(max(species) + 1)):
        raise ValueError("Invalid species indices")
    q, t, scores = (data[k] for k in ("hit_queries", "hit_targets", "hit_scores"))
    if len({len(q), len(t), len(scores)}) != 1:
        raise ValueError("Hit lengths differ")
    if any(type(i) is not int or not 0 <= i < len(names) for i in q + t):
        raise ValueError("Invalid gene index")
    pairs = list(zip(q, t))
    if pairs != sorted(set(pairs)):
        raise ValueError("Repeated or unordered directed hits")
    if any(isinstance(s, bool) or not isinstance(s, (int, float)) or not math.isfinite(s) or s <= 0 for s in scores):
        raise ValueError("Invalid graph score")


def write_partition(path, clusters, names):
    flat = [int(i) for group in clusters for i in group]
    if any(not group for group in clusters) or sorted(flat) != list(range(len(names))):
        raise ValueError("Partition must contain every gene exactly once")
    groups = sorted(tuple(sorted(names[int(i)] for i in group)) for group in clusters)
    with path.open("x") as stream:
        for group in groups:
            stream.write("\t".join(group) + "\n")


def run(numeric, output):
    if not sys.flags.isolated:
        raise RuntimeError("Use the installed Python with -I")
    import numpy as np
    from orthohmm import accuracy, externals, refinement

    modules = [accuracy, externals, refinement]
    if any(not Path(m.__file__).resolve().is_relative_to(Path(sys.prefix).resolve()) for m in modules):
        raise RuntimeError("Production module outside installed prefix")
    before = record(numeric)
    data = json.loads(numeric.read_text())
    validate_numeric(data)
    output.mkdir(parents=True, exist_ok=False)
    names = data["gene_names"]
    species, q, t = [np.asarray(data[k], dtype=np.int32) for k in ("gene_to_species", "hit_queries", "hit_targets")]
    scores = np.asarray(data["hit_scores"], dtype=np.float64)
    checkpoint = accuracy.write_accuracy_checkpoint(str(output), names, species, q, t, scores)
    times = {}
    started = time.monotonic()
    edges = accuracy.build_rbnh_edges(names, species, q, t, scores)
    np.savez(output / "rbnh_edges.npz", sources=edges.sources, targets=edges.targets, weights=edges.weights)
    externals.execute_leiden(.1, str(output), edges=edges, include_isolates=True, seed=4)
    clustered = output / "orthohmm_working_res/orthohmm_edges_clustered.txt"
    initial = accuracy.read_index_clusters(str(clustered), names)
    write_partition(output / "initial.tsv", initial, names)
    times["initial_graph_and_clustering_s"] = time.monotonic() - started
    started = time.monotonic()
    singleton = accuracy.build_singleton_assignment_edges(names, initial, q, t, scores)
    combined = accuracy.combine_edges(edges, singleton)
    np.savez(output / "multipass_edges.npz", sources=combined.sources, targets=combined.targets, weights=combined.weights)
    externals.execute_leiden(.1, str(output), edges=combined, include_isolates=True, seed=4)
    multipass = accuracy.read_index_clusters(str(clustered), names)
    write_partition(output / "multipass.tsv", multipass, names)
    times["singleton_graph_and_clustering_s"] = time.monotonic() - started
    started = time.monotonic()
    rq, rt, rs = ([], [], []) if len(np.unique(species)) >= refinement.DEFAULT_COPY_SPLIT_MIN_DATASET_SPECIES else (q, t, scores)
    final = refinement.refine_cluster_indices(multipass, rq, rt, rs, combined.sources,
                                             combined.targets, species, rbnh_scores=combined.weights)
    write_partition(output / "final.tsv", final, names)
    times["refinement_s"] = time.monotonic() - started
    if record(numeric) != before:
        raise ValueError("Numeric input changed during inference")
    receipt = dict(status="native_graph_completed_pending_independent_readback", numeric=before,
                   source=record(__file__), modules=[record(m.__file__) for m in modules],
                   checkpoint_manifest=record(checkpoint / "manifest.json"),
                   outputs=[record(output / name) for name in ("rbnh_edges.npz", "multipass_edges.npz", "initial.tsv", "multipass.tsv", "final.tsv")],
                   settings=dict(cpm_resolution=.1, leiden_seed=4, include_isolates=True,
                                 refinement=True, profile_expansion=False, candidate_expansion=False, phylogeny=False),
                   genes=len(names), species=len(set(species)), hits=len(scores),
                   groups=len(final), timings=times, isolated=True, truth_loaded=False,
                   executable=sys.executable, prefix=sys.prefix)
    with (output / "receipt.json").open("x") as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True)
        stream.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--numeric", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.numeric.resolve(), args.output.absolute())
