"""Validate immutable numeric search checkpoints in bounded-memory chunks."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
from orthohmm.accuracy import load_accuracy_checkpoint
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.orthobench_stage_diagnostics import file_provenance


def check_arrays(names, species, queries, targets, scores, chunk_size=1000000):
    if not names or len(set(names)) != len(names) or any(not n or n.strip() != n for n in names):
        raise ValueError("Empty, duplicate or malformed gene identifiers")
    if species.shape != (len(names),) or species.dtype != np.dtype("int32") or np.any(species < 0):
        raise ValueError("Invalid species array")
    if queries.ndim != 1 or targets.shape != queries.shape or scores.shape != queries.shape:
        raise ValueError("Hit array shapes differ")
    if queries.dtype != np.dtype("int32") or targets.dtype != np.dtype("int32") or scores.dtype != np.dtype("float64"):
        raise ValueError("Unexpected native hit dtypes")
    if chunk_size < 1:
        raise ValueError("Positive chunk size required")
    self_hits, nonpositive = 0, 0
    minimum = maximum = None
    for start in range(0, len(scores), chunk_size):
        q, t, s = queries[start:start + chunk_size], targets[start:start + chunk_size], scores[start:start + chunk_size]
        if np.any(q < 0) or np.any(t < 0) or np.any(q >= len(names)) or np.any(t >= len(names)):
            raise ValueError("Hit index outside gene table")
        if not np.isfinite(s).all():
            raise ValueError("Nonfinite hit score")
        self_hits += int(np.count_nonzero(q == t))
        nonpositive += int(np.count_nonzero(s <= 0))
        minimum = float(s.min()) if minimum is None else min(minimum, float(s.min()))
        maximum = float(s.max()) if maximum is None else max(maximum, float(s.max()))
    return {"genes": len(names), "hits": len(scores), "species": len(np.unique(species)),
            "gene_names_lexically_sorted": names == sorted(names), "self_hits": self_hits,
            "nonpositive_scores": nonpositive, "minimum_score": minimum, "maximum_score": maximum}


def audit(path, expected_hash):
    manifest_path = path / "manifest.json"
    manifest = read_frozen(manifest_path, expected_hash)
    expected = {"gene_names.txt", "gene_to_species.npy", "hit_queries.npy", "hit_targets.npy", "hit_scores.npy"}
    if set(manifest["files"]) != expected or {p.name for p in path.iterdir()} != expected | {"manifest.json"}:
        raise ValueError("Unexpected checkpoint file inventory")
    arrays = load_accuracy_checkpoint(path, verify=True)
    summary = check_arrays(*arrays)
    if summary["genes"] != manifest["genes"] or summary["hits"] != manifest["hits"]:
        raise ValueError("Manifest counts differ")
    if file_provenance(manifest_path)["sha256"] != expected_hash:
        raise ValueError("Manifest changed during audit")
    return {"schema_version": 1, "status": "numeric_checkpoint_verified", "manifest": file_provenance(manifest_path),
            "summary": summary, "auditor": file_provenance(Path(__file__)), "accuracy_evaluated": False,
            "limitations": "Checks hashes and numeric integrity, not FASTA/source equivalence, hit completeness or duplicate ordered pairs."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.checkpoint.resolve(), args.sha256)
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(result["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
