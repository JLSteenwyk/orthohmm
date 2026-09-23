"""Read-only failure localization and array validation; never reruns clustering."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def arrays(payload, chunk_size=100000):
    if chunk_size < 1:
        raise ValueError("Invalid chunk size")
    names = (payload / "gene_names.txt").read_text().splitlines()
    if not names or len(names) != len(set(names)) or any(not n for n in names):
        raise ValueError("Empty or duplicate gene identities")
    src, dst, weights = [np.load(payload / (name + ".npy"), mmap_mode="r", allow_pickle=False)
                         for name in ("sources", "targets", "weights")]
    if (any(a.ndim != 1 or a.shape != src.shape for a in (src, dst, weights))
            or [a.dtype for a in (src, dst, weights)] != [np.dtype("int32"), np.dtype("int32"), np.dtype("float64")]
            or not len(src)):
        raise ValueError("Invalid array shape or dtype")
    digest = hashlib.sha256()
    negative = zero = self_edges = 0
    minimum, maximum = float("inf"), float("-inf")
    for start in range(0, len(src), chunk_size):
        a, b, w = (x[start:start+chunk_size] for x in (src, dst, weights))
        if any(np.any(x < 0) or np.any(x >= len(names)) for x in (a, b)) or not np.isfinite(w).all():
            raise ValueError("Invalid endpoints or nonfinite weights")
        digest.update(np.column_stack((a, b)).tobytes(order="C"))
        negative += int(np.count_nonzero(w < 0))
        zero += int(np.count_nonzero(w == 0))
        self_edges += int(np.count_nonzero(a == b))
        minimum, maximum = min(minimum, float(w.min())), max(maximum, float(w.max()))
    return dict(vertices=len(names), edges=len(src), constructor_int32_bytes_sha256=digest.hexdigest(),
        negative_weights=negative, zero_weights=zero, self_edges=self_edges,
        minimum_weight=minimum, maximum_weight=maximum)


def audit(root, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    stage = root.resolve() / "benchmarks/results/qfo_parameter_cpm_replay_v3/cpm_high/clustering/cluster_3_profile_expanded"
    payload = stage / "payload"
    paths = [stage / "execution.json", stage / "payload_manifest.json", payload / "worker_before.json",
             payload / "constructor_adapter.json", payload / "metadata.json"]
    execution, manifest, before, adapter, metadata = [json.loads(p.read_text()) for p in paths]
    tracked = [record(p) for p in paths] + [record(__file__), *manifest["inputs"],
        *before["modules"].values(), *before["native_libraries"], before["python"]]
    for item in tracked:
        check(item)
    if (execution["status"] != "failed" or execution["index"] != 3 or "SIGSEGV" not in execution["error"]
            or before["status"] != "before_native_clustering" or before["metadata"] != metadata
            or metadata["cpm_resolution"] != .12 or metadata["seed"] != 4
            or metadata["include_isolates"] is not True or adapter["format"] != "python_pairs"
            or len(adapter["calls"]) != 1 or adapter["calls"][0]["status"] != "before_constructor"
            or (payload / "native_boundary.json").exists()):
        raise ValueError("Failure evidence no longer matches diagnostic target")
    observed = arrays(payload)
    call = adapter["calls"][0]
    if (observed["vertices"] != 984137 or observed["edges"] != 25501180
            or call["n"] != observed["vertices"] or call["shape"] != [observed["edges"], 2]
            or call["ordered_input_bytes_sha256"] != observed["constructor_int32_bytes_sha256"]
            or call["directed"] is not False):
        raise ValueError("Constructor input disagrees with saved arrays")
    for item in tracked:
        check(item)
    report = dict(status="failed_high_cpm_payload_checked_without_execution", observed=observed,
        checked_records=tracked, versions=before["versions"], environment=before["environment"],
        cpu_affinity=before["cpu_affinity"], accuracy_evaluated=False, rerun_authorized=False,
        publication_ready=False, limitations=[
            "Last marker is before graph construction; no optimizer marker exists. This does not identify the faulting native instruction.",
            "Current files match saved identities; this cannot prove historical process memory was uncorrupted.",
            "No graph construction, optimizer, recovery or accuracy scoring runs in this audit."])
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    audit(args.root, args.output)
