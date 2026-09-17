"""Isolate direct graph construction and weight assignment; no optimizer calls."""

import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

ADMISSION_SHA = "499f90e06b6da311356c80f435149ba0cf132f75b13e2eacf720e02ea8f2382e"
MODES = ("minimal_imports", "frozen_imports")
FORMATS = ("numpy", "python_pairs")


def edge_argument(edges, edge_format):
    if edge_format == "numpy":
        return edges
    if edge_format == "python_pairs":
        return ((int(a), int(b)) for a, b in edges)
    raise ValueError("Unknown edge input format")


def planned_workers(compare_formats=False):
    return [(index, mode, edge_format) for index in range(3)
            for mode, edge_format in ([("minimal_imports", fmt) for fmt in FORMATS]
                                     if compare_formats else [(mode, "numpy") for mode in MODES])]


def worker(root, payload, mode, edge_format="numpy"):
    from repeat_qfo_saved_graph import record, mapped_libraries, ENV_KEYS, set_worker_affinity
    from probe_leiden_boundary import endpoint_differences, graph_fingerprint, saved_fingerprint
    if mode not in MODES:
        raise ValueError("Unknown import mode")
    inherited = set_worker_affinity([min(os.sched_getaffinity(0))])
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    sys.path.insert(0, str(launcher))
    if mode == "frozen_imports":
        from orthohmm import leiden_worker
        import leidenalg
        if Path(leiden_worker.__file__).resolve() != launcher / "orthohmm/leiden_worker.py":
            raise ValueError("Wrong frozen worker import")
        def forbidden(*args, **kwargs):
            raise RuntimeError("Optimizer must not run in direct construction diagnostic")
        leidenalg.find_partition = forbidden
    import igraph
    import numpy as np
    if mode == "minimal_imports" and any(name.split(".")[0] in {"orthohmm", "leidenalg"} for name in sys.modules):
        raise ValueError("Minimal worker imported scientific worker or optimizer")
    inputs = [record(payload / name) for name in ("gene_names.txt", "sources.npy", "targets.npy", "weights.npy")]
    names = (payload / "gene_names.txt").read_text().splitlines()
    arrays = [np.load(payload / (name + ".npy"), mmap_mode="r", allow_pickle=False)
              for name in ("sources", "targets", "weights")]
    if (len(set(names)) != len(names) or not names or any(a.ndim != 1 for a in arrays)
            or len({len(a) for a in arrays}) != 1 or arrays[0].dtype != np.int32 or arrays[1].dtype != np.int32
            or arrays[2].dtype != np.float64 or not np.isfinite(arrays[2]).all()
            or any(np.any((a < 0) | (a >= len(names))) for a in arrays[:2])):
        raise ValueError("Invalid direct-construction inputs")
    used_ids = np.arange(len(names), dtype=np.int32)
    local_sources = np.searchsorted(used_ids, arrays[0]).astype(np.int32)
    local_targets = np.searchsorted(used_ids, arrays[1]).astype(np.int32)
    edges = np.column_stack((local_sources, local_targets))
    weights = np.asarray(arrays[2], dtype=np.float64)
    snapshot = {"mode": mode, "edge_format": edge_format, "source": record(__file__), "inputs": inputs,
        "inherited_cpu_affinity": inherited, "cpu_affinity": sorted(os.sched_getaffinity(0)),
        "cwd": str(Path.cwd()), "host": platform.node(), "platform": platform.platform(),
        "python": record(sys.executable), "environment": {key: os.environ.get(key) for key in ENV_KEYS},
        "versions": {name: importlib.metadata.version(name) for name in ("numpy", "igraph")},
        "modules": {name: record(module.__file__) for name, module in list(sys.modules.items())
            if name.split(".")[0] in {"numpy", "igraph", "leidenalg", "orthohmm"} and getattr(module, "__file__", None)},
        "native_libraries": [record(path) for path in mapped_libraries(Path("/proc/self/maps").read_text())],
        "accuracy_evaluated": False, "optimizer_called": False}
    (payload / "worker_before.json").write_text(json.dumps(snapshot, indent=2, sort_keys=True) + "\n")
    graph = igraph.Graph(n=len(names), edges=edge_argument(edges, edge_format), directed=False)
    before = {"vertices": graph.vcount(), "edges": graph.ecount(), "directed": graph.is_directed(),
              "differences": endpoint_differences(graph, payload, edges)}
    (payload / "before_weights.json").write_text(json.dumps(before, indent=2, sort_keys=True) + "\n")
    graph.es["weight"] = weights
    after = {"fingerprint": graph_fingerprint(graph), "differences": endpoint_differences(graph, payload, edges)}
    saved = saved_fingerprint(payload)
    if [record(payload / name) for name in ("gene_names.txt", "sources.npy", "targets.npy", "weights.npy")] != inputs:
        raise ValueError("Saved inputs changed during direct construction")
    result = {"status": "direct_construction_observed", "mode": mode, "edge_format": edge_format, "before_weights": before,
              "after_weights": after, "saved": saved, "snapshot": record(payload / "worker_before.json"),
              "accuracy_evaluated": False, "optimizer_called": False}
    (payload / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--worker-payload", type=Path)
    parser.add_argument("--mode", choices=MODES)
    parser.add_argument("--edge-format", choices=FORMATS, default="numpy")
    parser.add_argument("--compare-formats", action="store_true", help="Compare NumPy and Python-pair input using minimal imports")
    args = parser.parse_args()
    root, output = args.root.resolve(), args.output.resolve()
    if args.worker_payload:
        worker(root, args.worker_payload.resolve(), args.mode, args.edge_format)
        return
    if output.exists():
        raise FileExistsError(output)
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
    from benchmark_tools.run_simulation_methods import read_frozen
    from benchmark_tools.verify_qfo_replay_launcher import verify
    admission_path = root / "benchmark_tools/results/qfo_construction_verified_20260916.json"
    admission = read_frozen(admission_path, ADMISSION_SHA)
    if admission["status"] != "construction_observations_verified":
        raise ValueError("Previous construction panel not admitted")
    for item in admission["provenance_checked"]:
        check(item)
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    frozen = root / "benchmarks/work/publication_method_native_v2"
    runtime_path = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    runtime = verify(frozen, launcher, runtime_path)
    output.mkdir(parents=True)
    report = {"status": "running", "job_id": os.environ.get("SLURM_JOB_ID"), "source": record(__file__),
        "admission": record(admission_path), "runtime": runtime,
        "graph_inputs": admission["native_report"]["graph_inputs"], "workers": [],
        "planned_workers": planned_workers(args.compare_formats), "compare_formats": args.compare_formats,
        "accuracy_evaluated": False, "optimizer_called": False,
        "limitations": ["Direct construction omits frozen worker bookkeeping; before-weight observation changes allocation/timing.",
            "Three fresh workers per selected arm; original int32 array retained for comparisons; not proof of general correctness.",
            "Python-pair input changes conversion and allocation, not a proven fix or isolated causal mechanism.",
            "No optimizer, partitions, accuracy scores or default changes."]}
    env = os.environ.copy()
    env.update(PYTHONPATH=str(launcher), PYTHONHASHSEED="0", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    try:
        for index, mode, edge_format in planned_workers(args.compare_formats):
            label = f"{edge_format}_{index}" if args.compare_formats else f"{mode}_{index}"
            payload = output / label / "payload"
            payload.mkdir(parents=True)
            for item in report["graph_inputs"]:
                source = Path(item["path"])
                (payload / source.name.removeprefix("rbnh_")).symlink_to(source)
            command = [sys.executable, str(Path(__file__).resolve()), "--root", str(root), "--output", str(output),
                       "--worker-payload", str(payload), "--mode", mode, "--edge-format", edge_format]
            with (payload.parent / "worker.log").open("x") as log:
                completed = subprocess.run(command, cwd=launcher, env=env, stdout=log, stderr=subprocess.STDOUT)
            if completed.returncode:
                raise RuntimeError("Direct-construction worker failed: " + str(payload))
            result = json.loads((payload / "result.json").read_text())
            report["workers"].append({"index": index, "mode": mode, "edge_format": edge_format, "exit_code": completed.returncode,
                "result": result, "observation": record(payload / "result.json")})
            (output / "progress.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        if verify(frozen, launcher, runtime_path) != runtime:
            raise ValueError("Frozen runtime changed")
        for item in admission["provenance_checked"]:
            check(item)
        report["status"] = "six_direct_workers_complete_unscored"
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        (output / "results.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
