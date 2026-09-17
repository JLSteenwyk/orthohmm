"""Bounded fresh-worker graph construction tests; never optimize or score."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

PREVIOUS_SHA = "f0f247e17c63235df36038a9f8a978e1d30b137fcdadc4ea0c2a1768f78b8b50"


def converted_edges(edges, mode):
    import numpy as np
    if mode == "original_int32":
        if edges.dtype != np.int32:
            raise ValueError("Expected original int32 constructor array")
        return edges
    if mode == "explicit_int64":
        return np.array(edges, dtype=np.int64, order="C", copy=True)
    raise ValueError("Unknown constructor mode")


def construction_worker(root, payload, mode):
    from repeat_qfo_saved_graph import worker, set_worker_affinity
    from probe_leiden_boundary import graph_fingerprint, saved_fingerprint, endpoint_differences
    set_worker_affinity([min(os.sched_getaffinity(0))])
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    sys.path.insert(0, str(launcher))
    import igraph
    import leidenalg
    constructor = igraph.Graph.__init__
    captured = {}
    class Observed(Exception):
        pass
    def construct(graph, *args, **kwargs):
        original = kwargs["edges"]
        converted = converted_edges(original, mode)
        constructor(graph, *args, **{**kwargs, "edges": converted})
        captured.update(original=original, converted=converted)
    def observe(graph, *args, **kwargs):
        result = {"status": "construction_observed_without_optimization", "mode": mode,
            "native": graph_fingerprint(graph), "saved": saved_fingerprint(payload),
            "differences": endpoint_differences(graph, payload, captured["original"]),
            "converted_differences": endpoint_differences(graph, payload, captured["converted"]),
            "converted_dtype": str(captured["converted"].dtype), "optimizer_called": False,
            "accuracy_evaluated": False, "witnesses": []}
        for example in result["differences"]["native_vs_saved"]["examples"]:
            edge = graph.es[example["edge_index"]]
            expected = example["right"]
            result["witnesses"].append({"edge_index": example["edge_index"],
                "tuple": list(edge.tuple), "source_target": [edge.source, edge.target],
                "expected_pair_edge_id": graph.get_eid(*expected, directed=False, error=False)})
        (payload / "construction.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        raise Observed()
    igraph.Graph.__init__ = construct
    leidenalg.find_partition = observe
    try:
        worker(launcher, payload)
    except Observed:
        sys.stdout.flush()
        sys.stderr.flush()
        os._exit(0)
    raise RuntimeError("Construction observer was not reached")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--worker-payload", type=Path)
    parser.add_argument("--mode", choices=("original_int32", "explicit_int64"))
    args = parser.parse_args()
    root, output = args.root.resolve(), args.output.resolve()
    if args.worker_payload:
        construction_worker(root, args.worker_payload.resolve(), args.mode)
        return
    if output.exists():
        raise FileExistsError(output)
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
    from benchmark_tools.verify_qfo_replay_launcher import verify
    from benchmark_tools.run_simulation_methods import read_frozen
    previous_path = root / "benchmarks/results/qfo_native_boundary_v2/results.json"
    previous = read_frozen(previous_path, PREVIOUS_SHA)
    for item in previous["graph_inputs"]:
        check(item)
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    frozen = root / "benchmarks/work/publication_method_native_v2"
    runtime = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    before = verify(frozen, launcher, runtime)
    output.mkdir(parents=True)
    report = {"status": "running", "job_id": os.environ.get("SLURM_JOB_ID"), "source": record(__file__),
        "previous": record(previous_path), "graph_inputs": previous["graph_inputs"], "runtime": before,
        "optimizer_called": False, "accuracy_evaluated": False, "workers": [],
        "limitations": ["Three alternating fresh workers per dtype; construction-only diagnostic, not proof of general determinism.",
            "Explicit int64 conversion changes allocation as well as dtype; no default fix inferred from this panel."]}
    env = os.environ.copy()
    env.update(PYTHONPATH=str(launcher), PYTHONHASHSEED="0", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    try:
        for index in range(3):
            for mode in ("original_int32", "explicit_int64"):
                directory = output / f"{mode}_{index}"
                payload = directory / "payload"
                payload.mkdir(parents=True)
                (directory / "orthohmm_working_res").mkdir()
                for item in previous["graph_inputs"]:
                    source = Path(item["path"])
                    name = source.name.removeprefix("rbnh_")
                    (payload / name).symlink_to(source)
                (payload / "metadata.json").write_text(json.dumps({"cpm_resolution": .1, "seed": 4,
                    "include_isolates": True, "output_directory": str(directory)}) + "\n")
                command = [sys.executable, str(Path(__file__).resolve()), "--root", str(root), "--output", str(output),
                           "--worker-payload", str(payload), "--mode", mode]
                with (directory / "worker.log").open("x") as log:
                    run = subprocess.run(command, cwd=launcher, env=env, stdout=log, stderr=subprocess.STDOUT)
                if run.returncode:
                    raise RuntimeError("Construction worker failed: " + str(directory))
                result = json.loads((payload / "construction.json").read_text())
                report["workers"].append({"mode": mode, "index": index, "result": result,
                    "snapshot": record(payload / "worker_before.json"), "observation": record(payload / "construction.json")})
                (output / "progress.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        if verify(frozen, launcher, runtime) != before:
            raise ValueError("Frozen runtime changed")
        for item in previous["graph_inputs"]:
            check(item)
        report["status"] = "six_construction_workers_complete_unscored"
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        (output / "results.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
