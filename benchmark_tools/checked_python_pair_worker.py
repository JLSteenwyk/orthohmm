"""Run one frozen saved-graph worker with checked Python-pair construction."""

import argparse
from contextlib import contextmanager
import hashlib
import json
from pathlib import Path
import sys


def require_admission(report):
    expected = [(i, fmt) for i in range(3) for fmt in ("numpy", "python_pairs")]
    if (report["status"] != "constructor_format_observations_verified"
            or report["accuracy_evaluated"] is not False or report["optimizer_called"] is not False
            or [(row["index"], row["edge_format"]) for row in report["workers"]] != expected
            or any(row["mode"] != "minimal_imports" for row in report["workers"])):
        raise ValueError("Require complete independently admitted constructor-format panel")
    for row in report["workers"]:
        if row["edge_format"] == "python_pairs" and (row["before_mismatched_edges"] != 0 or row["after_mismatched_edges"] != 0):
            raise ValueError("Python-pair construction did not preserve the saved graph in every admitted repeat")


@contextmanager
def python_pair_constructor(module, destination):
    import numpy as np
    original = module.Graph.__init__
    calls = []

    def flush():
        destination.write_text(json.dumps({"format": "python_pairs", "calls": calls,
            "accuracy_evaluated": False}, indent=2, sort_keys=True) + "\n")

    def construct(graph, *args, **kwargs):
        edges = kwargs.get("edges")
        if not isinstance(edges, np.ndarray):
            return original(graph, *args, **kwargs)
        if calls or args or edges.ndim != 2 or edges.shape[1] != 2 or edges.dtype != np.int32:
            raise ValueError("Expected one keyword-only int32 edge-array constructor")
        digest = hashlib.sha256(edges.tobytes(order="C")).hexdigest()
        row = {"status": "before_constructor", "shape": list(edges.shape), "dtype": str(edges.dtype),
               "ordered_input_bytes_sha256": digest, "n": kwargs.get("n"), "directed": kwargs.get("directed")}
        calls.append(row)
        flush()
        pairs = ((int(a), int(b)) for a, b in edges)
        original(graph, **{**kwargs, "edges": pairs})
        if hashlib.sha256(edges.tobytes(order="C")).hexdigest() != digest:
            raise ValueError("Constructor input changed")
        row["status"] = "constructor_returned"
        flush()

    module.Graph.__init__ = construct
    try:
        yield calls
    finally:
        module.Graph.__init__ = original


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--payload", required=True, type=Path)
    parser.add_argument("--admission", required=True, type=Path)
    parser.add_argument("--admission-sha256", required=True)
    args = parser.parse_args()
    root, payload = args.root.resolve(), args.payload.resolve()
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from benchmark_tools.run_simulation_methods import read_frozen
    from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
    from repeat_qfo_saved_graph import worker, set_worker_affinity
    import os
    set_worker_affinity([min(os.sched_getaffinity(0))])
    admission = read_frozen(args.admission, args.admission_sha256)
    require_admission(admission)
    for item in admission["provenance_checked"]:
        check(item)
    observed = [record(payload / name) for name in ("gene_names.txt", "sources.npy", "targets.npy", "weights.npy")]
    if observed != admission["native_report"]["graph_inputs"]:
        raise ValueError("Worker payload differs from admitted saved graph")
    if (payload / "native_boundary.json").exists() or (payload / "constructor_adapter.json").exists():
        raise FileExistsError("Refuse to reuse a worker payload")
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    sys.path.insert(0, str(launcher))
    import igraph
    (payload / "checked_worker_provenance.json").write_text(json.dumps({"source": record(__file__),
        "admission": record(args.admission), "inputs": observed,
        "helpers": [record(Path(__file__).with_name(name)) for name in ("repeat_qfo_saved_graph.py", "probe_leiden_boundary.py")],
        "scope": "One initial-graph diagnostic; complete native graph gate before and after unchanged optimizer; no accuracy selection"},
        indent=2, sort_keys=True) + "\n")
    with python_pair_constructor(igraph, payload / "constructor_adapter.json"):
        # The frozen worker exits the process only after the observer records the optimizer result.
        worker(launcher, payload, native_boundary=True)
    raise RuntimeError("Frozen worker unexpectedly returned")


if __name__ == "__main__":
    main()
