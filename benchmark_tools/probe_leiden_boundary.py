"""Observe native graph contents without changing optimizer arguments."""

from contextlib import contextmanager
import hashlib
import inspect
import json


def graph_fingerprint(graph, chunk_size=100000):
    import numpy as np

    if chunk_size < 1:
        raise ValueError("Chunk size must be positive")
    endpoints, weights = hashlib.sha256(), hashlib.sha256()
    for start in range(0, graph.ecount(), chunk_size):
        edges = graph.es[start:start + chunk_size]
        pairs = np.asarray([edge.tuple for edge in edges], dtype="<i8").reshape(-1, 2)
        if not graph.is_directed():
            pairs.sort(axis=1)
        values = np.asarray(edges["weight"], dtype="<f8")
        if not np.isfinite(values).all():
            raise ValueError("Nonfinite native weights")
        endpoints.update(pairs.tobytes())
        weights.update(values.tobytes())
    return {"vertices": graph.vcount(), "edges": graph.ecount(),
            "directed": graph.is_directed(), "ordered_endpoints_sha256": endpoints.hexdigest(),
            "ordered_weights_sha256": weights.hexdigest()}


def saved_fingerprint(payload, chunk_size=100000):
    import numpy as np

    if chunk_size < 1:
        raise ValueError("Chunk size must be positive")
    arrays = [np.load(payload / (name + ".npy"), mmap_mode="r", allow_pickle=False)
              for name in ("sources", "targets", "weights")]
    if any(array.ndim != 1 or len(array) != len(arrays[0]) for array in arrays):
        raise ValueError("Invalid saved graph shapes")
    endpoints, weights = hashlib.sha256(), hashlib.sha256()
    for start in range(0, len(arrays[0]), chunk_size):
        pairs = np.column_stack([array[start:start + chunk_size] for array in arrays[:2]]).astype("<i8")
        pairs.sort(axis=1)
        values = np.asarray(arrays[2][start:start + chunk_size], dtype="<f8")
        if not np.isfinite(values).all():
            raise ValueError("Nonfinite saved weights")
        endpoints.update(pairs.tobytes())
        weights.update(values.tobytes())
    with (payload / "gene_names.txt").open() as handle:
        vertices = sum(1 for _ in handle)
    return {"vertices": vertices, "edges": len(arrays[0]), "directed": False,
            "ordered_endpoints_sha256": endpoints.hexdigest(),
            "ordered_weights_sha256": weights.hexdigest()}


def endpoint_differences(graph, payload, constructor_edges=None, chunk_size=100000):
    """Locate endpoint changes on either side of the igraph constructor."""
    import numpy as np

    sources, targets = [np.load(payload / (name + ".npy"), mmap_mode="r", allow_pickle=False)
                        for name in ("sources", "targets")]
    if len(sources) != graph.ecount() or len(targets) != graph.ecount():
        raise ValueError("Cannot compare differently sized endpoint arrays")
    if constructor_edges is not None and constructor_edges.shape != (graph.ecount(), 2):
        raise ValueError("Unexpected constructor array shape")
    report = {"native_vs_saved": {"different_edges": 0, "examples": []}}
    if constructor_edges is not None:
        report.update(constructor_dtype=str(constructor_edges.dtype),
                      constructor_c_contiguous=bool(constructor_edges.flags.c_contiguous),
                      constructor_vs_saved={"different_edges": 0, "examples": []},
                      native_vs_constructor={"different_edges": 0, "examples": []})
    def compare(label, left, right, start):
        different = np.flatnonzero(np.any(left != right, axis=1))
        result = report[label]
        result["different_edges"] += len(different)
        for index in different[:max(0, 20 - len(result["examples"]))]:
            result["examples"].append({"edge_index": start + int(index),
                                       "left": left[index].tolist(), "right": right[index].tolist()})
    for start in range(0, graph.ecount(), chunk_size):
        expected = np.column_stack((sources[start:start + chunk_size], targets[start:start + chunk_size])).astype("<i8")
        observed = np.asarray([edge.tuple for edge in graph.es[start:start + chunk_size]], dtype="<i8").reshape(-1, 2)
        expected.sort(axis=1)
        observed.sort(axis=1)
        compare("native_vs_saved", observed, expected, start)
        if constructor_edges is not None:
            constructed = np.array(constructor_edges[start:start + chunk_size], dtype="<i8", copy=True)
            constructed.sort(axis=1)
            compare("constructor_vs_saved", constructed, expected, start)
            compare("native_vs_constructor", observed, constructed, start)
    return report


@contextmanager
def observe_partition(module, payload):
    original = module.find_partition
    calls = []
    destination = payload / "native_boundary.json"

    def flush():
        destination.write_text(json.dumps({"calls": calls, "accuracy_evaluated": False},
                                          indent=2, sort_keys=True) + "\n")

    def observed(*args, **kwargs):
        bound = inspect.signature(original).bind(*args, **kwargs)
        bound.apply_defaults()
        arguments = dict(bound.arguments)
        graph = arguments.pop("graph")
        partition_type = arguments.pop("partition_type")
        arguments["partition_type"] = partition_type.__module__ + "." + partition_type.__qualname__
        row = {"arguments": arguments, "before": graph_fingerprint(graph),
               "saved": saved_fingerprint(payload), "status": "before_optimizer"}
        calls.append(row)
        flush()
        if row["before"] != row["saved"]:
            frame = inspect.currentframe().f_back
            try:
                # Observe the frozen caller's actual constructor input without changing it.
                constructor = frame.f_locals.get("graph_edges")
                row["caller"] = {"function": frame.f_code.co_name, "file": frame.f_code.co_filename}
                row["endpoint_differences"] = endpoint_differences(graph, payload, constructor)
            finally:
                del frame
            row["status"] = "native_graph_mismatch_before_optimizer"
            flush()
            raise ValueError("Native graph differs from ordered saved graph")
        try:
            result = original(*args, **kwargs)
        except BaseException as error:
            row.update(status="optimizer_failed", error_type=type(error).__name__)
            flush()
            raise
        row.update(after=graph_fingerprint(graph), status="optimizer_returned")
        flush()
        if row["after"] != row["before"]:
            raise ValueError("Optimizer changed graph contents")
        return result

    module.find_partition = observed
    try:
        yield
    finally:
        module.find_partition = original
