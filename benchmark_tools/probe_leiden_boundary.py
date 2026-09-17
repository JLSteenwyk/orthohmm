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
