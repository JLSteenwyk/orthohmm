"""Observe a frozen worker's graph at optimizer entry, then stop without clustering."""

from contextlib import contextmanager
import inspect
import json
import os

from benchmark_tools.probe_leiden_boundary import endpoint_differences, graph_fingerprint, saved_fingerprint


class DiagnosticStop(Exception):
    """Expected termination before any optimizer invocation."""


@contextmanager
def stop_before_partition(module, payload):
    original = module.find_partition
    destination = payload / "preoptimizer_stop.json"
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    calls = []

    def observe(graph, *args, **kwargs):
        if calls:
            raise ValueError("Repeated optimizer-entry observation")
        calls.append(True)
        frame = inspect.currentframe().f_back
        try:
            constructor = frame.f_locals.get("graph_edges")
            if constructor is None:
                raise ValueError("Missing frozen worker constructor array")
            differences = endpoint_differences(graph, payload, constructor)
            fingerprint = graph_fingerprint(graph)
            saved = saved_fingerprint(payload)
            if fingerprint != saved or any(differences[name]["different_edges"] for name in (
                    "native_vs_saved", "constructor_vs_saved", "native_vs_constructor")):
                raise ValueError("Graph differs at diagnostic optimizer boundary")
            result = dict(status="stopped_before_optimizer", optimizer_called=False,
                accuracy_evaluated=False, publication_ready=False, fingerprint=fingerprint,
                saved=saved, differences=differences,
                limitation="Observed frozen-worker entry only; not a historical crash reproduction or clustering result.")
            with destination.open("x") as stream:
                json.dump(result, stream, indent=2, sort_keys=True)
                stream.write("\n")
                stream.flush()
                os.fsync(stream.fileno())
        finally:
            del frame
        raise DiagnosticStop("Graph checked; optimizer intentionally not invoked")

    module.find_partition = observe
    try:
        yield calls
    finally:
        module.find_partition = original
