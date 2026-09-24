"""Require exact constructor/graph parity before a single optimizer call."""

from contextlib import contextmanager
import inspect
import json
import os
from pathlib import Path

from benchmark_tools.probe_leiden_boundary import endpoint_differences, graph_fingerprint, saved_fingerprint


@contextmanager
def require_exact_constructor(module, payload):
    """Wrap the existing boundary observer without altering optimizer arguments."""
    original = module.find_partition
    destination = payload / "constructor_parity.json"
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    called = False

    def guarded(graph, *args, **kwargs):
        nonlocal called
        if called:
            raise ValueError("Repeated optimizer invocation")
        called = True
        frame = inspect.currentframe().f_back
        try:
            constructor = frame.f_locals.get("graph_edges")
            if constructor is None:
                raise ValueError("Missing frozen caller constructor array")
            differences = endpoint_differences(graph, payload, constructor)
            observed, saved = graph_fingerprint(graph), saved_fingerprint(payload)
            if observed != saved or any(differences[name]["different_edges"] for name in (
                    "native_vs_saved", "constructor_vs_saved", "native_vs_constructor")):
                raise ValueError("Recovery graph or constructor differs from saved input")
            report = dict(status="constructor_parity_verified_before_optimizer", differences=differences,
                          fingerprint=observed, saved=saved, accuracy_evaluated=False,
                          limitation="Pre-call check only; optimizer completion requires separate boundary evidence.")
            with destination.open("x") as stream:
                json.dump(report, stream, indent=2, sort_keys=True)
                stream.write("\n")
                stream.flush()
                os.fsync(stream.fileno())
        finally:
            del frame
        return original(graph, *args, **kwargs)

    module.find_partition = guarded
    try:
        yield
    finally:
        module.find_partition = original


def run_frozen_worker(launcher, payload):
    """Internal child entry; the recovery parent must first validate admission."""
    from benchmark_tools.checked_python_pair_worker import python_pair_constructor
    from benchmark_tools.probe_leiden_boundary import observe_partition
    from benchmark_tools.repeat_qfo_saved_graph import worker, validate_clustering_metadata

    if Path.cwd() != launcher:
        raise ValueError("Wrong recovery worker working directory")
    metadata = json.loads((payload / "metadata.json").read_text())
    validate_clustering_metadata(metadata, .12)
    if metadata["output_directory"] != str(payload.parent):
        raise ValueError("Recovery output is outside its fresh directory")
    for name in ("constructor_parity.json", "constructor_adapter.json", "native_boundary.json", "worker_before.json"):
        destination = payload / name
        if destination.exists() or destination.is_symlink():
            raise FileExistsError(destination)
    (payload.parent / "orthohmm_working_res").mkdir(exist_ok=False)
    import igraph
    import leidenalg
    with python_pair_constructor(igraph, payload / "constructor_adapter.json"):
        with observe_partition(leidenalg, payload):
            with require_exact_constructor(leidenalg, payload):
                worker(launcher, payload, requested_affinity=[min(os.sched_getaffinity(0))],
                       native_boundary=False, expected_cpm_resolution=.12)
