from types import SimpleNamespace

import numpy as np
import pytest

from orthohmm.helpers import IndexedEdges
from benchmark_tools.capture_qfo_replay_graph import CaptureComplete, Recorder, fingerprint


def edges():
    return IndexedEdges(["a", "b"], np.array([0], dtype=np.int32), np.array([1], dtype=np.int32), np.array([1.]))


def test_capture_stops_before_second_native_cluster_and_preserves_arrays(tmp_path):
    calls = []
    working = tmp_path / "replay/orthohmm_working_res"
    working.mkdir(parents=True)
    graph = edges()

    def cluster(*args, **kwargs):
        calls.append("cluster")
        assert kwargs["edges"] is graph
        (working / "orthohmm_edges_clustered.txt").write_text("a b\n")

    def singleton(*args, **kwargs):
        calls.append("singleton")
        return graph

    recorder = Recorder(SimpleNamespace(execute_leiden=cluster, build_singleton_assignment_edges=singleton), tmp_path)
    recorder.cluster_call(.1, tmp_path / "replay", graph, True, 4)
    assert recorder.singleton_call("unchanged") is graph
    with pytest.raises(CaptureComplete):
        recorder.cluster_call(.1, tmp_path / "replay", graph, True, 4)
    assert calls == ["cluster", "singleton"]
    assert recorder.report["status"] == "captured_before_second_clustering"
    assert (tmp_path / "initial_partition.txt").read_text() == "a b\n"
    for key in ("sources", "targets", "weights"):
        assert fingerprint(np.load(tmp_path / ("rbnh_" + key + ".npy"))) == recorder.report["rbnh_arrays"][key]
    assert recorder.report["singleton_edges"] == 1


@pytest.mark.parametrize("resolution,isolate,seed", [(.2, True, 4), (.1, False, 4), (.1, True, 0)])
def test_parameter_drift_rejected_before_native_call(tmp_path, resolution, isolate, seed):
    def forbidden(*args, **kwargs):
        pytest.fail("Unexpected native call")
    recorder = Recorder(SimpleNamespace(execute_leiden=forbidden, build_singleton_assignment_edges=forbidden), tmp_path)
    with pytest.raises(ValueError):
        recorder.cluster_call(resolution, tmp_path, edges(), isolate, seed)


def test_out_of_order_singleton_call_rejected(tmp_path):
    module = SimpleNamespace(execute_leiden=lambda *a: None, build_singleton_assignment_edges=lambda *a: edges())
    with pytest.raises(ValueError):
        Recorder(module, tmp_path).singleton_call()
