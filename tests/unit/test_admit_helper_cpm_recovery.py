import copy
import hashlib
import io
import json
from pathlib import Path
import struct

import numpy as np
import pytest

from benchmark_tools import admit_helper_cpm_recovery as admission
from benchmark_tools.probe_leiden_boundary import saved_fingerprint


def graph_fixture(tmp_path):
    payload = tmp_path / "payload"
    payload.mkdir()
    (payload / "gene_names.txt").write_text("a\nb\nc\nd\ne\n")
    arrays = [np.array([3, 0, 2, 1], dtype="<i4"), np.array([1, 2, 2, 4], dtype="<i4"),
              np.array([.1, .2, .3, .4], dtype="<f8")]
    for name, data in zip(("sources", "targets", "weights"), arrays):
        np.save(payload / (name + ".npy"), data, allow_pickle=False)
    return payload, arrays


@pytest.mark.parametrize("chunk_size", [1, 2, 10])
def test_stdlib_graph_fingerprint_matches_existing_numpy_reader(tmp_path, chunk_size):
    payload, arrays = graph_fixture(tmp_path)
    graph, constructor = admission.stream_graph(payload, 5, chunk_size)
    assert graph == saved_fingerprint(payload, chunk_size=2)
    expected = np.column_stack(arrays[:2]).astype("<i4").tobytes()
    assert constructor == hashlib.sha256(expected).hexdigest()


@pytest.mark.parametrize("bad", ["negative", "out_of_range", "nan", "inf", "dtype", "shape",
                               "length", "truncated", "trailing", "chunk", "vertices"])
def test_invalid_graph_cannot_be_admitted(tmp_path, bad):
    payload, arrays = graph_fixture(tmp_path)
    if bad == "negative": arrays[0][0] = -1
    elif bad == "out_of_range": arrays[1][0] = 5
    elif bad == "nan": arrays[2][0] = np.nan
    elif bad == "inf": arrays[2][0] = np.inf
    elif bad == "dtype": arrays[0] = arrays[0].astype("<i8")
    elif bad == "shape": arrays[0] = arrays[0].reshape(2, 2)
    elif bad == "length": arrays[2] = arrays[2][:-1]
    for name, data in zip(("sources", "targets", "weights"), arrays):
        np.save(payload / (name + ".npy"), data, allow_pickle=False)
    if bad == "truncated":
        path = payload / "weights.npy"
        path.write_bytes(path.read_bytes()[:-1])
    elif bad == "trailing":
        path = payload / "weights.npy"
        path.write_bytes(path.read_bytes() + b"extra")
    with pytest.raises(ValueError):
        admission.stream_graph(payload, 0 if bad == "vertices" else 5, 0 if bad == "chunk" else 2)


@pytest.mark.parametrize("bad", ["magic", "version", "length", "header", "large", "fortran", "bool_shape"])
def test_invalid_npy_header_rejected(bad):
    header = dict(descr="<i4", fortran_order=bad == "fortran", shape=(True if bad == "bool_shape" else 4,))
    raw = repr(header).encode()
    data = b"\x93NUMPY\x01\x00" + struct.pack("<H", len(raw)) + raw
    if bad == "magic": data = b"wrong!" + data[6:]
    elif bad == "version": data = data[:6] + b"\x03\x00" + data[8:]
    elif bad == "length": data = data[:9]
    elif bad == "header": data = data[:-1]
    elif bad == "large": data = b"\x93NUMPY\x02\x00" + struct.pack("<I", 65537)
    with pytest.raises(ValueError):
        admission.npy_header(io.BytesIO(data), "<i4")


def test_version_two_header_accepted():
    raw = repr(dict(descr="<i4", fortran_order=False, shape=(4,))).encode()
    data = b"\x93NUMPY\x02\x00" + struct.pack("<I", len(raw)) + raw
    assert admission.npy_header(io.BytesIO(data), "<i4") == 4


def optimizer_fixture(tmp_path, monkeypatch):
    payload, arrays = graph_fixture(tmp_path)
    graph, constructor = admission.stream_graph(payload, 5)
    monkeypatch.setattr(admission, "CONSTRUCTOR_SHA", constructor)
    python = tmp_path / "python"
    python.write_bytes(b"python")
    monkeypatch.setattr(admission, "HISTORICAL_PYTHON", str(python))
    launcher = tmp_path / "source"
    (launcher / "orthohmm").mkdir(parents=True)
    modules = {}
    for name in ("orthohmm.leiden_worker", "orthohmm.externals", "orthohmm.helpers"):
        path = launcher / (name.replace(".", "/") + ".py")
        path.write_bytes(b"pass\n")
        modules[name] = admission.record(path)
    snapshot = dict(status="before_native_clustering", accuracy_evaluated=False,
        metadata=dict(cpm_resolution=.12, seed=4, include_isolates=True, output_directory=str(payload.parent)),
        cwd=str(launcher), cpu_affinity=[2], native_libraries=[dict(path="/native")],
        environment=dict(PYTHONPATH=str(launcher), PYTHONHASHSEED="0", OMP_NUM_THREADS="1",
                         OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1"),
        inputs=[admission.record(payload / name) for name in
                ("gene_names.txt", "sources.npy", "targets.npy", "weights.npy")],
        modules=modules, python=admission.record(python))
    arguments = dict(initial_membership=None, weights="weight", n_iterations=2, max_comm_size=0,
        seed=4, kwargs={"resolution_parameter": .12}, partition_type="leidenalg.VertexPartition.CPMVertexPartition")
    boundary = dict(accuracy_evaluated=False, calls=[dict(arguments=arguments, before=graph, saved=graph,
                   after=graph, status="optimizer_returned")])
    adapter = dict(format="python_pairs", accuracy_evaluated=False, calls=[dict(status="constructor_returned",
        shape=[4, 2], dtype="int32", n=5, directed=False, ordered_input_bytes_sha256=constructor)])
    parity = dict(status="constructor_parity_verified_before_optimizer", accuracy_evaluated=False,
        fingerprint=graph, saved=graph, differences=dict(constructor_dtype="int32", constructor_c_contiguous=True,
        **{key: dict(different_edges=0, examples=[]) for key in
           ("native_vs_saved", "constructor_vs_saved", "native_vs_constructor")}))
    return snapshot, boundary, adapter, parity, graph, payload, launcher


def test_exact_original_optimizer_contract(tmp_path, monkeypatch):
    admission.optimizer_contract(*optimizer_fixture(tmp_path, monkeypatch))


@pytest.mark.parametrize("bad", ["status", "metadata", "environment", "affinity", "libraries", "inputs",
                               "python", "source", "extra_call", "failed_call", "arguments", "after_graph",
                               "adapter", "parity", "different_edges"])
def test_optimizer_contract_rejects_drift_or_partial_evidence(tmp_path, monkeypatch, bad):
    values = optimizer_fixture(tmp_path, monkeypatch)
    snapshot, boundary, adapter, parity = values[:4]
    if bad == "status": snapshot["status"] = "failed"
    elif bad == "metadata": snapshot["metadata"]["cpm_resolution"] = .1
    elif bad == "environment": snapshot["environment"]["OMP_NUM_THREADS"] = "32"
    elif bad == "affinity": snapshot["cpu_affinity"] = [1, 2]
    elif bad == "libraries": snapshot["native_libraries"] = []
    elif bad == "inputs": snapshot["inputs"].pop()
    elif bad == "python": snapshot["python"]["sha256"] = "wrong"
    elif bad == "source": snapshot["modules"]["orthohmm.externals"]["sha256"] = "wrong"
    elif bad == "extra_call": boundary["calls"].append(copy.deepcopy(boundary["calls"][0]))
    elif bad == "failed_call": boundary["calls"][0]["status"] = "optimizer_failed"
    elif bad == "arguments": boundary["calls"][0]["arguments"]["seed"] = 7
    elif bad == "after_graph": boundary["calls"][0]["after"] = dict(values[4], edges=3)
    elif bad == "adapter": adapter["calls"][0]["dtype"] = "int64"
    elif bad == "parity": parity["status"] = "failed"
    else: parity["differences"]["native_vs_saved"]["different_edges"] = 1
    with pytest.raises(ValueError):
        admission.optimizer_contract(*values)


@pytest.mark.parametrize("text,groups", [("a b\nc\n", 2), ("a\nb c\n", 2)])
def test_complete_partition_membership(tmp_path, text, groups):
    path = tmp_path / "partition"
    path.write_text(text)
    assert admission.partition_coverage(path, dict(a=0, b=1, c=2), groups) == dict(genes=3, groups=2)


@pytest.mark.parametrize("text,groups", [("a b\n", 2), ("a b c\n", 2), ("a b\nb c\n", 2),
                                       ("a a\nb c\n", 2), ("a b\nunknown\n", 2)])
def test_incomplete_duplicate_unknown_or_wrong_count_partition(tmp_path, text, groups):
    path = tmp_path / "partition"
    path.write_text(text)
    with pytest.raises(ValueError):
        admission.partition_coverage(path, dict(a=0, b=1, c=2), groups)


@pytest.mark.parametrize("bad", [None, "completed", "cpu", "duplicate", "missing"])
def test_historical_failures_must_remain_failed(bad):
    row = "22155|FAILED|1:0|2|bizon\n"
    if bad == "completed": row = row.replace("FAILED", "COMPLETED")
    elif bad == "cpu": row = row.replace("|2|", "|1|")
    elif bad == "duplicate": row *= 2
    elif bad == "missing": row = ""
    text = "JobID|State|ExitCode|AllocCPUS|NodeList\n" + row
    if bad:
        with pytest.raises(ValueError): admission.scheduler_failure(text, "22155", 2)
    else: assert admission.scheduler_failure(text, "22155", 2)["State"] == "FAILED"


def test_exact_historical_contracts_load_without_import_side_effects(monkeypatch):
    path = Path(admission.__file__).with_name("admit_cpm_checkpoint_recovery.py")
    monkeypatch.setattr(admission, "CONTRACT_SHA", admission.record(path)["sha256"])
    completed, parent = admission.original_contracts(path)
    row = "123|COMPLETED|0:0|1|64G|bizon\n"
    text = "JobID|State|ExitCode|AllocCPUS|ReqMem|NodeList\n" + row
    assert completed(text, "123")["State"] == "COMPLETED"
    assert callable(parent)
    with pytest.raises(ValueError): completed(text.replace("COMPLETED", "FAILED"), "123")


def test_changed_historical_contract_is_rejected(tmp_path):
    path = tmp_path / "contract.py"
    path.write_text("raise RuntimeError('must never execute')\n")
    with pytest.raises(ValueError, match="Changed"):
        admission.original_contracts(path)


def test_failed_preflight_retained_without_native_attempt(tmp_path, monkeypatch):
    def fail(*_): raise ValueError("changed evidence")
    monkeypatch.setattr(admission, "_admit", fail)
    output = tmp_path / "output"
    with pytest.raises(ValueError): admission.admit(tmp_path, output, "protocol")
    report = json.loads((output / "status.json").read_bytes())
    assert report["status"] == "helper_recovery_preflight_failed"
    assert report["seed_admitted"] is False and report["native_attempts"] == 0
    before = (output / "status.json").read_bytes()
    with pytest.raises(FileExistsError): admission.admit(tmp_path, output, "protocol")
    assert (output / "status.json").read_bytes() == before
