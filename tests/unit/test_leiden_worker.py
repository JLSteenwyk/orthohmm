"""Validate the short-lived worker's payload and no-pickle boundary."""

import json

import numpy as np
import pytest

from orthohmm import leiden_worker


@pytest.fixture
def payload(tmp_path):
    metadata = dict(cpm_resolution=.1, output_directory=str(tmp_path / "output"),
                    include_isolates=True, seed=4)
    (tmp_path / "metadata.json").write_text(json.dumps(metadata))
    (tmp_path / "gene_names.txt").write_text("second\nfirst\nisolate\n")
    for name, data in (("sources", np.array([0], dtype=np.int32)),
                       ("targets", np.array([1], dtype=np.int32)),
                       ("weights", np.array([1.], dtype=np.float64))):
        np.save(tmp_path / (name + ".npy"), data)
    return tmp_path, metadata


def test_worker_preserves_order_options_and_read_only_payload(payload, monkeypatch):
    directory, metadata = payload
    calls = []
    monkeypatch.setattr(leiden_worker, "_execute_leiden_in_process", lambda *a, **k: calls.append((a, k)))

    def exit_worker(code):
        raise SystemExit(code)

    monkeypatch.setattr(leiden_worker.os, "_exit", exit_worker)
    with pytest.raises(SystemExit) as error:
        leiden_worker.main([str(directory)])
    assert error.value.code == 0 and len(calls) == 1
    arguments, options = calls[0]
    assert arguments == (metadata["cpm_resolution"], metadata["output_directory"])
    assert options["include_isolates"] is True and options["seed"] == 4
    edges = options["edges"]
    assert edges.gene_names == ["second", "first", "isolate"]
    for actual, expected in ((edges.sources, [0]), (edges.targets, [1]), (edges.weights, [1.])):
        np.testing.assert_array_equal(actual, expected)
        assert isinstance(actual, np.memmap) and not actual.flags.writeable


@pytest.mark.parametrize("arguments", [[], ["one", "two"]])
def test_argument_errors_do_not_attempt_payload_reads(arguments, monkeypatch):
    monkeypatch.setattr(leiden_worker.np, "load", lambda *a, **k: pytest.fail("must validate arguments first"))
    with pytest.raises(SystemExit, match="usage:"):
        leiden_worker.main(arguments)


def test_missing_payload_is_not_reported_as_success(tmp_path, monkeypatch):
    monkeypatch.setattr(leiden_worker, "_execute_leiden_in_process", lambda *a, **k: pytest.fail("must not cluster"))
    with pytest.raises(FileNotFoundError):
        leiden_worker.main([str(tmp_path)])


@pytest.mark.parametrize("name", ["sources", "targets", "weights"])
def test_object_arrays_cannot_enable_pickle_loading(payload, monkeypatch, name):
    directory, _ = payload
    np.save(directory / (name + ".npy"), np.array(["object"], dtype=object))
    monkeypatch.setattr(leiden_worker, "_execute_leiden_in_process", lambda *a, **k: pytest.fail("must not cluster"))
    with pytest.raises(ValueError, match="Python objects|Object arrays|object"):
        leiden_worker.main([str(directory)])
