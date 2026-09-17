import json
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pytest

from benchmark_tools.checked_replay_payload_worker import FILES, validate_payload
from benchmark_tools.checked_replay_interceptor import CheckedReplaySubprocess
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def payload_fixture(tmp_path):
    payload = tmp_path / "original"
    payload.mkdir()
    (payload / "gene_names.txt").write_text("a\nb\nc\n")
    for name, array in (("sources", np.array([0], dtype=np.int32)), ("targets", np.array([1], dtype=np.int32)),
                        ("weights", np.array([1.]))):
        np.save(payload / (name + ".npy"), array)
    output = tmp_path / "replay"
    (output / "orthohmm_working_res").mkdir(parents=True)
    (payload / "metadata.json").write_text(json.dumps({"seed": 4, "cpm_resolution": .1,
        "include_isolates": True, "output_directory": str(output)}))
    manifest = {"stage": "initial", "index": 0, "accuracy_evaluated": False,
                "inputs": [record(payload / name) for name in FILES], "output_directory": str(output)}
    return payload, output, manifest


@pytest.mark.parametrize("problem", [None, "index", "names", "settings", "observed", "input"])
def test_payload_gate(tmp_path, problem):
    payload, output, manifest = payload_fixture(tmp_path)
    admitted_names = dict(manifest["inputs"][0])
    if problem == "index":
        manifest["index"] = 1
    elif problem == "names":
        admitted_names["sha256"] = "wrong"
    elif problem == "settings":
        manifest["output_directory"] = str(tmp_path / "elsewhere")
    elif problem == "observed":
        (payload / "worker_before.json").write_text("{}")
    elif problem == "input":
        (payload / "gene_names.txt").write_text("a\nb\nd\n")
    if problem:
        with pytest.raises((ValueError, FileExistsError)):
            validate_payload(manifest, payload, admitted_names)
    else:
        assert validate_payload(manifest, payload, admitted_names)["output_directory"] == str(output)


@pytest.mark.parametrize("fail", [False, True])
def test_interceptor_preserves_payload_and_failure_without_retry(tmp_path, fail):
    payload, output, _ = payload_fixture(tmp_path)
    calls = []
    def run(command, **kwargs):
        calls.append(command)
        if fail:
            raise subprocess.CalledProcessError(1, command)
        (output / "orthohmm_working_res/orthohmm_edges_clustered.txt").write_text("a b\nc\n")
        return SimpleNamespace(returncode=0)
    original = SimpleNamespace(run=run, STDOUT=subprocess.STDOUT)
    def validate(retained, manifest):
        assert retained != payload
        assert manifest["inputs"] == [record(retained / name) for name in FILES]
        return {"status": "fixture_checked"}
    proxy = CheckedReplaySubprocess(original, tmp_path, tmp_path / "observations", tmp_path / "worker.py", validate)
    command = [sys.executable, "-m", "orthohmm.leiden_worker", str(payload)]
    if fail:
        with pytest.raises(subprocess.CalledProcessError):
            proxy.run(command, check=True)
        with pytest.raises(ValueError):
            proxy.run(command, check=True)
        assert proxy.calls[0]["status"] == "failed"
    else:
        assert proxy.run(command, check=True).returncode == 0
        assert Path(proxy.calls[0]["partition"]["path"]).read_text() == "a b\nc\n"
        assert proxy.calls[0]["status"] == "checked"
    assert len(calls) == 1
    directory = tmp_path / "observations/cluster_0_initial"
    assert (directory / "execution.json").exists()
    assert all((directory / "payload" / name).exists() for name in FILES)


def test_unrelated_subprocess_passes_through(tmp_path):
    calls = []
    original = SimpleNamespace(run=lambda *a, **k: calls.append((a, k)))
    proxy = CheckedReplaySubprocess(original, tmp_path, tmp_path / "unused", tmp_path / "worker.py", None)
    proxy.run(["other", "--flag"], cwd=tmp_path)
    assert calls == [((["other", "--flag"],), {"cwd": tmp_path})]
    assert proxy.calls == []


def test_exact_four_stage_inventory(tmp_path):
    payload, output, _ = payload_fixture(tmp_path)
    def run(*args, **kwargs):
        (output / "orthohmm_working_res/orthohmm_edges_clustered.txt").write_text("a b\nc\n")
        return SimpleNamespace(returncode=0)
    proxy = CheckedReplaySubprocess(SimpleNamespace(run=run, STDOUT=subprocess.STDOUT), tmp_path,
        tmp_path / "observations", tmp_path / "worker.py", lambda p, m: {"stage": m["stage"]})
    command = [sys.executable, "-m", "orthohmm.leiden_worker", str(payload)]
    for _ in range(4):
        proxy.run(command, check=True)
    assert [r["stage"] for r in proxy.calls] == ["initial", "multipass", "profile_base", "profile_expanded"]
    with pytest.raises(ValueError):
        proxy.run(command, check=True)
    assert len(proxy.calls) == 4


def test_profile_parent_threads_do_not_leak_into_clustering(tmp_path, monkeypatch):
    payload, output, _ = payload_fixture(tmp_path)
    for key, value in (("OMP_NUM_THREADS", "32"), ("OPENBLAS_NUM_THREADS", "8"), ("MKL_NUM_THREADS", "4"),
                       ("PYTHONPATH", "/unchanged/launcher"), ("PRESERVE_TEST_VARIABLE", "yes")):
        monkeypatch.setenv(key, value)
    parent = dict(os.environ)
    children = []
    def run(command, **kwargs):
        children.append(kwargs["env"])
        assert kwargs["env"] is not os.environ
        for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
            assert kwargs["env"][key] == "1"
        assert kwargs["env"]["PYTHONPATH"] == "/unchanged/launcher"
        assert kwargs["env"]["PRESERVE_TEST_VARIABLE"] == "yes"
        (output / "orthohmm_working_res/orthohmm_edges_clustered.txt").write_text("a b\nc\n")
        return SimpleNamespace(returncode=0)
    proxy = CheckedReplaySubprocess(SimpleNamespace(run=run, STDOUT=subprocess.STDOUT), tmp_path,
        tmp_path / "observations", tmp_path / "worker.py", lambda p, m: {})
    for _ in range(4):
        proxy.run([sys.executable, "-m", "orthohmm.leiden_worker", str(payload)], check=True)
    assert dict(os.environ) == parent
    assert len(children) == 4
    assert all(row["thread_environment"]["inherited"]["OMP_NUM_THREADS"] == "32" for row in proxy.calls)
