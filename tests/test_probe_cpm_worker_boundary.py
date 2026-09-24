import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

import benchmark_tools.probe_cpm_worker_boundary as module


def result():
    return dict(status="stopped_before_optimizer", optimizer_called=False, accuracy_evaluated=False,
        fingerprint={"fixture": True}, saved={"fixture": True}, differences={
            k: {"different_edges": 0} for k in ("native_vs_saved", "constructor_vs_saved", "native_vs_constructor")})


@pytest.mark.parametrize("change", [None, "optimizer", "fingerprint", "constructor", "mismatch"])
def test_required_observation(change):
    observed = result()
    adapter = dict(format="python_pairs", calls=[{"status": "constructor_returned"}])
    if change == "optimizer":
        observed["optimizer_called"] = True
    elif change == "fingerprint":
        observed["fingerprint"] = {}
    elif change == "constructor":
        adapter["calls"][0]["status"] = "before_constructor"
    elif change == "mismatch":
        observed["differences"]["native_vs_saved"]["different_edges"] = 1
    if change:
        with pytest.raises(ValueError):
            module.require_stop(observed, adapter, {"fixture": True})
    else:
        module.require_stop(observed, adapter, {"fixture": True})


def test_existing_output_preserved(tmp_path):
    with pytest.raises(FileExistsError):
        module.run(tmp_path, tmp_path)


def test_small_frozen_worker_stops_without_clusters(tmp_path):
    root = Path(module.__file__).resolve().parent.parent
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    if not (launcher / "orthohmm/leiden_worker.py").is_file():
        pytest.skip("Requires retained frozen native worker")
    payload = tmp_path / "payload"
    payload.mkdir()
    (payload / "gene_names.txt").write_text("a\nb\nc\n")
    for name, data in (("sources", np.array([0, 1], dtype=np.int32)),
                       ("targets", np.array([1, 2], dtype=np.int32)),
                       ("weights", np.array([.5, .8], dtype=np.float64))):
        np.save(payload / (name + ".npy"), data)
    (payload / "metadata.json").write_text(json.dumps(dict(cpm_resolution=.12, seed=4,
        include_isolates=True, output_directory=str(tmp_path))))
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               PYTHONPATH=str(launcher), PYTHONHASHSEED="0")
    completed = subprocess.run([sys.executable, "-B", str(Path(module.__file__).resolve()),
        "--root", str(root), "--output", str(tmp_path), "--worker"], cwd=launcher, env=env,
        capture_output=True, text=True, timeout=90)
    assert completed.returncode == 0, completed.stderr
    observed = json.loads((payload / "preoptimizer_stop.json").read_text())
    adapter = json.loads((payload / "constructor_adapter.json").read_text())
    module.require_stop(observed, adapter, observed["saved"])
    assert observed["fingerprint"]["vertices"] == 3
    assert (payload / "worker_before.json").is_file()
    assert not (tmp_path / "orthohmm_working_res").exists()
