import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from benchmark_tools import probe_qfo_direct_graph as probe


@pytest.mark.parametrize("mode", probe.MODES)
def test_fresh_worker_before_after_weights(tmp_path, mode):
    launcher = tmp_path / "benchmarks/work/publication_qfo_replay_native_v1"
    package = launcher / "orthohmm"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text("")
    (package / "leiden_worker.py").write_text("def main(*args):\n    raise AssertionError('Worker must not execute')\n")
    payload = tmp_path / "payload"
    payload.mkdir()
    (payload / "gene_names.txt").write_text("a\nb\nc\nisolate\n")
    np.save(payload / "sources.npy", np.array([0, 1, 2], dtype=np.int32))
    np.save(payload / "targets.npy", np.array([1, 2, 2], dtype=np.int32))
    np.save(payload / "weights.npy", np.array([.1, .2, .3], dtype=np.float64))
    env = {**os.environ, "PYTHONPATH": str(launcher), "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
    run = subprocess.run([sys.executable, str(Path(probe.__file__).resolve()), "--root", str(tmp_path),
        "--output", str(tmp_path / "unused"), "--worker-payload", str(payload), "--mode", mode],
        cwd=launcher, env=env, text=True, capture_output=True)
    assert run.returncode == 0, run.stderr
    result = json.loads((payload / "result.json").read_text())
    before = result["before_weights"]
    assert before["vertices"] == 4 and before["edges"] == 3 and before["directed"] is False
    assert result["after_weights"]["fingerprint"] == result["saved"]
    for stage in (before, result["after_weights"]):
        for key in ("native_vs_saved", "constructor_vs_saved", "native_vs_constructor"):
            assert stage["differences"][key] == {"different_edges": 0, "examples": []}
    assert result["optimizer_called"] is False and result["accuracy_evaluated"] is False
    snapshot = json.loads((payload / "worker_before.json").read_text())
    assert len(snapshot["cpu_affinity"]) == 1
    assert ("orthohmm.leiden_worker" in snapshot["modules"]) == (mode == "frozen_imports")
    assert not (tmp_path / "unused").exists()
    assert not list(tmp_path.rglob("*clustered*"))
