import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import pytest

from benchmark_tools import probe_qfo_construction as probe


def test_conversion_preserves_values_and_original():
    edges = np.array([[2, 1], [0, 0]], dtype=np.int32)
    assert probe.converted_edges(edges, "original_int32") is edges
    converted = probe.converted_edges(edges, "explicit_int64")
    assert converted.dtype == np.int64 and converted.flags.c_contiguous
    assert np.array_equal(converted, edges)
    converted[0, 0] = 7
    assert edges[0, 0] == 2
    with pytest.raises(ValueError):
        probe.converted_edges(converted, "original_int32")
    with pytest.raises(ValueError):
        probe.converted_edges(edges, "unknown")


@pytest.mark.parametrize("mode", ["original_int32", "explicit_int64"])
def test_fresh_worker_constructs_without_optimization(tmp_path, mode):
    launcher = tmp_path / "benchmarks/work/publication_qfo_replay_native_v1"
    package = launcher / "orthohmm"
    package.mkdir(parents=True)
    source = Path(probe.__file__).resolve().parent.parent / "orthohmm"
    for name in ("__init__.py", "externals.py", "helpers.py", "files.py", "leiden_worker.py"):
        shutil.copyfile(source / name, package / name)
    directory = tmp_path / "output/worker"
    payload = directory / "payload"
    payload.mkdir(parents=True)
    (directory / "orthohmm_working_res").mkdir()
    (payload / "gene_names.txt").write_text("a\nb\nc\n")
    np.save(payload / "sources.npy", np.array([1, 0], dtype=np.int32))
    np.save(payload / "targets.npy", np.array([0, 0], dtype=np.int32))
    np.save(payload / "weights.npy", np.array([1., 2.]))
    (payload / "metadata.json").write_text(json.dumps({"cpm_resolution": .1, "seed": 4,
        "include_isolates": True, "output_directory": str(directory)}))
    command = [sys.executable, probe.__file__, "--root", str(tmp_path), "--output", str(tmp_path / "output"),
               "--worker-payload", str(payload), "--mode", mode]
    run = subprocess.run(command, cwd=launcher, env={**os.environ, "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}, capture_output=True, text=True, timeout=30)
    assert run.returncode == 0, run.stderr
    result = json.loads((payload / "construction.json").read_text())
    assert result["native"] == result["saved"]
    assert result["native"]["vertices"] == 3
    assert result["optimizer_called"] is False
    assert result["differences"]["native_vs_constructor"]["different_edges"] == 0
    assert result["converted_differences"]["constructor_vs_saved"]["different_edges"] == 0
    assert not (directory / "orthohmm_working_res/orthohmm_edges_clustered.txt").exists()
