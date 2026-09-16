from pathlib import Path
import json
import os
import shutil
import subprocess
import sys

import numpy as np
import pytest

from benchmark_tools.repeat_qfo_saved_graph import check_worker, mapped_libraries, record
from benchmark_tools import repeat_qfo_saved_graph as diagnostic


def test_mapped_library_inventory_ignores_addresses_and_duplicates():
    maps = """001-002 r--p 0000 08:01 1 /lib/libx.so.1
002-003 r-xp 0010 08:01 1 /lib/libx.so.1
003-004 r--p 0000 08:01 2 /lib/another library.so
004-005 rw-p 0000 00:00 0 [heap]
005-006 r--p 0000 08:01 3 /data/weights.npy
"""
    assert mapped_libraries(maps) == [Path("/lib/another library.so"), Path("/lib/libx.so.1")]


def test_deleted_loaded_library_rejected():
    with pytest.raises(ValueError, match="deleted"):
        mapped_libraries("001-002 r-xp 0000 08:01 1 /lib/libx.so (deleted)")


def test_record_follows_symlink_to_preserved_input(tmp_path):
    source = tmp_path / "saved.npy"
    source.write_bytes(b"example")
    link = tmp_path / "payload.npy"
    link.symlink_to(source)
    assert record(source) == record(link)


@pytest.mark.parametrize("problem", [None, "seed", "cwd", "source", "environment", "libraries", "affinity"])
def test_worker_identity_guards(problem):
    launcher = Path("/launcher")
    payload = Path("/out/repeat_0/payload")
    overrides = {"PYTHONHASHSEED": "0"}
    snapshot = {"status": "before_native_clustering", "accuracy_evaluated": False,
        "metadata": {"cpm_resolution": .1, "seed": 4, "include_isolates": True, "output_directory": "/out/repeat_0"},
        "cwd": str(launcher), "environment": dict(overrides), "cpu_affinity": [0], "native_libraries": [{}],
        "modules": {name: {"path": str(launcher / (name.replace(".", "/") + ".py"))}
                    for name in ("orthohmm.leiden_worker", "orthohmm.externals", "orthohmm.helpers")}}
    if problem == "seed":
        snapshot["metadata"]["seed"] = 0
    elif problem == "cwd":
        snapshot["cwd"] = "/other"
    elif problem == "source":
        snapshot["modules"]["orthohmm.externals"]["path"] = "/other/externals.py"
    elif problem == "environment":
        snapshot["environment"]["PYTHONHASHSEED"] = "1"
    elif problem == "libraries":
        snapshot["native_libraries"] = []
    elif problem == "affinity":
        snapshot["cpu_affinity"] = []
    if problem:
        with pytest.raises(ValueError):
            check_worker(snapshot, launcher, payload, overrides)
    else:
        check_worker(snapshot, launcher, payload, overrides)


@pytest.mark.skipif(not Path("/proc/self/maps").exists(), reason="Linux loaded-library diagnostic")
@pytest.mark.parametrize("explicit_affinity", [False, True])
def test_instrumented_real_worker_exits_and_preserves_isolate(tmp_path, explicit_affinity):
    root = tmp_path
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    package = launcher / "orthohmm"
    package.mkdir(parents=True)
    source = Path(diagnostic.__file__).resolve().parent.parent / "orthohmm"
    for name in ("__init__.py", "externals.py", "helpers.py", "files.py", "leiden_worker.py"):
        shutil.copyfile(source / name, package / name)
    output = root / "output"
    directory = output / "repeat_0"
    payload = directory / "payload"
    payload.mkdir(parents=True)
    (directory / "orthohmm_working_res").mkdir()
    (payload / "gene_names.txt").write_text("a\nb\nc\n")
    np.save(payload / "sources.npy", np.array([0], dtype=np.int32), allow_pickle=False)
    np.save(payload / "targets.npy", np.array([1], dtype=np.int32), allow_pickle=False)
    np.save(payload / "weights.npy", np.array([1.], dtype=np.float64), allow_pickle=False)
    (payload / "metadata.json").write_text(json.dumps({"cpm_resolution": .1, "seed": 4,
        "include_isolates": True, "output_directory": str(directory)}))
    overrides = {"PYTHONPATH": str(launcher), "PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1",
                 "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    inherited = sorted(os.sched_getaffinity(0))
    command = [sys.executable, diagnostic.__file__, "--root", str(root), "--output", str(output),
               "--worker-payload", str(payload)]
    if explicit_affinity:
        command += ["--cpu-affinity", str(inherited[0])]
    run = subprocess.run(command, cwd=launcher, env={**os.environ, **overrides},
                         capture_output=True, text=True, timeout=30)
    assert run.returncode == 0, run.stderr
    snapshot = json.loads((payload / "worker_before.json").read_text())
    check_worker(snapshot, launcher, payload, overrides)
    assert snapshot["inherited_cpu_affinity"] == inherited
    assert snapshot["requested_cpu_affinity"] == (inherited[:1] if explicit_affinity else None)
    assert snapshot["cpu_affinity"] == (inherited[:1] if explicit_affinity else inherited)
    assert sorted(os.sched_getaffinity(0)) == inherited
    assert any("_c_leiden" in name for name in snapshot["modules"])
    assert any("_igraph" in name for name in snapshot["modules"])
    assert snapshot["inputs"] == [record(payload / name) for name in ("gene_names.txt", "sources.npy", "targets.npy", "weights.npy")]
    assert set((directory / "orthohmm_working_res/orthohmm_edges_clustered.txt").read_text().splitlines()) == {"a b", "c"}


def test_affinity_panel_is_prespecified_alternating_subset():
    arms = diagnostic.affinity_arms(range(80, 8, -2))
    assert [name for name, _ in arms] == ["one_cpu_0", "all_cpus_0", "one_cpu_1", "all_cpus_1"]
    assert arms[0][1] == arms[2][1] == [10]
    assert arms[1][1] == arms[3][1] == list(range(10, 74, 2))


def test_affinity_panel_rejects_insufficient_distinct_cpus():
    with pytest.raises(ValueError, match="32 allocated"):
        diagnostic.affinity_arms([0] * 32)


@pytest.mark.parametrize("requested", [[], [1, 1], [9], [-1]])
def test_invalid_affinity_does_not_call_setter(monkeypatch, requested):
    monkeypatch.setattr(os, "sched_getaffinity", lambda pid: {1, 3})
    def forbidden(*args):
        pytest.fail("Invalid CPU request reached setter")
    monkeypatch.setattr(os, "sched_setaffinity", forbidden)
    with pytest.raises(ValueError, match="outside inherited"):
        diagnostic.set_worker_affinity(requested)


def test_affinity_setter_verifies_actual_result(monkeypatch):
    monkeypatch.setattr(os, "sched_getaffinity", lambda pid: {1, 3})
    monkeypatch.setattr(os, "sched_setaffinity", lambda pid, cpus: None)
    with pytest.raises(ValueError, match="Actual worker affinity"):
        diagnostic.set_worker_affinity([1])
