import json
import os
from pathlib import Path
import subprocess

import pytest

from benchmark_tools import stage_orthomcl_native_inputs as module


def fixture(tmp_path, monkeypatch):
    source = tmp_path / "source"
    source.mkdir()
    for key, name in module.NAMES.items():
        (source / name).write_text("one: A\ntwo: B\n" if key == "species" else "sample\n")
    records = {key: module.record(source / name) for key, name in module.NAMES.items()}
    def native(argv, directory, env, name):
        (directory / "indexes.stdout").write_text(json.dumps({"status": "test_verified"}))
    monkeypatch.setattr(module, "native_step", native)
    return records, {"status": "test_verified"}


def test_native_layout_uses_independent_copies(tmp_path, monkeypatch):
    records, indexes = fixture(tmp_path, monkeypatch)
    directory = tmp_path / "inputs"
    report = module.stage(records, directory, 2, 2, indexes)
    assert report["status"] == "native_inputs_staged_and_indexes_verified"
    for key, item in report["staged"].items():
        path = Path(item["path"])
        assert path.name == module.NAMES[key]
        assert not path.is_symlink()
        assert path.stat().st_ino != Path(records[key]["path"]).stat().st_ino
        module.check(item)
    assert not list(directory.glob("pair_*"))
    with pytest.raises(FileExistsError):
        module.stage(records, directory, 2, 2, indexes)


@pytest.mark.parametrize("problem", ["missing", "alias", "empty", "changed", "path", "count"])
def test_invalid_sources_rejected_before_staging(tmp_path, monkeypatch, problem):
    records, indexes = fixture(tmp_path, monkeypatch)
    directory = tmp_path / "inputs"
    proteins = 2
    if problem == "missing":
        records.pop("bpo")
    elif problem == "alias":
        records["offsets"] = records["bpo"]
    elif problem == "empty":
        records["bpo"]["bytes"] = 0
    elif problem == "changed":
        Path(records["bpo"]["path"]).write_text("changed")
    elif problem == "path":
        directory = tmp_path / "unsafe path"
    else:
        proteins = True
    with pytest.raises(ValueError):
        module.stage(records, directory, proteins, 2, indexes)
    assert not directory.exists()


@pytest.mark.parametrize("problem", ["proteins", "species", "indexes", "copy", "native"])
def test_partial_stage_failure_is_retained(tmp_path, monkeypatch, problem):
    records, indexes = fixture(tmp_path, monkeypatch)
    directory = tmp_path / "inputs"
    if problem == "indexes":
        indexes = {}
    elif problem == "copy":
        monkeypatch.setattr(module.shutil, "copyfile", lambda src, dest: Path(dest).write_text("corrupt"))
    elif problem == "native":
        def fail(*args):
            raise ValueError("native failure")
        monkeypatch.setattr(module, "native_step", fail)
    with pytest.raises(ValueError):
        module.stage(records, directory, 3 if problem == "proteins" else 2,
                     3 if problem == "species" else 2, indexes)
    report = json.loads((directory / "staging.json").read_text())
    assert report["status"] == "failed"
    assert report["accuracy_admitted"] is report["publication_ready"] is False
    for item in records.values():
        module.check(item)


def test_dangling_destination_symlink_rejected(tmp_path, monkeypatch):
    records, indexes = fixture(tmp_path, monkeypatch)
    directory = tmp_path / "inputs"
    directory.symlink_to(tmp_path / "missing")
    with pytest.raises(FileExistsError):
        module.stage(records, directory, 2, 2, indexes)


@pytest.mark.skipif(os.environ.get("ORTHOHMM_LEGACY_BLAST_SMOKE") != "1",
                    reason="Installed dedicated Python and native fixture required")
def test_installed_staged_native_inference(tmp_path):
    root = Path(__file__).resolve().parents[2]
    output = tmp_path / "probe"
    command = [str(root / "benchmarks/work/orthomcl_python_env_20260918/bin/python"),
               "-I", "-B", "-X", "pycache_prefix=" + str(tmp_path / "absent_cache"),
               str(root / "benchmark_tools/probe_orthomcl_staged_inference.py"),
               "--root", str(root), "--output", str(output)]
    subprocess.run(command, env={**module.environment(), "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"},
                   check=True, capture_output=True, timeout=120)
    report = json.loads((output / "report.json").read_text())
    assert report["status"] == "staged_native_fixture_partition_and_input_preservation_verified"
    assert report["groups"] == 12 and report["grouped_proteins"] == 41
    assert report["staging"]["index_validation"]["records"] == 379
