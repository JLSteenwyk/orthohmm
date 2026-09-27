import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from benchmark_tools.run_integrated_publication_workflow import install_commands, record, relocate_assets, stage, validate_assets, validate_data


def fixture_manifest(tmp_path):
    def item(name):
        path = tmp_path / name
        path.write_text(name)
        return record(path)
    return dict(dataset="installation_fixture", genes=16,
                fasta=[item(f"S{i}.fa") for i in range(4)],
                references=[item(f"R{i}.txt") for i in range(3)], uncertain=[])


def test_fixture_scope(tmp_path):
    validate_data(fixture_manifest(tmp_path))


def test_fixture_cannot_be_claimed_as_full_benchmark(tmp_path):
    data = fixture_manifest(tmp_path)
    data["dataset"] = "orthobench"
    with pytest.raises(ValueError, match="scope or dimensions"):
        validate_data(data)


def test_changed_input_rejected(tmp_path):
    data = fixture_manifest(tmp_path)
    Path(data["fasta"][0]["path"]).write_text("changed")
    with pytest.raises(ValueError, match="Changed pinned"):
        validate_data(data)


def test_duplicate_input_rejected(tmp_path):
    data = fixture_manifest(tmp_path)
    data["fasta"][1] = data["fasta"][0]
    with pytest.raises(ValueError, match="Duplicate input"):
        validate_data(data)


def test_installs_are_offline_and_hash_required(tmp_path):
    commands = install_commands(Path("/base"), Path("/installer"), Path("/wheels"), Path("/lock"), tmp_path / "env")
    assert "--without-pip" in commands[0]
    assert all(x in commands[1] for x in ("--no-index", "--require-hashes", "--only-binary=:all:"))
    assert commands[2][-2:] == ["pip", "check"]


def test_completed_stage_records_success(tmp_path):
    result = stage(tmp_path, "success", [sys.executable, "-I", "-c", "print('complete')"], dict(os.environ), 10)
    assert result["returncode"] == 0
    assert (tmp_path / "success.log").read_text().strip() == "complete"


def test_failed_stage_does_not_retry(tmp_path):
    with pytest.raises(RuntimeError, match="without retry"):
        stage(tmp_path, "failure", [sys.executable, "-I", "-c", "raise SystemExit(9)"], dict(os.environ), 10)
    assert json.loads((tmp_path / "failure_finished.json").read_text())["returncode"] == 9
    with pytest.raises(FileExistsError):
        stage(tmp_path, "failure", [sys.executable, "-I", "-c", "pass"], dict(os.environ), 10)


def test_timeout_is_preserved(tmp_path):
    with pytest.raises(subprocess.TimeoutExpired):
        stage(tmp_path, "timeout", [sys.executable, "-I", "-c", "import time; time.sleep(10)"], dict(os.environ), 0.05)
    assert json.loads((tmp_path / "timeout_failed.json").read_text()) == dict(timeout=True, retry=False)


def test_changed_assets_manifest_rejected(tmp_path):
    (tmp_path / "copied_assets.json").write_text("[]")
    with pytest.raises(ValueError, match="native-asset manifest"):
        validate_assets(tmp_path, tmp_path)


def test_asset_relocation_never_overwrites(tmp_path):
    with pytest.raises(FileExistsError):
        relocate_assets(tmp_path, tmp_path)


def test_absolute_mafft_links_become_relative_only_in_new_copy(tmp_path, monkeypatch):
    from benchmark_tools import run_integrated_publication_workflow as module
    source, output = tmp_path / "original", tmp_path / "relocated"
    for name in ("wheels", "benchmark_tools", "mafft/bin", "mafft/libexec/mafft"):
        (source / name).mkdir(parents=True)
    for name in ("FastTree", "copied_assets.json", "manifest.json"):
        (source / name).write_text("fixture")
    for name in ("mafft-distance", "mafft-profile"):
        target = source / "mafft/libexec/mafft" / name
        target.write_text(name)
        (source / "mafft/bin" / name).symlink_to(target)
    original_record = module.record
    def fake_record(path):
        value = original_record(path)
        if Path(path) == source / "copied_assets.json":
            value["sha256"] = "8f8be7f1609d549da79a4ec4231e937a08833a7b8acc415d03e66b16274d5067"
        return value
    monkeypatch.setattr(module, "record", fake_record)
    changes = relocate_assets(source, output)
    assert len(changes) == 2
    for name in ("mafft-distance", "mafft-profile"):
        assert os.path.isabs(os.readlink(source / "mafft/bin" / name))
        assert os.readlink(output / "mafft/bin" / name) == "../libexec/mafft/" + name
        assert (output / "mafft/bin" / name).read_text() == name
