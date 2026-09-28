import importlib
import json
from pathlib import Path

import pytest


@pytest.fixture
def module(monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[2] / "benchmark_tools"))
    return importlib.import_module("run_native_fas_path_panel")


def test_one_file_preserves_counts_and_fingerprint(module, monkeypatch, tmp_path):
    annotation = tmp_path / "annotations.json"
    annotation.write_text("{}")
    source = module.fingerprint(annotation)
    result = dict(annotation=source, summary=dict(proteins=3, rejected=1, maximum_paths=10**16),
                  options=dict(paths_limit=10**15), greedyfas_version="fixture")
    monkeypatch.setattr(module, "run", lambda _: result)
    output = tmp_path / "output"
    output.mkdir()
    row = module.one_file((source, output))
    assert json.loads(Path(row["result"]["path"]).read_text()) == result
    assert row["result"] == module.fingerprint(Path(row["result"]["path"]))
    assert row["summary"]["rejected"] == 1
    with pytest.raises(FileExistsError):
        module.one_file((source, output))


def test_one_file_rejects_changed_input(module, monkeypatch, tmp_path):
    path = tmp_path / "input.json"
    path.write_text("old")
    source = module.fingerprint(path)
    path.write_text("new")
    monkeypatch.setattr(module, "run", lambda _: dict(annotation=module.fingerprint(path)))
    with pytest.raises(ValueError, match="frozen environment"):
        module.one_file((source, tmp_path))


@pytest.mark.parametrize("duplicate", [False, True])
def test_panel_requires_nonempty_unique_inventory(module, tmp_path, duplicate):
    env = tmp_path / "environment.json"
    records = [dict(path="/fas_annotations/a.json")] * 2 if duplicate else []
    env.write_text(json.dumps(dict(reference_files=records)))
    out = tmp_path / "output"
    with pytest.raises(ValueError, match="inventory"):
        module.panel(env, out)
    assert not out.exists()


def test_panel_checks_frozen_bytes_before_creating_output(module, tmp_path):
    directory = tmp_path / "fas_annotations"
    directory.mkdir()
    annotation = directory / "a.json"
    annotation.write_text("{}")
    ref = module.fingerprint(annotation)
    annotation.write_text("changed")
    env = tmp_path / "environment.json"
    env.write_text(json.dumps(dict(reference_files=[ref])))
    out = tmp_path / "output"
    with pytest.raises(ValueError, match="changed"):
        module.panel(env, out)
    assert not out.exists()
