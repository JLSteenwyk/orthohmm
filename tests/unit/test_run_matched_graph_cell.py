import json
from pathlib import Path
import sys

import pytest

from benchmark_tools import run_matched_graph_cell as runner


@pytest.fixture
def prepared(tmp_path, monkeypatch):
    numeric = tmp_path / "numeric.json"
    numeric.write_text("{}")
    metadata = tmp_path / "metadata.json"
    metadata.write_text("{}")
    staging = tmp_path / "runtime/staging.json"
    staging.parent.mkdir()
    staging.write_text("{}")
    python = staging.parent / "venv_clean/bin/python"
    python.parent.mkdir(parents=True)
    python.symlink_to(sys.executable)
    install = tmp_path / "installation.json"
    install.write_text(json.dumps(dict(checked_records=[runner.record(staging)])))
    manifest = tmp_path / "manifest.json"
    cell = dict(numeric=runner.record(numeric), sources=[], arm="hmm", seed=20261106, condition="fixture")
    manifest.write_text(json.dumps(dict(cells=[cell] * 70, protocol=runner.record(metadata),
                                       search_result=runner.record(metadata), graph_settings={"fixture": True})))
    monkeypatch.setattr(runner, "MANIFEST_SHA", runner.record(manifest)["sha256"])
    monkeypatch.setattr(runner, "INSTALL_SHA", runner.record(install)["sha256"])
    return manifest, install


def test_failure_retained_no_retry(prepared, tmp_path, monkeypatch):
    calls = []
    def fail(*args):
        calls.append(args)
        raise RuntimeError("native failed")
    monkeypatch.setattr(runner, "execute", fail)
    output = tmp_path / "output"
    with pytest.raises(RuntimeError, match="native failed"):
        runner.run(*prepared, 0, output)
    receipt = json.loads((output / "execution.json").read_text())
    assert receipt["status"] == "failed" and receipt["attempt"] == 1
    assert "native failed" in receipt["error"]
    with pytest.raises(FileExistsError):
        runner.run(*prepared, 0, output)
    assert len(calls) == 1


def test_success_and_private_copy(prepared, tmp_path, monkeypatch):
    def fake(command, output, name, env, deadline, stages):
        assert "-I" in command and "PYTHONPATH" not in env
        private = Path(command[command.index("--numeric") + 1])
        assert private.parent == output
        graph = Path(command[command.index("--output") + 1])
        graph.mkdir()
        checkpoint = graph / "checkpoint.json"
        checkpoint.write_text("{}")
        (graph / "receipt.json").write_text(json.dumps(dict(settings={"fixture": True},
            numeric=runner.record(private), outputs=[], checkpoint_manifest=runner.record(checkpoint))))
    monkeypatch.setattr(runner, "execute", fake)
    output = tmp_path / "output"
    runner.run(*prepared, 0, output)
    receipt = json.loads((output / "execution.json").read_text())
    assert receipt["status"] == "native_completed_pending_independent_readback"
    assert receipt["native_receipt"]["path"].endswith("graph/receipt.json")


def test_changed_input_and_bad_index(prepared, tmp_path):
    with pytest.raises(ValueError, match="index"):
        runner.run(*prepared, -1, tmp_path / "bad_index")
    assert not (tmp_path / "bad_index").exists()
    (tmp_path / "numeric.json").write_text("changed")
    with pytest.raises(ValueError):
        runner.run(*prepared, 0, tmp_path / "changed")
    assert json.loads((tmp_path / "changed/execution.json").read_text())["status"] == "failed"


def test_unknown_manifest_rejected(prepared, tmp_path):
    prepared[0].write_text(prepared[0].read_text() + " ")
    with pytest.raises(ValueError, match="Unrecognized"):
        runner.run(*prepared, 0, tmp_path / "output")
    assert not (tmp_path / "output").exists()
