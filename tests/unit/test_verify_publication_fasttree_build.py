import json
from types import SimpleNamespace

import pytest

from benchmark_tools import verify_publication_fasttree_build as module
from benchmark_tools.acquire_publication_fasttree import identity


def test_refuse_existing_output(tmp_path):
    with pytest.raises(FileExistsError):
        module.run(tmp_path, tmp_path / "missing", tmp_path / "mafft", tmp_path)


@pytest.mark.parametrize("data", [dict(status="complete"), dict(status="failed", error="native failure"),
    dict(status="failed", error="Unexpected built FastTree identity/precision", builds=[]),
    dict(status="failed", error="Unexpected built FastTree identity/precision", builds=[{}, {}], repeat_build_byte_equal=False)])
def test_wrong_prior_state(tmp_path, data):
    report = tmp_path / "prior.json"
    report.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        module.run(tmp_path, report, tmp_path / "mafft", tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_fixture_failure_retains_report_without_rebuild(tmp_path, monkeypatch):
    binary = tmp_path / "binary"
    binary.write_text("test fixture bytes, never executed")
    item = identity(binary)
    prior = dict(status="failed", error="Unexpected built FastTree identity/precision",
        builds=[item, item], repeat_build_byte_equal=True, sources=[item], compiler=dict(binary=item),
        linker=item, compile_environment=dict(PATH="/usr/bin:/bin"))
    path = tmp_path / "prior.json"
    path.write_text(json.dumps(prior))
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k: SimpleNamespace(
        returncode=0, stdout="", stderr="FastTree 2.2.0 Double precision:\n"))
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: "test metadata")
    calls = []

    def fail(command, *args, **kwargs):
        calls.append(command)
        raise RuntimeError("injected fixture failure")

    monkeypatch.setattr(module, "execute", fail)
    output = tmp_path / "out"
    with pytest.raises(RuntimeError, match="injected"):
        module.run(tmp_path, path, binary, output)
    assert len(calls) == 1
    assert "benchmark_tools.verify_frozen_phylogeny_install" in calls[0]
    result = json.loads((output / "report.json").read_text())
    assert result["status"] == "failed"
    assert result["builds_repeated"] is False
    assert json.loads(path.read_text()) == prior
