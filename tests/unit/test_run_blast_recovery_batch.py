import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import run_blast_recovery_batch as module


def test_only_query_and_output_change():
    original = ["blastall", "-p", "blastp", "-i", "all.fa", "-d", "full.db", "-a", "180", "-o", "old", "-e", "1e-5"]
    changed = module.command_for(original, "subset.fa", "new")
    assert [(i, a, b) for i, (a, b) in enumerate(zip(original, changed)) if a != b] == [(4, "all.fa", "subset.fa"), (10, "old", "new")]


@pytest.mark.parametrize("command", [["blastall"], ["-i", "a", "-o", "b", "-d", "c", "-a", "1"],
                                      ["-i", "a", "-i", "b", "-o", "c", "-d", "d", "-a", "180"]])
def test_malformed_command_rejected(command):
    with pytest.raises(ValueError):
        module.command_for(command, "q", "out")


def fixture(tmp_path, monkeypatch, returncode=0, missing=False):
    directory = tmp_path / "batch"
    preflight = dict(directory=str(directory), index=0, command=["blastall"], cwd=str(tmp_path),
                     environment={}, checked_records=[], plan_path="plan", runtime_path="runtime")
    def native(command, **kwargs):
        (directory / "blast.time.txt").write_text("timing\n")
        kwargs["stdout"].write(b"diagnostic\n")
        if not missing:
            (directory / "hits.blast.partial").write_bytes(b"native bytes\n")
        return SimpleNamespace(returncode=returncode)
    monkeypatch.setattr(module.subprocess, "run", native)
    monkeypatch.setattr(module, "verify", lambda *args: {})
    return preflight


def test_success_durable_status_does_not_admit_search(tmp_path, monkeypatch):
    preflight = fixture(tmp_path, monkeypatch)
    result = module.execute(preflight, "job", "array")
    directory = Path(preflight["directory"])
    assert result == json.loads((directory / "status.json").read_text())
    assert result["status"] == "native_batch_completed_pending_admission"
    assert result["reuse_authorized"] is result["search_admitted"] is False
    assert (directory / "hits.blast").read_bytes() == b"native bytes\n"
    assert not (directory / "hits.blast.partial").exists()
    with pytest.raises(FileExistsError):
        module.execute(preflight, "job", "array")


@pytest.mark.parametrize("failure", ["exit", "missing", "verification", "interrupt"])
def test_failures_preserved_without_success_marker(tmp_path, monkeypatch, failure):
    preflight = fixture(tmp_path, monkeypatch, returncode=1 if failure == "exit" else 0, missing=failure == "missing")
    if failure == "verification":
        def changed(*args):
            raise ValueError("Changed runtime")
        monkeypatch.setattr(module, "verify", changed)
    elif failure == "interrupt":
        def interrupted(*args, **kwargs):
            raise KeyboardInterrupt()
        monkeypatch.setattr(module.subprocess, "run", interrupted)
    with pytest.raises((RuntimeError, ValueError, KeyboardInterrupt)):
        module.execute(preflight, "job", "array")
    directory = Path(preflight["directory"])
    assert json.loads((directory / "status.json").read_text())["status"] == "failed_or_interrupted"
    assert not (directory / "hits.blast").exists()


def test_atomic_status_syncs_file_then_directory(tmp_path, monkeypatch):
    seen = []
    fsync = module.os.fsync
    def sync(fd):
        seen.append(fd)
        return fsync(fd)
    monkeypatch.setattr(module.os, "fsync", sync)
    path = tmp_path / "status.json"
    module.save_status(path, {"status": "test"})
    assert len(seen) == 2
    assert json.loads(path.read_text()) == {"status": "test"}
    assert not path.with_suffix(".tmp").exists()


def test_unscheduled_run_rejected_before_preparation(tmp_path, monkeypatch):
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
    monkeypatch.setattr(module, "prepare", lambda *args: pytest.fail("Prepared without allocation"))
    with pytest.raises(ValueError, match="scheduled"):
        module.run(tmp_path, 0)
