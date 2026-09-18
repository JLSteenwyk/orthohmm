import json
import os
from pathlib import Path
import subprocess

import pytest

from benchmark_tools import run_qfo_corrected_blast as runner


@pytest.fixture
def execution(tmp_path, monkeypatch):
    root = tmp_path / "output"
    work = root / "work"
    work.mkdir(parents=True)
    (work / "all.fa").write_text(">gene\nACDE\n")
    (work / "all.gg").write_text("species: gene\n")
    plan = {"output_root": str(root), "search_commands": {"formatdb": ["formatdb"], "blast": ["blast"]}}
    plan_path = tmp_path / "plan.json"
    plan_path.write_text(json.dumps(plan))
    runtime = tmp_path / "runtime.json"
    runtime.write_text("{}")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "180")
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    monkeypatch.setattr(runner, "verify", lambda *args: plan)
    calls = []

    def execute(argv, **kwargs):
        stage = argv[-1]
        calls.append(stage)
        Path(argv[3]).write_text("test timing")
        if stage == "formatdb":
            for suffix in ("phr", "pin", "psq"):
                (work / ("all.fa." + suffix)).write_text("test database")
        else:
            (work / "all.blast.partial").write_text("test hit")
        return subprocess.CompletedProcess(argv, 0)

    monkeypatch.setattr(runner.subprocess, "run", execute)
    return plan, plan_path, runtime, calls, execute


def test_stage_order_and_no_implicit_admission(execution):
    plan, path, runtime, calls, _ = execution
    result = runner.run(path, runtime)
    assert calls == ["formatdb", "blast"]
    assert result["status"] == "search_exited_zero_pending_query_and_database_admission"
    assert result["search_admitted"] is False and result["accuracy_admitted"] is False
    runner.check(result["blast_output"])
    with pytest.raises(FileExistsError):
        runner.run(path, runtime)


def test_check_only_does_not_execute(execution):
    _, path, runtime, calls, _ = execution
    assert runner.run(path, runtime, True)["status"] == "preflight_passed_no_search"
    assert not calls


@pytest.mark.parametrize("stage", ["formatdb", "blast"])
def test_native_failure_preserved(execution, monkeypatch, stage):
    plan, path, runtime, calls, original = execution

    def execute(argv, **kwargs):
        result = original(argv, **kwargs)
        if argv[-1] == stage:
            result.returncode = 1
        return result

    monkeypatch.setattr(runner.subprocess, "run", execute)
    with pytest.raises(RuntimeError, match="exited 1"):
        runner.run(path, runtime)
    status = json.loads((Path(plan["output_root"]) / "search_execution/status.json").read_text())
    assert status["status"] == "failed"
    assert calls == (["formatdb"] if stage == "formatdb" else ["formatdb", "blast"])


def test_missing_database_blocks_search(execution, monkeypatch):
    plan, path, runtime, calls, original = execution

    def execute(argv, **kwargs):
        result = original(argv, **kwargs)
        (Path(plan["output_root"]) / "work/all.fa.pin").unlink()
        return result

    monkeypatch.setattr(runner.subprocess, "run", execute)
    with pytest.raises(ValueError, match="Incomplete formatted"):
        runner.run(path, runtime)
    assert calls == ["formatdb"]


def test_wrong_allocation(execution, monkeypatch):
    _, path, runtime, _, _ = execution
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "32")
    with pytest.raises(ValueError, match="allocation"):
        runner.run(path, runtime)


def test_environment_excludes_search_and_loader_overrides(monkeypatch):
    monkeypatch.setenv("BLASTMAT", "/foreign")
    monkeypatch.setenv("NCBI", "/foreign")
    monkeypatch.setenv("LD_PRELOAD", "/foreign")
    assert not {"BLASTMAT", "NCBI", "LD_PRELOAD"} & runner.environment().keys()


@pytest.mark.skipif(os.environ.get("ORTHOHMM_LEGACY_BLAST_SMOKE") != "1", reason="Opt-in installed legacy engine smoke")
def test_installed_legacy_engine_smoke(tmp_path):
    from benchmark_tools.prepare_qfo_corrected_orthomcl import SOFTWARE, commands
    sequence = "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEQHFKGLVLIAFSQYLQQCPFDEHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDK"
    (tmp_path / "all.fa").write_text(">smoke_a\n" + sequence + "\n>smoke_b\n" + sequence + "\n")
    argv = commands(tmp_path, SOFTWARE / "blast-2.2.13/bin/blastall", SOFTWARE / "blast-2.2.13/bin/formatdb", 1)
    for stage in ("formatdb", "blast"):
        subprocess.run(argv[stage], cwd=tmp_path, env=runner.environment(), check=True,
                       capture_output=True, text=True, timeout=30)
    lines = (tmp_path / "all.blast.partial").read_text().splitlines()
    pairs = {tuple(line.split("\t")[:2]) for line in lines}
    assert pairs == {(a, b) for a in ("smoke_a", "smoke_b") for b in ("smoke_a", "smoke_b")}
    assert all(len(line.split("\t")) == 12 for line in lines)
