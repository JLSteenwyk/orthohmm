import json
from pathlib import Path
import subprocess

import pytest

from benchmark_tools import run_qfo_sequence_search_control as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_prepare_qfo_sequence_search_control import fixture


@pytest.mark.parametrize("problem", [None, "status", "genes", "version", "command", "path", "order"])
def test_verify_frozen_plan(tmp_path, monkeypatch, problem):
    stage, inventory = fixture()
    path = tmp_path / "manifest.json"
    report = {"status": "corrected_qfo_search_control_prepared_unrun", "accuracy_evaluated": False,
        "execution_authorized": False, "genes": 984137, "proteomes": 78, "source": {}, "helpers": [],
        "queries": {"path": str(tmp_path / "queries.fasta")}, "gene_metadata": {"path": str(tmp_path / "gene_metadata.json")},
        "diamond": {"path": "/diamond"}, "diamond_version": "diamond version 2.1.11",
        "checked_records": [{"path": "/stage"}, {"path": "/inventory"}], "inputs": stage["input_fastas"]}
    report["searches"] = module.search_plan(report["inputs"], Path("/diamond"), tmp_path / "queries.fasta", tmp_path)
    monkeypatch.setattr(module, "read_frozen", lambda p, sha: {str(path): report, "/stage": stage, "/inventory": inventory}[str(p)])
    monkeypatch.setattr(module, "check", lambda item: None)
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: "diamond version 2.1.11\n")
    if problem == "status":
        report["execution_authorized"] = True
    elif problem == "genes":
        report["genes"] = 976504
    elif problem == "version":
        report["diamond_version"] = "diamond version 2.0.13"
    elif problem == "command":
        report["searches"][0]["search"].append("--fast")
    elif problem == "path":
        report["queries"]["path"] = "/elsewhere"
    elif problem == "order":
        report["inputs"] = list(reversed(report["inputs"]))
    if problem:
        with pytest.raises(ValueError):
            module.verify_plan(path)
    else:
        assert len(module.verify_plan(path)["searches"]) == 78


@pytest.mark.parametrize("failure", [None, "makedb", "search", "postflight"])
def test_run_retains_outcomes_and_rejects_reuse(tmp_path, monkeypatch, failure):
    manifest = tmp_path / "manifest.json"
    manifest.write_text("{}")
    searches = []
    for i in range(2):
        directory = tmp_path / f"target_{i:02d}"
        directory.mkdir()
        searches.append({"index": i, "output": str(directory / "hits.tsv"),
                         "makedb": ["fixture", "makedb"], "search": ["fixture", "search"]})
    calls = []
    def verify(path):
        calls.append(path)
        if failure == "postflight" and len(calls) > 1:
            raise ValueError("Changed evidence")
        return {"searches": searches}
    monkeypatch.setattr(module, "verify_plan", verify)
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "32")
    monkeypatch.setenv("SLURM_JOB_NODELIST", "bizon")
    def phase(argv, directory, name, env):
        assert env["OMP_NUM_THREADS"] == "1"
        for filename in (name + ".log", name + ".time.log", "hits.tsv" if name == "search" else "target.dmnd"):
            (directory / filename).write_text("fixture\n")
        return {"argv": argv, "exit_code": int(failure == name), "wall_s": .1,
                "log": record(directory / (name + ".log")), "gnu_time": record(directory / (name + ".time.log"))}
    monkeypatch.setattr(module, "run_phase", phase)
    if failure:
        with pytest.raises((RuntimeError, ValueError)):
            module.run(manifest)
    else:
        module.run(manifest)
    result = json.loads((tmp_path / "execution.json").read_text())
    assert result["status"] == ("failed" if failure else "complete_pending_numeric_validation")
    assert result["numeric_validated"] is False and result["accuracy_evaluated"] is False
    assert len(result["targets"]) == (1 if failure in ("makedb", "search") else 2)
    monkeypatch.setattr(module, "verify_plan", lambda path: {"searches": searches})
    with pytest.raises(FileExistsError):
        module.run(manifest)


def test_preflight_and_missing_allocation(tmp_path, monkeypatch):
    path = tmp_path / "manifest.json"
    monkeypatch.setattr(module, "verify_plan", lambda p: {"searches": []})
    assert module.run(path, check_only=True)["status"] == "corrected_search_preflight_verified"
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    with pytest.raises(ValueError):
        module.run(path)
    assert not (tmp_path / "execution.json").exists()


def test_batch_resources():
    path = Path(module.__file__).parent / "results/qfo_sequence_search_batch_20260918.sh"
    subprocess.run(["bash", "-n", str(path)], check=True)
    for setting in ("--nodelist=bizon", "--cpus-per-task=32", "--mem=192G", "--time=7-00:00:00", "--no-requeue"):
        assert "#SBATCH " + setting in path.read_text()
