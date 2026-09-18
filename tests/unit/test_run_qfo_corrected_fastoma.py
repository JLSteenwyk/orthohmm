from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import run_qfo_corrected_fastoma as runner
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def staged(tmp_path):
    directory = tmp_path / "benchmarks/work/qfo_corrected_fastoma_inputs_20260918"
    (directory / "proteome").mkdir(parents=True)
    source = tmp_path / "source"
    source.mkdir()
    copies = []
    for i in range(79):
        original = source / (f"sp{i}.fasta" if i < 78 else "SpeciesTree_rooted_node_labels.txt")
        target = directory / "proteome" / f"sp{i}.fa" if i < 78 else directory / "species_tree.nwk"
        original.write_text(f"input{i}\n")
        target.write_bytes(original.read_bytes())
        copies.append({"source": record(original), "staged": record(target)})
    return {"status": "corrected_fastoma_inputs_staged_pending_launch_freeze",
            "execution_authorized": False, "accuracy_evaluated": False, "publication_ready": False,
            "input_proteins": 984137, "input_species": 78, "input_directory": str(directory),
            "copies": copies, "checked_records": [r["source"] for r in copies]}


def test_stage_exact_membership(tmp_path):
    stage = staged(tmp_path)
    assert runner.validate_stage(stage, tmp_path) == Path(stage["input_directory"])


@pytest.mark.parametrize("change", ["species", "proteins", "status", "authorized", "missing",
                                   "source", "checksum", "duplicate", "extra", "symlink", "directory"])
def test_invalid_staging_rejected(tmp_path, change):
    stage = staged(tmp_path)
    if change in ("species", "proteins"):
        stage["input_" + change] -= 1
    elif change == "status":
        stage["status"] = "pending"
    elif change == "authorized":
        stage["execution_authorized"] = True
    elif change == "missing":
        stage["copies"].pop()
    elif change == "source":
        stage["checked_records"].pop()
    elif change == "checksum":
        stage["copies"][0]["staged"]["sha256"] = "changed"
    elif change == "duplicate":
        stage["copies"][0] = deepcopy(stage["copies"][1])
    elif change == "directory":
        stage["input_directory"] = str(tmp_path)
    elif change == "extra":
        (Path(stage["input_directory"]) / "extra").touch()
    else:
        target = Path(stage["copies"][0]["staged"]["path"])
        target.unlink()
        target.symlink_to(stage["copies"][0]["source"]["path"])
    with pytest.raises(ValueError):
        runner.validate_stage(stage, tmp_path)


def test_command_is_fresh_and_fixed_config():
    argv = runner.command(Path("/root"), Path("/stage"), Path("/out"), Path("/workflow"))
    assert argv[:5] == ["/home/bizon/bin/nextflow", "-C", "/root/benchmark_tools/fastoma_corrected_execution.config",
                        "run", "/workflow/FastOMA.nf"]
    assert "-resume" not in argv and "-profile" not in argv
    assert argv[argv.index("--max_cpus") + 1] == "180"
    assert argv[argv.index("--max_memory") + 1] == "700.GB"
    assert argv[argv.index("--input_folder") + 1] == "/stage"
    assert argv[-2:] == ["-work-dir", "/out/work"]


def test_environment_excludes_startup_and_container_overrides(monkeypatch):
    for key in ("BASH_ENV", "JAVA_TOOL_OPTIONS", "LD_PRELOAD", "DOCKER_HOST", "NXF_OPTS"):
        monkeypatch.setenv(key, "untrusted")
    env = runner.environment()
    assert not {"BASH_ENV", "JAVA_TOOL_OPTIONS", "LD_PRELOAD", "DOCKER_HOST"}.intersection(env)
    assert env["NXF_OPTS"] == "-Xmx1g" and env["NXF_OFFLINE"] == "true"


def setup_run(monkeypatch, tmp_path, exit_code=0, drift=False):
    output = tmp_path / "output"
    base = {"cwd": str(output / "run"), "output_root": str(output), "native_argv": ["nextflow", "run"],
            "environment": {"LC_ALL": "C"}}
    calls = []
    def preflight(*args):
        calls.append(True)
        report = deepcopy(base)
        if drift and len(calls) > 1:
            report["native_argv"] = ["changed"]
        return report, {}
    def native(argv, **kwargs):
        assert argv[:3] == ["/usr/bin/time", "-v", "-o"]
        (output / "time.txt").write_text("timing\n")
        (output / "run/trace.txt").write_text("trace\n")
        (output / "run/.nextflow.log").write_text("log\n")
        (output / "output").mkdir()
        (output / "output/orthologs.tsv.gz").write_bytes(b"mock, not real biological output")
        return SimpleNamespace(returncode=exit_code)
    monkeypatch.setattr(runner, "preflight", preflight)
    monkeypatch.setattr(runner.subprocess, "run", native)
    monkeypatch.setattr(runner.os, "uname", lambda: SimpleNamespace(nodename="bizon"))
    for key, value in {"SLURM_CPUS_PER_TASK": "180", "SLURM_MEM_PER_NODE": "737280", "SLURM_JOB_ID": "123"}.items():
        monkeypatch.setenv(key, value)
    for key in ("DOCKER_HOST", "DOCKER_CONTEXT"):
        monkeypatch.delenv(key, raising=False)
    return output, calls


def test_success_remains_unadmitted(tmp_path, monkeypatch):
    output, calls = setup_run(monkeypatch, tmp_path)
    result = runner.run(tmp_path, tmp_path / "stage.json", "sha", 1)
    assert len(calls) == 2
    assert result["status"] == "process_succeeded_pending_native_admission"
    assert result["native_outputs_validated"] is False and result["accuracy_admitted"] is False
    assert json.loads((output / "execution.json").read_text()) == result
    with pytest.raises(FileExistsError):
        runner.run(tmp_path, tmp_path / "stage.json", "sha", 1)


@pytest.mark.parametrize("exit_code,drift", [(7, False), (0, True)])
def test_failed_execution_is_retained(tmp_path, monkeypatch, exit_code, drift):
    output, _ = setup_run(monkeypatch, tmp_path, exit_code, drift)
    with pytest.raises((ValueError, RuntimeError)):
        runner.run(tmp_path, tmp_path / "stage.json", "sha", 1)
    report = json.loads((output / "execution.json").read_text())
    assert report["status"] == "failed" and report["exit_code"] == exit_code
    assert report["accuracy_admitted"] is False


def test_check_only_creates_no_output(tmp_path, monkeypatch):
    output, calls = setup_run(monkeypatch, tmp_path)
    monkeypatch.delenv("SLURM_JOB_ID")
    result = runner.run(tmp_path, tmp_path / "stage.json", "sha", 1, check_only=True)
    assert result["status"] == "fastoma_preflight_passed_no_inference"
    assert not output.exists() and len(calls) == 1


@pytest.mark.parametrize("key", ["SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE", "SLURM_JOB_ID"])
def test_wrong_allocation_stops_before_preflight(tmp_path, monkeypatch, key):
    output, calls = setup_run(monkeypatch, tmp_path)
    monkeypatch.delenv(key)
    with pytest.raises(ValueError, match="allocation"):
        runner.run(tmp_path, tmp_path / "stage.json", "sha", 1)
    assert calls == [] and not output.exists()
