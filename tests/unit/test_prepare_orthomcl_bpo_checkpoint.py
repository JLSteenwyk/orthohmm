import json
import os
from pathlib import Path
import subprocess

import pytest

from benchmark_tools import prepare_orthomcl_bpo_checkpoint as module
from benchmark_tools.probe_orthomcl_bpo_parity import fixture


def setup(tmp_path, monkeypatch, problem=None):
    results = tmp_path / "benchmark_tools/results"
    results.mkdir(parents=True)
    for name in ("qfo_corrected_orthomcl_perl_runtime_20260918.json",
                 "qfo_corrected_orthomcl_system_helpers_20260918.json"):
        (results / name).write_text("{}")
    monkeypatch.setattr(module, "read_frozen", lambda *args: {})
    monkeypatch.setattr(module, "verify", lambda manifest: None)
    fasta, blast = fixture(tmp_path)

    def native(argv, directory, env, name):
        if problem == name:
            raise ValueError("injected native failure")
        if name == "build":
            (directory / "indexes").mkdir()
            for index in ("all_bpo.idx", "all_bpo.se"):
                (directory / "indexes" / index).write_bytes(b"mock index")
        else:
            result = {"status": "native_bpo_indexes_verified", "records": 6,
                      "queries": 3, "offset_entries_including_eof": 7,
                      "bpo_bytes": (directory / "all.bpo").stat().st_size}
            if problem in result:
                result[problem] = "bad"
            (directory / "validate.stdout").write_text(json.dumps(result))
            if problem == "mutation":
                blast.write_text(blast.read_text() + "changed\n")
    monkeypatch.setattr(module, "native_step", native)
    return fasta, blast


def test_checkpoint_content_and_index_workflow(tmp_path, monkeypatch):
    fasta, blast = setup(tmp_path, monkeypatch)
    output = tmp_path / "checkpoint"
    report = module.prepare(tmp_path, blast, fasta, output)
    assert report["status"] == "bpo_checkpoint_content_and_indexes_verified"
    assert report["content"]["bpo_pair_records"] == 6
    assert report["content"]["cutoff_excluded_pair_blocks"] == 1
    assert report["search_admitted"] is report["accuracy_admitted"] is False
    assert len(report["commands"]) == 2
    assert "run_orthomcl_perl_script.pl" in report["commands"][0]["argv"][2]
    for item in report["checked_records"] + report["outputs"]:
        module.check(item)
    with pytest.raises(FileExistsError):
        module.prepare(tmp_path, blast, fasta, output)


@pytest.mark.parametrize("problem", ["build", "validate", "status", "records", "bpo_bytes",
                                      "offset_entries_including_eof", "mutation"])
def test_failures_are_retained_without_success(tmp_path, monkeypatch, problem):
    fasta, blast = setup(tmp_path, monkeypatch, problem)
    output = tmp_path / "checkpoint"
    with pytest.raises((ValueError, RuntimeError)):
        module.prepare(tmp_path, blast, fasta, output)
    report = json.loads((output / "report.json").read_text())
    assert report["status"] == "failed"
    assert report["search_admitted"] is report["accuracy_admitted"] is False
    assert (output / "all.bpo").exists()


def test_content_failure_prevents_native_execution(tmp_path, monkeypatch):
    fasta, blast = setup(tmp_path, monkeypatch)
    blast.write_text(blast.read_text().replace("80.00", "79.00"))
    output = tmp_path / "checkpoint"
    with pytest.raises(ValueError, match="accounting"):
        module.prepare(tmp_path, blast, fasta, output)
    report = json.loads((output / "report.json").read_text())
    assert report["commands"] == [] and report["status"] == "failed"
    assert not (output / "indexes").exists()


def test_broken_symlink_is_not_overwritten(tmp_path):
    output = tmp_path / "checkpoint"
    output.symlink_to(tmp_path / "missing")
    with pytest.raises(FileExistsError):
        module.prepare(tmp_path, tmp_path / "blast", tmp_path / "fasta", output)


@pytest.mark.parametrize("code,diagnostics", [(1, b""), (0, b"warning\n"), (0, b"")])
def test_native_return_code_and_diagnostics(tmp_path, monkeypatch, code, diagnostics):
    def run(argv, **kwargs):
        kwargs["stderr"].write(diagnostics)
        return subprocess.CompletedProcess(argv, code)
    monkeypatch.setattr(module.subprocess, "run", run)
    if code or diagnostics:
        with pytest.raises(ValueError):
            module.native_step(["test"], tmp_path, {}, "build")
    else:
        module.native_step(["test"], tmp_path, {}, "build")


@pytest.mark.skipif(os.environ.get("ORTHOHMM_LEGACY_BLAST_SMOKE") != "1",
                    reason="Opt-in installed frozen OrthoMCL runtime")
def test_installed_checkpoint(tmp_path):
    fasta, blast = fixture(tmp_path)
    root = Path(__file__).resolve().parents[2]
    report = module.prepare(root, blast, fasta, tmp_path / "checkpoint")
    assert report["status"] == "bpo_checkpoint_content_and_indexes_verified"
    assert report["index_validation"]["records"] == 6
    assert report["index_validation"]["queries"] == 3
