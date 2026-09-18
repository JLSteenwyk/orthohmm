import copy
import json
import os
from pathlib import Path

import pytest

from benchmark_tools import admit_qfo_corrected_bpo as module


def preparation():
    scheduler = {"JobIDRaw": "1", "State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon",
                 "AllocCPUS": "2", "ReqMem": "64G"}
    runtime = {"status": "dedicated_bpo_python_runtime_verified", "manifest": {"runtime": 1}, "mapped_files": []}
    report = {"status": "corrected_bpo_checkpoint_prepared_pending_admission", "job_id": "1",
              "accuracy_admitted": False, "publication_ready": False, "source": {"source": 1},
              "checkpoint": {"checkpoint": 1}, "started_epoch": 1, "finished_epoch": 2,
              "runtime_before": runtime, "runtime_after": copy.deepcopy(runtime)}
    return report, scheduler


def validate(report, scheduler):
    module.validate_preparation(report, scheduler, {"source": 1}, {"checkpoint": 1}, {"runtime": 1})


def test_preparation_identity():
    validate(*preparation())


@pytest.mark.parametrize("key,value", [("status", "running"), ("job_id", "2"),
    ("accuracy_admitted", True), ("publication_ready", True), ("source", {}), ("checkpoint", {}),
    ("started_epoch", float("nan")), ("finished_epoch", 0), ("started_epoch", True)])
def test_invalid_preparation(tmp_path, key, value):
    report, scheduler = preparation()
    report[key] = value
    with pytest.raises(ValueError):
        validate(report, scheduler)


@pytest.mark.parametrize("key", ["State", "ExitCode", "NodeList", "AllocCPUS", "ReqMem"])
def test_wrong_scheduler(key):
    report, scheduler = preparation()
    scheduler[key] = "wrong"
    with pytest.raises(ValueError):
        validate(report, scheduler)


@pytest.mark.parametrize("stage", ["runtime_before", "runtime_after"])
def test_changed_runtime_binding(stage):
    report, scheduler = preparation()
    report[stage]["manifest"] = {}
    with pytest.raises(ValueError):
        validate(report, scheduler)


def checkpoint(tmp_path):
    directory, executor = tmp_path / "checkpoint", tmp_path / "executor"
    directory.mkdir()
    for name in module.OUTPUT_NAMES:
        path = directory / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("" if name.endswith("stderr") else "x")
    prefix = [str(module.TOOL / "venv_orthomcl/bin/perl"), "-I" + str(module.TOOL),
              str(executor / "benchmark_tools/run_orthomcl_perl_script.pl")]
    report = {"status": "bpo_checkpoint_content_and_indexes_verified", "search_admitted": False,
              "accuracy_admitted": False, "publication_ready": False, "environment": module.environment(),
              "checked_records": [], "outputs": [module.record(directory / n) for n in sorted(module.OUTPUT_NAMES)],
              "commands": [
                  {"stage": "build", "cwd": str(directory), "argv": prefix + [str(executor / "benchmark_tools/build_orthomcl_bpo_indexes.pl"), str(directory / "all.bpo"), str(directory / "indexes")]},
                  {"stage": "validate", "cwd": str(directory), "argv": prefix + [str(executor / "benchmark_tools/validate_orthomcl_bpo_indexes.pl"), str(directory / "all.bpo"), str(directory / "indexes/all_bpo.idx"), str(directory / "indexes/all_bpo.se")]}]}
    return report, directory, executor


def test_checkpoint_inventory_and_commands(tmp_path):
    module.validate_checkpoint(*checkpoint(tmp_path))


@pytest.mark.parametrize("problem", ["status", "environment", "duplicate", "missing", "path", "command", "cwd", "stderr", "mutation"])
def test_wrong_checkpoint_rejected(tmp_path, problem):
    report, directory, executor = checkpoint(tmp_path)
    if problem in {"status", "environment"}:
        report[problem] = "wrong"
    elif problem == "duplicate":
        report["outputs"].append(report["outputs"][0])
    elif problem == "missing":
        report["outputs"].pop()
    elif problem == "path":
        report["outputs"][0]["path"] += ".partial"
    elif problem == "command":
        report["commands"][0]["argv"].append("--changed")
    elif problem == "cwd":
        report["commands"][1]["cwd"] = str(tmp_path)
    elif problem == "stderr":
        (directory / "build.stderr").write_text("warning")
        report["outputs"] = [module.record(directory / n) for n in sorted(module.OUTPUT_NAMES)]
    else:
        (directory / "all.bpo").write_text("changed")
    with pytest.raises(ValueError):
        module.validate_checkpoint(report, directory, executor)


@pytest.mark.skipif(os.environ.get("ORTHOHMM_LEGACY_BLAST_SMOKE") != "1", reason="Installed native runtime")
def test_independent_complete_recheck(tmp_path):
    from benchmark_tools.prepare_orthomcl_bpo_checkpoint import prepare
    from benchmark_tools.probe_orthomcl_bpo_parity import fixture
    fasta, blast = fixture(tmp_path)
    directory = tmp_path / "checkpoint"
    expected = prepare(Path(__file__).resolve().parents[2], blast, fasta, directory)
    result = module.recheck_content(blast, fasta, directory, tmp_path / "recheck",
                                    expected["content"], expected["index_validation"])
    assert result["content"]["bpo_pair_records"] == 6
    assert result["index_validation"]["queries"] == 3
    with pytest.raises(ValueError, match="content counts"):
        module.recheck_content(blast, fasta, directory, tmp_path / "bad_content", {}, expected["index_validation"])
    with pytest.raises(ValueError, match="index summary"):
        module.recheck_content(blast, fasta, directory, tmp_path / "bad_indexes", expected["content"], {})


@pytest.mark.parametrize("problem", [None, "recheck", "runtime_after"])
def test_admission_orchestration_retains_terminal_result(tmp_path, monkeypatch, problem):
    root = tmp_path
    base = root / "benchmarks/results/qfo_corrected_orthomcl_v1/bpo_preparation"
    (base / "checkpoint").mkdir(parents=True)
    fake_record = lambda path: {"path": str(path), "bytes": 1, "sha256": "x"}
    monkeypatch.setattr(module, "record", fake_record)
    monkeypatch.setattr(module, "check", lambda item: None)
    monkeypatch.setattr(module, "validate_checkpoint", lambda *args: None)
    monkeypatch.setattr(module, "read_frozen", lambda *args: {})
    monkeypatch.setattr(module, "verify", lambda value: None)
    count = 0
    def runtime(root):
        nonlocal count
        count += 1
        if problem == "runtime_after" and count == 2:
            raise ValueError("runtime changed")
        return {"manifest": {"runtime": 1}}
    monkeypatch.setattr(module, "verify_runtime", runtime)
    report, scheduler = preparation()
    checkpoint = {"content": {"input_proteins": 984137}, "index_validation": {},
                  "checked_records": [], "outputs": []}
    (base / "checkpoint/report.json").write_text(json.dumps(checkpoint))
    executor = root / "benchmarks/work/publication_qfo_corrected_bpo_v1"
    admission_path = root / "benchmarks/work/qfo_corrected_blast_admission_20260918/report.json"
    report.update(source=fake_record(executor / "benchmark_tools/prepare_qfo_corrected_bpo.py"),
                  checkpoint=fake_record(base / "checkpoint/report.json"),
                  checked_records=[fake_record(admission_path)], content=checkpoint["content"],
                  index_validation={}, query_coverage={"input_proteins": 984137},
                  admission_scheduler={"JobIDRaw": "2"}, admitted_search_scheduler={"JobIDRaw": "3"})
    (base / "report.json").write_text(json.dumps(report))
    def query(argv, **kwargs):
        if argv[0] == "git":
            return module.EXECUTOR
        return "|".join(scheduler) + "\n" + "|".join(scheduler.values()) + "\n"
    monkeypatch.setattr(module.subprocess, "check_output", query)
    monkeypatch.setattr(module.subprocess, "run", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "verify_admission", lambda *args: (
        {"scheduler": report["admitted_search_scheduler"], "query_coverage": report["query_coverage"]},
        [fake_record(root / "blast"), fake_record(root / "fasta")], [], report["admission_scheduler"], ""))
    def recheck(*args):
        if problem == "recheck":
            raise ValueError("recheck failed")
        return {"outputs": [], "content": checkpoint["content"], "index_validation": {}}
    monkeypatch.setattr(module, "recheck_content", recheck)
    output = root / "admission"
    if problem:
        with pytest.raises(ValueError):
            module.admit(root, 1, output)
        result = json.loads((output / "report.json").read_text())
        assert result["status"] == "failed"
    else:
        result = module.admit(root, 1, output)
        assert result["status"] == "corrected_orthomcl_bpo_checkpoint_admitted"
        assert len(result["native_inputs"]) == 3
        with pytest.raises(FileExistsError):
            module.admit(root, 1, output)
    assert result["accuracy_admitted"] is result["publication_ready"] is False
