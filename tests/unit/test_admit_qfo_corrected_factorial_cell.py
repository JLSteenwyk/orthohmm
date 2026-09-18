from copy import deepcopy
from types import SimpleNamespace

import pytest

import benchmark_tools.admit_qfo_corrected_factorial_cell as module


def fixture():
    scheduler = {"JobIDRaw": "123", "State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon", "AllocCPUS": "32"}
    status = {"status": "finished_pending_native_validation", "failed_methods": [],
              "accuracy_evaluated": False, "native_outputs_validated": False,
              "provenance": {"slurm_job_id": "123", "slurm_array_task_id": "0", "cell_index": 0}}
    post = {"status": "corrected_inputs_runtime_sources_reverified", "failed_methods": [],
            "native_outputs_validated": False, "accuracy_evaluated": False}
    return scheduler, status, post


def test_success_requires_corrected_postflight():
    scheduler, status, post = fixture()
    module.require_success(scheduler, status, post, 0)
    post["status"] = "inputs_runtime_sources_reverified"
    with pytest.raises(ValueError):
        module.require_success(scheduler, status, post, 0)


@pytest.mark.parametrize("key,value", [("State", "RUNNING"), ("ExitCode", "1:0"),
    ("NodeList", "spark-7ff0"), ("AllocCPUS", "20"), ("JobIDRaw", "999")])
def test_bad_scheduler(key, value):
    scheduler, status, post = fixture()
    scheduler[key] = value
    with pytest.raises(ValueError):
        module.require_success(scheduler, status, post, 0)


@pytest.mark.parametrize("key,value", [("status", "running"), ("failed_methods", ["cell"]),
    ("accuracy_evaluated", True), ("native_outputs_validated", True)])
def test_unfinished_or_scored_status(key, value):
    scheduler, status, post = fixture()
    status[key] = value
    with pytest.raises(ValueError):
        module.require_success(scheduler, status, post, 0)


@pytest.mark.parametrize("key,value", [("slurm_job_id", "999"), ("slurm_array_task_id", "1"),
    ("cell_index", 1)])
def test_wrong_task(key, value):
    scheduler, status, post = fixture()
    status["provenance"][key] = value
    with pytest.raises(ValueError):
        module.require_success(scheduler, status, post, 0)


def test_provenance_exact_inventory():
    expected = {"candidate_admission": {"sha256": "corrected"}, "executed_argv": ["frozen"], "cwd": "/launcher"}
    manifest = {"input_fastas": [{"path": "/corrected/a.fa", "sha256": "input", "bytes": 1}]}
    cell = {"label": "p0_c0_r1"}
    status = {"provenance": expected, "verified_inputs": {"status": "ready", "inputs": [
        {**manifest["input_fastas"][0], "absolute_path": "/corrected/a.fa"}]},
        "methods": {cell["label"]: {}}, "dataset": cell["label"]}
    module.validate_provenance(status, expected, manifest, cell)
    for key, value in [("provenance", {**expected, "cwd": "/wrong"}),
                       ("verified_inputs", {"status": "ready", "inputs": []}),
                       ("methods", {"wrong": {}}), ("dataset", "wrong")]:
        changed = deepcopy(status)
        changed[key] = value
        with pytest.raises(ValueError):
            module.validate_provenance(changed, expected, manifest, cell)


@pytest.mark.parametrize("state,exit_code", [("RUNNING", "0:0"), ("FAILED", "1:0"), ("CANCELLED", "0:0")])
def test_terminal_gate_precedes_any_output_access(monkeypatch, tmp_path, state, exit_code):
    accounting = f"JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n123|{state}|{exit_code}|01:00:00|bizon|32\n"
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: accounting)
    monkeypatch.setattr(module, "verify_admission", lambda *a: pytest.fail("Must not read partial outputs"))
    with pytest.raises(ValueError):
        module.admit(tmp_path, 0, "123", tmp_path / "missing", "0" * 64, "122")


@pytest.mark.parametrize("index", [-1, 4, True, "0"])
def test_invalid_index_rejected_without_io(tmp_path, index):
    with pytest.raises(ValueError):
        module.admit(tmp_path, index, "123", tmp_path / "missing", "0" * 64, "122")


def test_full_corrected_gene_universe(monkeypatch, tmp_path):
    native = {"input_proteomes": [{"filename": f"sp{i}.fa", "taxon": f"sp{i}"} for i in range(78)]}
    manifest = {"input_fastas": [{"path": f"/input/sp{i}.fa"} for i in range(78)]}
    def parse(path, fmt):
        i = int(path.split("sp")[1].split(".")[0])
        return (SimpleNamespace(id=f"g{j}") for j in range(i, 984137, 78))
    monkeypatch.setattr(module.SeqIO, "parse", parse)
    candidate = tmp_path / "partition.txt"
    candidate.write_text(" ".join(f"g{i}" for i in range(984137)) + "\n")
    owners, candidates = module.gene_ownership(manifest, native, candidate)
    assert len(owners) == len(candidates) == 984137
    assert owners["g79"] == "sp1"
    candidate.write_text("g0 g1\n")
    with pytest.raises(ValueError, match="universe"):
        module.gene_ownership(manifest, native, candidate)


def test_duplicate_taxa_and_candidate_members(monkeypatch, tmp_path):
    native = {"input_proteomes": [{"filename": f"sp{i}.fa", "taxon": f"sp{i}"} for i in range(78)]}
    manifest = {"input_fastas": [{"path": "/input/sp0.fa"}]}
    monkeypatch.setattr(module.SeqIO, "parse", lambda *a: [SimpleNamespace(id="g0")])
    candidate = tmp_path / "partition.txt"
    candidate.write_text("g0 g0\n")
    with pytest.raises(ValueError, match="Duplicate candidate"):
        module.gene_ownership(manifest, native, candidate)
    native["input_proteomes"][1]["taxon"] = "sp0"
    with pytest.raises(ValueError, match="distinct"):
        module.gene_ownership(manifest, native, candidate)
