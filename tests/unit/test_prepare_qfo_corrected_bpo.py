import copy
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import prepare_qfo_corrected_bpo as module


def admission(root):
    work = root / "benchmarks/results/qfo_corrected_orthomcl_v1/work"
    return {"status": "corrected_orthomcl_search_evidence_verified", "search_admitted": True,
            "accuracy_admitted": False, "publication_ready": False,
            "downstream_execution_authorized": False,
            "scheduler": {"State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon",
                          "AllocCPUS": "180", "ReqMem": "900G"},
            "query_coverage": {"input_proteins": 984137, "logged_failed_queries": ["retained"],
                               "hsp_rows": 8, "distinct_directed_pairs": 7},
            "database_content": {"input_sequences": 984137},
            "checked_records": [{"path": str(work / name), "bytes": 10, "sha256": "x"}
                                for name in ("all.blast", "all.fa")]}


def test_admitted_inputs_preserve_failure_diagnostics(tmp_path):
    value = admission(tmp_path)
    before = copy.deepcopy(value)
    assert module.validate_admission(value, tmp_path) == value["checked_records"]
    assert value == before
    value["checked_records"].append(copy.deepcopy(value["checked_records"][1]))
    assert len(module.validate_admission(value, tmp_path)) == 2


@pytest.mark.parametrize("key,value", [("status", "running"), ("search_admitted", False),
    ("search_admitted", 1), ("accuracy_admitted", True), ("publication_ready", True),
    ("downstream_execution_authorized", True)])
def test_invalid_admission_flags(tmp_path, key, value):
    data = admission(tmp_path)
    data[key] = value
    with pytest.raises(ValueError):
        module.validate_admission(data, tmp_path)


@pytest.mark.parametrize("key,value", [("State", "RUNNING"), ("ExitCode", "1:0"),
    ("NodeList", "elsewhere"), ("AllocCPUS", "2"), ("ReqMem", "64G")])
def test_invalid_native_scheduler(tmp_path, key, value):
    data = admission(tmp_path)
    data["scheduler"][key] = value
    with pytest.raises(ValueError):
        module.validate_admission(data, tmp_path)


@pytest.mark.parametrize("problem", ["table_scope", "database_scope", "missing", "conflicting",
                                      "empty", "other_path"])
def test_wrong_input_universe_and_records(tmp_path, problem):
    data = admission(tmp_path)
    if problem == "table_scope":
        data["query_coverage"]["input_proteins"] -= 1
    elif problem == "database_scope":
        data["database_content"]["input_sequences"] -= 1
    elif problem == "missing":
        data["checked_records"].pop()
    elif problem == "conflicting":
        data["checked_records"].append({**data["checked_records"][0], "sha256": "y"})
    elif problem == "empty":
        data["checked_records"][0]["bytes"] = 0
    else:
        data["checked_records"][0]["path"] += ".partial"
    with pytest.raises(ValueError):
        module.validate_admission(data, tmp_path)


@pytest.mark.parametrize("problem", [None, "failure", "wrong_scope", "mutation", "counts", "status"])
def test_preparation_preserves_success_or_failure(tmp_path, monkeypatch, problem):
    monkeypatch.setattr(module, "verify_runtime", lambda root: {"status": "mock_verified"})
    for key, value in {"SLURM_JOB_ID": "123", "SLURM_CPUS_PER_TASK": "2", "SLURM_MEM_PER_NODE": "65536"}.items():
        monkeypatch.setenv(key, value)
    monkeypatch.setattr(module.os, "uname", lambda: SimpleNamespace(nodename="bizon"))
    source = tmp_path / "source"
    source.write_text("original")
    data = admission(tmp_path)
    checked = [module.record(source)]
    monkeypatch.setattr(module, "verify_admission", lambda *args: (data, data["checked_records"], checked, {}, ""))
    def checkpoint(root, blast, fasta, output):
        output.mkdir()
        if problem == "failure":
            raise ValueError("test failure")
        if problem == "mutation":
            source.write_text("changed")
        (output / "report.json").write_text("{}")
        return {"status": "failed" if problem == "status" else "bpo_checkpoint_content_and_indexes_verified",
                "content": {"input_proteins": 984136 if problem == "wrong_scope" else 984137,
                            "source_hsp_rows": 9 if problem == "counts" else 8, "source_pair_blocks": 7},
                "index_validation": {"records": 6}, "checked_records": [], "outputs": []}
    monkeypatch.setattr(module, "checkpoint", checkpoint)
    output = tmp_path / "benchmarks/results/qfo_corrected_orthomcl_v1/bpo_preparation"
    if problem:
        with pytest.raises(ValueError):
            module.prepare(tmp_path, tmp_path / "admission", "sha", 99)
        result = json.loads((output / "report.json").read_text())
        assert result["status"] == "failed"
    else:
        result = module.prepare(tmp_path, tmp_path / "admission", "sha", 99)
        assert result["status"] == "corrected_bpo_checkpoint_prepared_pending_admission"
        assert result["query_coverage"] == data["query_coverage"]
        with pytest.raises(FileExistsError):
            module.prepare(tmp_path, tmp_path / "admission", "sha", 99)
    assert result["accuracy_admitted"] is result["publication_ready"] is False


def test_unscheduled_execution_is_rejected(tmp_path, monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    with pytest.raises(ValueError, match="scheduled"):
        module.prepare(tmp_path, tmp_path / "admission", "sha", 99)


@pytest.mark.parametrize("problem", [None, "commit", "source", "bytes", "job", "allocation", "sha"])
def test_admission_provenance_binding(tmp_path, monkeypatch, problem):
    data = admission(tmp_path)
    records = []
    for item in data["checked_records"]:
        path = Path(item["path"])
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("data")
        records.append(module.record(path))
    data["checked_records"] = records
    source = tmp_path / "benchmarks/work/publication_qfo_corrected_blast_admission_v1/benchmark_tools/admit_qfo_corrected_blast.py"
    source.parent.mkdir(parents=True)
    source.write_text("source")
    data["source"] = module.record(source)
    for name in ("database_audit", "table_audit"):
        path = tmp_path / name
        path.write_text("{}")
        data[name] = module.record(path)
    if problem == "source":
        data["source"]["sha256"] = "bad"
    path = tmp_path / "admission.json"
    path.write_text(json.dumps(data))
    sha = module.record(path)["sha256"]
    if problem == "bytes":
        Path(records[0]["path"]).write_text("changed")
    def query(argv, **kwargs):
        if argv[0] == "git":
            return "bad" if problem == "commit" else module.ADMITTER
        return ("JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS|ReqMem\n"
                + ("100" if problem == "job" else "99") + "|COMPLETED|0:0|00:01|bizon|"
                + ("3" if problem == "allocation" else "2") + "|64G\n")
    monkeypatch.setattr(module.subprocess, "check_output", query)
    monkeypatch.setattr(module.subprocess, "run", lambda *args, **kwargs: None)
    if problem:
        with pytest.raises((ValueError, RuntimeError)):
            module.verify_admission(tmp_path, path, "bad" if problem == "sha" else sha, 99)
    else:
        result, inputs, checked, scheduler, _ = module.verify_admission(tmp_path, path, sha, 99)
        assert result == data and inputs == records and scheduler["JobIDRaw"] == "99"
        for item in checked:
            module.check(item)
