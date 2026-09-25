from copy import deepcopy
import json
from pathlib import Path

import pytest

import benchmark_tools.admit_blast_recovery_search as module


def accounting():
    return "JobID|State|ExitCode|NodeList|AllocCPUS|ReqMem|Elapsed\n22150|COMPLETED|0:0|bizon|2|64G|00:10:00\n"


def test_completed_merge_allocation():
    assert module.completed_merge(accounting())["JobID"] == "22150"


@pytest.mark.parametrize("old,new", [("COMPLETED", "RUNNING"), ("0:0", "1:0"),
    ("bizon", "other"), ("|2|", "|1|"), ("64G", "8G"), ("22150", "21713")])
def test_wrong_merge_execution(old, new):
    with pytest.raises(ValueError):
        module.completed_merge(accounting().replace(old, new))


def test_duplicate_scheduler_record():
    text = accounting()
    with pytest.raises(ValueError):
        module.completed_merge(text + text.splitlines()[1] + "\n")


def status(directory):
    source = dict(path="frozen/driver.py", bytes=1, sha256="fixture")
    candidate = dict(status="merged_candidate_requires_full_admission",
        path=str(directory / "table/all.blast.candidate"), bytes=10, sha256="fixture",
        search_admitted=False, reuse_authorized=False, publication_ready=False)
    return dict(status="merged_candidate_pending_full_table_admission", job_id="22150",
        source=source, candidate=candidate, selected_log={"path": str(directory / "selected.blast.log")},
        search_admitted=False, reuse_authorized=False, publication_ready=False)


def test_merge_output_binding():
    directory = Path("/fixture")
    report = status(directory)
    assert module.validate_merge(report, report["source"], directory)["bytes"] == 10


@pytest.mark.parametrize("change", ["job", "source", "status", "admitted", "candidate_path",
    "candidate_status", "candidate_admitted", "log"])
def test_wrong_merge_evidence(change):
    directory = Path("/fixture")
    report = status(directory)
    source = deepcopy(report["source"])
    if change == "job":
        report["job_id"] = "21713"
    elif change == "source":
        report["source"]["sha256"] = "changed"
    elif change == "status":
        report["status"] = "merge_failed_preserved"
    elif change == "admitted":
        report["search_admitted"] = True
    elif change == "candidate_path":
        report["candidate"]["path"] = "/other/all.blast"
    elif change == "candidate_status":
        report["candidate"]["status"] = "partial"
    elif change == "candidate_admitted":
        report["candidate"]["reuse_authorized"] = True
    elif change == "log":
        report["selected_log"]["path"] = "/other/blast.log"
    with pytest.raises(ValueError):
        module.validate_merge(report, source, directory)


def test_pending_merge_rejected_before_artifacts(tmp_path, monkeypatch):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
        accounting().replace("COMPLETED", "PENDING"))
    with pytest.raises(ValueError, match="completed frozen"):
        module.admit(tmp_path, tmp_path / "admission")
    assert not (tmp_path / "admission").exists()


def test_existing_output_preserved(tmp_path):
    output = tmp_path / "admission"
    output.mkdir()
    with pytest.raises(FileExistsError):
        module.admit(tmp_path, output)


def orchestration(tmp_path, monkeypatch, failure=None, replacement=False):
    executor = tmp_path / ("benchmarks/work/blast_replacement_merge_v1_20260925/benchmark_tools" if replacement
                           else "benchmarks/work/blast_recovery_merge_v1_20260923/benchmark_tools")
    executor.mkdir(parents=True)
    local = Path(module.__file__).with_name("run_blast_recovery_merge.py")
    (executor / local.name).write_bytes(local.read_bytes())
    directory = tmp_path / ("benchmarks/results/qfo_blast_replacement_merge_v1" if replacement
                            else "benchmarks/results/qfo_blast_recovery_merge_v1")
    (directory / "table").mkdir(parents=True)
    candidate = directory / "table/all.blast.candidate"
    candidate.write_text("fixture candidate\n")
    log = directory / "selected.blast.log"
    log.write_text("")
    report = status(directory)
    if replacement:
        report["job_id"] = module.REPLACEMENT_MERGE_JOB
    report["source"] = module.record(executor / local.name)
    report["candidate"].update(module.record(candidate), rows=3, query_blocks=2)
    report.update(selected_log=module.record(log), checked_inputs=[], diagnostics={}, replay_coverage={})
    (directory / "status.json").write_text(json.dumps(report))
    results = tmp_path / "benchmark_tools/results"
    results.mkdir(parents=True)
    for name in ("qfo_corrected_orthomcl_prepared_20260918.json", "qfo_corrected_legacy_blast_runtime_20260918.json"):
        (results / name).write_text("{}\n")
    def subprocess_output(command, **kwargs):
        if command[0] == "sacct":
            return accounting().replace(module.MERGE_JOB, module.REPLACEMENT_MERGE_JOB) if replacement else accounting()
        return (module.REPLACEMENT_MERGE_COMMIT if replacement else module.MERGE_COMMIT) + "\n"
    monkeypatch.setattr(module.subprocess, "check_output", subprocess_output)
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k: None)
    calls = []
    def prepare(root, output, replacement=False):
        calls.append("prepare")
        if failure == "prepare":
            raise ValueError("Injected prerequisite failure")
        return dict(selected_diagnostics={}, panel={"coverage": {}, "admissions": []},
            prefix_audit={"blocks": {"path": "fixture-prefix"}}, records=[])
    def verify(*args):
        calls.append("verify")
        return {"output_root": str(tmp_path / ("changed" if failure == "plan" and calls.count("verify") > 1 else "native"))}
    def database(fasta, runtime, output):
        calls.append("database")
        output.mkdir()
        (output / "report.json").write_text("{}\n")
        return dict(status="database_exact_sequence_parity_verified", checked_records=[], outputs=[],
            content={"input_sequences": 1 if failure == "database" else 984137})
    def table(blast, fasta, log_path, prefix, batches, output):
        calls.append("table")
        output.write_text("{}\n")
        if failure == "table":
            raise ValueError("Injected table failure")
        if failure == "mutation":
            candidate.write_text("changed after initial validation\n")
        return dict(content={"input_proteins": 984137, "hsp_rows": 4 if failure == "rows" else 3,
            "diagnostics": [{"gene": "x", "query_failed": False, "messages": []}] if failure == "diagnostics" else []},
            query_blocks=2, checked_records=[])
    monkeypatch.setattr(module, "prepare", prepare)
    monkeypatch.setattr(module, "verify", verify)
    monkeypatch.setattr(module, "audit_database", database)
    monkeypatch.setattr(module, "audit_candidate", table)
    return calls


def test_successful_orchestration_admits_search_not_downstream(tmp_path, monkeypatch):
    calls = orchestration(tmp_path, monkeypatch)
    result = module.admit(tmp_path, tmp_path / "admission")
    assert calls == ["prepare", "verify", "database", "table", "verify"]
    assert result["status"] == "recovered_orthomcl_search_evidence_verified"
    assert result["search_admitted"] is True
    assert result["accuracy_admitted"] is result["publication_ready"] is result["downstream_execution_authorized"] is False
    assert json.loads((tmp_path / "admission/report.json").read_text()) == result


@pytest.mark.parametrize("failure", ["prepare", "database", "table", "rows", "diagnostics", "plan", "mutation"])
def test_orchestration_failure_never_admits(tmp_path, monkeypatch, failure):
    calls = orchestration(tmp_path, monkeypatch, failure)
    with pytest.raises(ValueError):
        module.admit(tmp_path, tmp_path / "admission")
    result = json.loads((tmp_path / "admission/report.json").read_text())
    assert result["status"] == "recovery_search_admission_failed"
    assert result["error_type"] == "ValueError"
    assert result["search_admitted"] is result["downstream_execution_authorized"] is False
    if failure in {"prepare", "database"}:
        assert "table" not in calls


@pytest.mark.parametrize("failure", [None, "representation", "table", "rows", "diagnostics", "mutation"])
def test_reviewed_native_mode_still_requires_whole_table(tmp_path, monkeypatch, failure):
    from benchmark_tools import reviewed_legacy_database
    calls = orchestration(tmp_path, monkeypatch, failure, replacement=True)
    def representation(root, database):
        calls.append("representation")
        if failure == "representation":
            raise ValueError("Unexpected transformation")
        return dict(exact_sequence_parity=False, checked_records=[], limitations=["Native O deletion retained"])
    monkeypatch.setattr(reviewed_legacy_database, "verify", representation)
    if failure:
        with pytest.raises(ValueError):
            module.admit(tmp_path, tmp_path / "admission", True, True)
        result = json.loads((tmp_path / "admission/report.json").read_text())
        assert result["search_admitted"] is False
    else:
        result = module.admit(tmp_path, tmp_path / "admission", True, True)
        assert result["status"] == "recovered_search_native_representation_verified"
        assert result["database_representation"]["exact_sequence_parity"] is False
        assert result["search_admitted"] is True
        assert result["downstream_execution_authorized"] is False
        assert "table" in calls
