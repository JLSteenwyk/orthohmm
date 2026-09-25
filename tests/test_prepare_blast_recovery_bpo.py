import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest

import benchmark_tools.prepare_blast_recovery_bpo as module


def test_isolated_script_entrypoint(tmp_path):
    result = subprocess.run([sys.executable, "-I", "-B", str(Path(module.__file__).resolve()), "--help"],
        cwd=tmp_path, capture_output=True, text=True, check=True)
    assert "--admission-sha256" in result.stdout


@pytest.mark.parametrize("problem", ["pending", "duplicate", "memory", "executor", "source"])
def test_provenance_gate(tmp_path, monkeypatch, problem):
    report = admission(tmp_path)
    path = tmp_path / "benchmarks/results/qfo_blast_recovery_search_admission_v1/report.json"
    path.parent.mkdir(parents=True)
    executor = tmp_path / "benchmarks/work/blast_recovery_search_admission_v1_20260923/benchmark_tools"
    executor.mkdir(parents=True)
    source = executor / "admit_blast_recovery_search.py"
    source.write_text("fixture source\n")
    report["checked_records"].append(module.record(source))
    if problem == "source":
        report["checked_records"].pop()
    path.write_text(json.dumps(report))
    row = "22151|COMPLETED|0:0|bizon|2|64G|00:10:00\n"
    if problem == "pending":
        row = row.replace("COMPLETED", "PENDING")
    if problem == "memory":
        row = row.replace("64G", "8G")
    if problem == "duplicate":
        row *= 2
    text = "JobID|State|ExitCode|NodeList|AllocCPUS|ReqMem|Elapsed\n" + row
    monkeypatch.setattr(module.subprocess, "check_output", lambda command, **kw:
        text if command[0] == "sacct" else "changed" if problem == "executor" else module.ADMITTER)
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **kw: None)
    with pytest.raises(ValueError):
        module.verify_admission(tmp_path, path, module.record(path)["sha256"])


def admission(root):
    inputs = [dict(path=str(root / name), bytes=1, sha256="fixture") for name in (
        "benchmarks/results/qfo_blast_recovery_merge_v1/table/all.blast.candidate",
        "benchmarks/results/qfo_corrected_orthomcl_v1/work/all.fa")]
    return dict(status="recovered_orthomcl_search_evidence_verified", search_admitted=True,
        accuracy_admitted=False, publication_ready=False, downstream_execution_authorized=False,
        scheduler=dict(JobID="22150", State="COMPLETED", ExitCode="0:0", NodeList="bizon", AllocCPUS="2", ReqMem="64G"),
        query_coverage=dict(input_proteins=984137, hsp_rows=3, distinct_directed_pairs=2),
        database_content=dict(input_sequences=984137, exact_sequence_parity=True),
        checked_records=inputs, candidate=dict(inputs[0]))


def test_recovered_admission_inputs(tmp_path):
    report = admission(tmp_path)
    assert module.validate_admission(report, tmp_path) == report["checked_records"]


def test_native_admission_path_identity(tmp_path):
    path = tmp_path / "benchmarks/results/qfo_blast_native_representation_admission_v1/report.json"
    replacement, job, commit, executor = module.admission_contract(tmp_path, path)
    assert replacement is True and job == "22166"
    assert commit == "42e6dcff5170d17c490273ab93229833b91a6b23"
    assert executor.name == "blast_native_representation_admission_v1_20260925"


@pytest.mark.parametrize("problem", [None, "status", "representation", "missing_record", "not_replacement", "universe"])
def test_native_representation_contract(tmp_path, monkeypatch, problem):
    from benchmark_tools import reviewed_legacy_database
    report = admission(tmp_path)
    report["status"] = "recovered_search_native_representation_verified"
    report["scheduler"]["JobID"] = "22162"
    report["database_content"]["exact_sequence_parity"] = False
    report["checked_records"][0]["path"] = str(tmp_path / "benchmarks/results/qfo_blast_replacement_merge_v1/table/all.blast.candidate")
    report["candidate"] = dict(report["checked_records"][0])
    helper = tmp_path / "benchmarks/work/blast_native_representation_admission_v1_20260925/benchmark_tools/reviewed_legacy_database.py"
    helper.parent.mkdir(parents=True)
    helper.write_text("fixture helper")
    expected = dict(status="reviewed_native_representation_verified_not_exact_parity",
                    exact_sequence_parity=False, transformations=["fixture"],
                    limitations=["fixture"], checked_records=[module.record(helper)])
    report["database_representation"] = json.loads(json.dumps(expected))
    report["checked_records"].extend(expected["checked_records"])
    monkeypatch.setattr(reviewed_legacy_database, "verify", lambda *a: expected)
    if problem == "status":
        report["status"] = "recovery_search_admission_failed"
    elif problem == "representation":
        report["database_representation"]["exact_sequence_parity"] = True
    elif problem == "missing_record":
        report["checked_records"].pop()
    elif problem == "universe":
        report["query_coverage"]["input_proteins"] = 1
    if problem:
        with pytest.raises(ValueError):
            module.validate_admission(report, tmp_path, problem != "not_replacement", True)
    else:
        assert len(module.validate_admission(report, tmp_path, True, True)) == 2
        with pytest.raises(ValueError):
            module.validate_admission(report, tmp_path, True)


@pytest.mark.parametrize("change", ["old_contract", "not_admitted", "wrong_job", "wrong_memory",
    "missing_input", "conflict", "wrong_candidate", "universe", "parity", "downstream"])
def test_reject_wrong_admission(tmp_path, change):
    report = admission(tmp_path)
    if change == "old_contract":
        report["status"] = "corrected_orthomcl_search_evidence_verified"
    elif change == "not_admitted":
        report["search_admitted"] = False
    elif change == "wrong_job":
        report["scheduler"]["JobID"] = "21713"
    elif change == "wrong_memory":
        report["scheduler"]["ReqMem"] = "900G"
    elif change == "missing_input":
        report["checked_records"].pop()
    elif change == "conflict":
        report["checked_records"].append(dict(report["checked_records"][0], sha256="changed"))
    elif change == "wrong_candidate":
        report["candidate"]["bytes"] = 2
    elif change == "universe":
        report["query_coverage"]["input_proteins"] = 10
    elif change == "parity":
        report["database_content"]["exact_sequence_parity"] = False
    elif change == "downstream":
        report["downstream_execution_authorized"] = True
    with pytest.raises(ValueError):
        module.validate_admission(report, tmp_path)


@pytest.mark.parametrize("failure", [None, "checkpoint", "count", "runtime"])
def test_preparation_orchestration(tmp_path, monkeypatch, failure):
    report = admission(tmp_path)
    path = tmp_path / "admission.json"
    path.write_text(json.dumps(report))
    (tmp_path / "benchmarks/results").mkdir(parents=True)
    monkeypatch.setenv("SLURM_JOB_ID", "fixture")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
    monkeypatch.setenv("SLURM_MEM_PER_NODE", "65536")
    monkeypatch.setattr(module.os, "uname", lambda: SimpleNamespace(nodename="bizon"))
    monkeypatch.setattr(module, "verify_admission", lambda *a:
        (report, report["checked_records"], [module.record(path)], {}, "fixture"))
    calls = []
    def runtime(root):
        calls.append("runtime")
        if failure == "runtime" and len(calls) == 2:
            raise ValueError("Changed runtime")
        return {"verified": True}
    monkeypatch.setattr(module, "verify_runtime", runtime)
    def checkpoint(root, blast, fasta, output):
        assert str(blast).endswith("all.blast.candidate")
        output.mkdir()
        (output / "report.json").write_text("{}\n")
        if failure == "checkpoint":
            raise ValueError("Injected checkpoint failure")
        return dict(status="bpo_checkpoint_content_and_indexes_verified",
            content=dict(input_proteins=984137, source_hsp_rows=4 if failure == "count" else 3, source_pair_blocks=2),
            checked_records=[], outputs=[], index_validation={"fixture": True})
    monkeypatch.setattr(module, "checkpoint", checkpoint)
    if failure:
        with pytest.raises(ValueError):
            module.prepare(tmp_path, path, "fixture")
    else:
        result = module.prepare(tmp_path, path, "fixture")
        assert result["status"] == "recovered_bpo_prepared_pending_admission"
    saved = json.loads((tmp_path / "benchmarks/results/qfo_blast_recovery_bpo_v1/report.json").read_text())
    assert saved["downstream_execution_authorized"] is saved["publication_ready"] is False
    if failure:
        assert saved["status"] == "recovered_bpo_preparation_failed"
