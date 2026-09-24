import json
from types import SimpleNamespace

import pytest

import benchmark_tools.prepare_blast_recovery_bpo as module


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
