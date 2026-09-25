from pathlib import Path

import pytest

from benchmark_tools.prepare_blast_recovery_bpo import admission_contract, validate_admission


@pytest.mark.parametrize("replacement", [False, True])
def test_exact_admission_path_selects_frozen_contract(replacement):
    root = Path("/project")
    name = "qfo_blast_replacement_search_admission_v1" if replacement else "qfo_blast_recovery_search_admission_v1"
    mode, job, commit, executor = admission_contract(root, root / f"benchmarks/results/{name}/report.json")
    assert mode is replacement
    assert job == ("22163" if replacement else "22151")
    assert len(commit) == 40
    assert executor.is_relative_to(root / "benchmarks/work")
    with pytest.raises(ValueError):
        admission_contract(root, root / f"benchmarks/results/{name}/other.json")


@pytest.mark.parametrize("replacement", [False, True])
@pytest.mark.parametrize("problem", [None, "merge", "path", "parity", "universe", "promoted", "conflict", "candidate"])
def test_admitted_inputs_cannot_mix_attempts(replacement, problem):
    root = Path("/project")
    merge = "qfo_blast_replacement_merge_v1" if replacement else "qfo_blast_recovery_merge_v1"
    candidate = dict(path=str(root / f"benchmarks/results/{merge}/table/all.blast.candidate"), bytes=100, sha256="table")
    fasta = dict(path=str(root / "benchmarks/results/qfo_corrected_orthomcl_v1/work/all.fa"), bytes=100, sha256="fasta")
    report = dict(status="recovered_orthomcl_search_evidence_verified", search_admitted=True,
                  accuracy_admitted=False, publication_ready=False, downstream_execution_authorized=False,
                  scheduler=dict(JobID="22162" if replacement else "22150", State="COMPLETED", ExitCode="0:0",
                                 NodeList="bizon", AllocCPUS="2", ReqMem="64G"),
                  query_coverage=dict(input_proteins=984137),
                  database_content=dict(input_sequences=984137, exact_sequence_parity=True),
                  checked_records=[candidate, fasta], candidate=candidate)
    if problem == "merge":
        report["scheduler"]["JobID"] = "22150" if replacement else "22162"
    elif problem == "path":
        report["checked_records"][0] = dict(candidate, path="/other")
    elif problem == "parity":
        report["database_content"]["exact_sequence_parity"] = False
    elif problem == "universe":
        report["query_coverage"]["input_proteins"] = 5000
    elif problem == "promoted":
        report["accuracy_admitted"] = True
    elif problem == "conflict":
        report["checked_records"].append(dict(candidate, sha256="changed"))
    elif problem == "candidate":
        report["candidate"] = dict(candidate, sha256="changed")
    if problem:
        with pytest.raises(ValueError):
            validate_admission(report, root, replacement)
    else:
        assert validate_admission(report, root, replacement) == [candidate, fasta]
        with pytest.raises(ValueError):
            validate_admission(report, root, not replacement)
