import json

import pytest

from benchmark_tools import summarize_blast_recovery_failures as module


def fixture(tmp_path):
    fasta, log = tmp_path / "queries_00.fa", tmp_path / "blast.log"
    fasta.write_text(">failed\nAAA\n>hit\nACDEFG\n>nohit\nGGGG\n")
    log.write_text("[blastall] WARNING: [blast] failed: SetUpBlastSearch failed.\n"
                   "[blastall] ERROR: [blast] failed: Blast: Query must be at least twice wordsize for two hit mode\n")
    accounting = "JobID|JobIDRaw|State|ExitCode|NodeList|AllocCPUS|Elapsed\n22103_0|22118|COMPLETED|0:0|bizon|180|00:01:00\n"
    genes, blocks = ["failed", "hit", "nohit"], [{"query": "hit"}]
    diagnostics = module.parse_diagnostics(log)
    report = dict(status="recovery_batch_execution_and_rows_verified", index=0, batch_admitted=True,
        executor_commit=module.EXECUTOR_COMMIT, scheduler=module.completed_task(accounting, 0),
        records=[module.record(fasta), module.record(log)], query_ids=genes, query_blocks=blocks,
        diagnostics=diagnostics, coverage=module.coverage(genes, blocks, diagnostics))
    path = tmp_path / "admission.json"
    path.write_text(json.dumps(report))
    return path, report, accounting


def test_partial_failures_not_nohits_or_whole_search(tmp_path):
    path, _, accounting = fixture(tmp_path)
    result = module.summarize([path], accounting)
    assert result["input_queries"] == 3
    assert result["failed_queries"] == result["queries_with_hits"] == result["no_hits_without_logged_failure"] == 1
    assert result["failed_query_lengths"] == dict(min=3, median=3, max=3)
    assert result["failed_queries_by_category"] == dict(setup_failure=1, short_query_failure=1)
    assert result["records"][0]["residue_counts"] == {"A": 3}
    assert all(result[key] is False for key in ("search_admitted", "accuracy_admitted", "publication_ready"))


@pytest.mark.parametrize("problem", ["status", "commit", "duplicate", "scheduler", "fasta", "log", "coverage", "genes", "record", "outgoing"])
def test_invalid_or_changed_admission_rejected(tmp_path, problem):
    path, report, accounting = fixture(tmp_path)
    paths = [path]
    if problem == "status":
        report["batch_admitted"] = False
    elif problem == "commit":
        report["executor_commit"] = "other"
    elif problem == "duplicate":
        paths *= 2
    elif problem == "scheduler":
        accounting = accounting.replace("COMPLETED", "FAILED")
    elif problem in ("fasta", "log"):
        (tmp_path / ("queries_00.fa" if problem == "fasta" else "blast.log")).write_text("changed")
    elif problem == "coverage":
        report["coverage"]["input_queries"] = 4
    elif problem == "genes":
        report["query_ids"].reverse()
    elif problem == "record":
        report["records"].pop()
    else:
        report["query_blocks"].append({"query": "failed"})
        report["coverage"] = module.coverage(report["query_ids"], report["query_blocks"], report["diagnostics"])
    path.write_text(json.dumps(report))
    with pytest.raises(ValueError):
        module.summarize(paths, accounting)


def test_empty_input_rejected():
    with pytest.raises(ValueError):
        module.summarize([], "")
