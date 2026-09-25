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


def replacement_fixture(tmp_path, monkeypatch):
    path, report, _ = fixture(tmp_path)
    fasta = tmp_path / "queries_14.fa"
    (tmp_path / "queries_00.fa").rename(fasta)
    original_dir = tmp_path / "interrupted"
    original_dir.mkdir()
    original_log = original_dir / "blast.log"
    original_log.write_text("Preserved interrupted log, not replacement diagnostics")
    history = dict(records=[module.record(original_log)], partial_rows_reused=False,
                   replacement_task="22160_14", original_scheduler={"State": "TIMEOUT"})
    monkeypatch.setattr(module, "interrupted_attempt", lambda root, accounting: history)
    accounting = "JobID|JobIDRaw|State|ExitCode|NodeList|AllocCPUS|Elapsed\n22160_14|22160|COMPLETED|0:0|bizon|180|00:01:00\n"
    report.update(index=14, replacement=dict(history),
                  scheduler=module.completed_task(accounting, 14, replacement=True),
                  records=[module.record(fasta), module.record(tmp_path / "blast.log"), *history["records"]])
    path.write_text(json.dumps(report))
    return path, report, accounting


def test_replacement_excludes_only_verified_interrupted_log(tmp_path, monkeypatch):
    path, _, accounting = replacement_fixture(tmp_path, monkeypatch)
    result = module.summarize([path], accounting, replacement_root=tmp_path)
    assert result["failed_queries"] == 1
    assert result["batches"][0]["scheduler"]["JobID"] == "22160_14"
    assert sum(r["path"].endswith("blast.log") for r in result["checked_records"]) == 2


@pytest.mark.parametrize("problem", ["root", "history", "missing", "duplicate", "changed", "index", "scheduler"])
def test_invalid_replacement_rejected(tmp_path, monkeypatch, problem):
    path, report, accounting = replacement_fixture(tmp_path, monkeypatch)
    root = tmp_path
    if problem == "root":
        root = None
    elif problem == "history":
        report["replacement"]["partial_rows_reused"] = True
    elif problem == "missing":
        report["records"].pop()
    elif problem == "duplicate":
        report["records"].append(report["records"][-1])
    elif problem == "changed":
        (tmp_path / "interrupted/blast.log").write_text("Changed")
    elif problem == "index":
        report["index"] = 0
    else:
        accounting = accounting.replace("COMPLETED", "FAILED")
    path.write_text(json.dumps(report))
    with pytest.raises(ValueError):
        module.summarize([path], accounting, replacement_root=root)
