import copy

import pytest

from benchmark_tools.admit_blast_recovery_batch import coverage, EXECUTOR_COMMIT
from benchmark_tools.verify_blast_recovery_panel import combine, unique_records, verify


def panel():
    reports = []
    expected = [["a", "b", "c"], ["d", "e"]]
    for index, genes in enumerate(expected):
        blocks = [dict(query=genes[0], rows=2)]
        diagnostics = {"b": dict(query_failed=True)} if index == 0 else {}
        reports.append(dict(index=index, status="recovery_batch_execution_and_rows_verified",
            executor_commit=EXECUTOR_COMMIT, batch_admitted=True, search_admitted=False,
            reuse_authorized=False, publication_ready=False, query_ids=genes,
            query_blocks=blocks, diagnostics=diagnostics, coverage=coverage(genes, blocks, diagnostics)))
    return reports, expected, ["a", "b", "c", "d", "e"]


def test_complete_panel_retains_failed_and_unexplained_no_hits():
    result = combine(*panel())
    assert result["replay_queries"] == 5
    assert result["hsp_rows"] == 4
    assert result["failed_queries"] == ["b"]
    assert result["no_hits_without_logged_failure"] == ["c", "e"]


@pytest.mark.parametrize("problem", ["missing", "order", "duplicate", "index", "bool_index",
    "unadmitted", "promoted", "coverage", "blocks", "foreign", "failed_hits", "executor"])
def test_reject_incomplete_or_changed_panel(problem):
    reports, expected, replay = copy.deepcopy(panel())
    if problem == "missing":
        reports.pop()
    elif problem == "order":
        replay.reverse()
    elif problem == "duplicate":
        expected[1][1] = "a"
        reports[1]["query_ids"] = expected[1]
        reports[1]["coverage"] = coverage(expected[1], reports[1]["query_blocks"], {})
        replay[-1] = "a"
    elif problem == "index":
        reports[1]["index"] = 0
    elif problem == "bool_index":
        reports[0]["index"] = False
    elif problem == "unadmitted":
        reports[0]["batch_admitted"] = False
    elif problem == "promoted":
        reports[0]["search_admitted"] = True
    elif problem == "coverage":
        reports[0]["coverage"]["input_queries"] = 6
    elif problem == "blocks":
        reports[0]["query_blocks"] *= 2
    elif problem == "foreign":
        reports[0]["diagnostics"]["z"] = dict(query_failed=True)
    elif problem == "failed_hits":
        reports[0]["query_blocks"].append(dict(query="b", rows=1))
        reports[0]["coverage"] = coverage(expected[0], reports[0]["query_blocks"], reports[0]["diagnostics"])
    elif problem == "executor":
        reports[0]["executor_commit"] = "changed"
    with pytest.raises(ValueError):
        combine(reports, expected, replay)


def test_provenance_deduplication_and_conflict():
    item = dict(path="/x", sha256="a", bytes=1)
    assert unique_records([item, dict(item)]) == [item]
    with pytest.raises(ValueError):
        unique_records([item, dict(item, sha256="b")])


def test_no_overwrite(tmp_path):
    with pytest.raises(FileExistsError):
        verify(tmp_path, tmp_path)
