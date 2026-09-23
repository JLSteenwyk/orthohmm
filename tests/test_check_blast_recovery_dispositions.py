from copy import deepcopy

import pytest

from benchmark_tools.check_blast_recovery_dispositions import reconcile


def fixture():
    batch = dict(query_ids=["b", "c", "d"], query_blocks=[dict(query="b", rows=2)],
        coverage=dict(queries_without_hits=["c", "d"], failed_queries=["c"]),
        diagnostics={"c": dict(query_failed=True)})
    content = dict(input_proteins=4, hsp_rows=5, queries_with_hits=2,
        queries_without_hits=2, failed_queries=1,
        queries_without_hits_and_without_logged_failure=1,
        failed_queries_with_query_hits=0, hsp_rows_above_1e_minus_5=0,
        query_ids_without_hits=["c", "d"],
        diagnostics=[dict(gene="c", query_failed=True, has_query_hits=False)])
    return [["a", "b", "c", "d"], ["a"], 3, [batch], content]


def test_exact_dispositions_keep_failed_and_unexplained_no_hits_separate():
    result = reconcile(*fixture())
    assert result["failed_queries"] == 1
    assert result["queries_without_hits_and_without_logged_failure"] == 1
    assert result["replay_hsp_rows"] == 2
    assert result["search_admitted"] is False


def test_actual_audit_and_batch_coverage_contract(tmp_path):
    from benchmark_tools.audit_orthomcl_search_table import audit_table
    from benchmark_tools.audit_orthomcl_blast import parse_diagnostics
    from benchmark_tools.admit_blast_recovery_batch import coverage, index_blocks

    fasta, table, replay, log = [tmp_path / name for name in ("all.fa", "all.blast", "replay.blast", "blast.log")]
    fasta.write_text("".join(f">{gene}\nACDEFGHIKL\n" for gene in "abcd"))
    def row(query, subject):
        return f"{query}\t{subject}\t100\t10\t0\t0\t1\t10\t1\t10\t1e-8\t50\n"
    # A failed query may still have incoming subject hits, which are not recovery.
    table.write_text(row("a", "c") + row("b", "b"))
    replay.write_text(row("b", "b"))
    log.write_text("[blastall] WARNING: [test] c: SetUpBlastSearch failed.\n")
    diagnostics = parse_diagnostics(log)
    blocks = index_blocks(replay, list("bcd"))
    batch = dict(query_ids=list("bcd"), query_blocks=blocks, diagnostics=diagnostics,
        coverage=coverage(list("bcd"), blocks, diagnostics))
    content = audit_table(table, fasta, log)
    assert content["failed_queries_with_subject_hits"] == 1
    result = reconcile(list("abcd"), ["a"], 1, [batch], content)
    assert result["failed_queries"] == 1 and result["hsp_rows"] == 2


@pytest.mark.parametrize("field", ["input_proteins", "hsp_rows", "queries_with_hits",
    "queries_without_hits", "failed_queries", "queries_without_hits_and_without_logged_failure",
    "failed_queries_with_query_hits", "hsp_rows_above_1e_minus_5"])
def test_wrong_total_rejected(field):
    args = fixture()
    args[-1][field] += 1
    with pytest.raises(ValueError, match="count disagrees"):
        reconcile(*args)


@pytest.mark.parametrize("change", ["missing_input", "duplicate_input", "overlap",
    "duplicate_batch", "swapped_absent", "changed_failure", "outgoing_failure",
    "bad_rows", "bad_retained_rows", "duplicate_diagnostic", "unknown_diagnostic",
    "bad_flag", "no_hit_partition", "batch_failure", "false_failure_flag"])
def test_identity_and_partition_errors(change):
    args = fixture()
    batch, content = args[3][0], args[4]
    if change == "missing_input":
        args[0].pop()
    elif change == "duplicate_input":
        args[0].append("a")
    elif change == "overlap":
        args[1].append("b")
    elif change == "duplicate_batch":
        args[3].append(deepcopy(batch))
    elif change == "swapped_absent":
        content["query_ids_without_hits"] = ["b", "d"]
    elif change == "changed_failure":
        content["diagnostics"][0]["gene"] = "d"
    elif change == "outgoing_failure":
        batch["coverage"]["failed_queries"] = ["b"]
        batch["diagnostics"] = {"b": dict(query_failed=True)}
    elif change == "bad_rows":
        batch["query_blocks"][0]["rows"] = True
    elif change == "bad_retained_rows":
        args[2] = 0
    elif change == "duplicate_diagnostic":
        content["diagnostics"] *= 2
    elif change == "unknown_diagnostic":
        content["diagnostics"][0]["gene"] = "unknown"
    elif change == "bad_flag":
        content["diagnostics"][0]["has_query_hits"] = True
    elif change == "no_hit_partition":
        batch["coverage"]["queries_without_hits"] = ["d"]
    elif change == "batch_failure":
        batch["diagnostics"] = {}
    elif change == "false_failure_flag":
        content["diagnostics"][0]["query_failed"] = "False"
    with pytest.raises(ValueError):
        reconcile(*args)
