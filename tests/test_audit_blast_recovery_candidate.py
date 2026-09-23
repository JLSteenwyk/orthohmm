import hashlib
import json

import pytest

from benchmark_tools.audit_blast_recovery_candidate import audit


def inputs(tmp_path, change=None):
    fasta, blast, log, prefix, batch, output = [tmp_path / name for name in
        ("all.fa", "all.blast", "blast.log", "prefix.jsonl", "batch.json", "audit.json")]
    fasta.write_text("".join(f">{gene}\nACDEFGHIKL\n" for gene in "abcd"))
    a = "a\tc\t100\t10\t0\t0\t1\t10\t1\t10\t1e-8\t50\n"
    b = a.replace("a\tc", "b\tb")
    if change == "cutoff":
        b = b.replace("1e-8", "1e-2")
    if change == "coordinates":
        b = b.replace("\t1\t10\t1\t10\t", "\t1\t11\t1\t10\t")
    def block(q, row, start=0, final=False):
        return dict(query=q, start=start, end=start+len(row), rows=1,
            sha256=hashlib.sha256(row.encode()).hexdigest(), final_observed_query=final)
    original = [block("a", a), block("b", b, len(a), True)]
    if change == "boundary":
        original[-1]["final_observed_query"] = False
    if change == "post_boundary":
        original.append(block("x", a, len(a)+len(b)))
    prefix.write_text("".join(json.dumps(row)+"\n" for row in original))
    report = dict(status="recovery_batch_execution_and_rows_verified", batch_admitted=True,
        query_ids=list("bcd"), query_blocks=[block("b", b)],
        coverage=dict(queries_without_hits=["c", "d"], failed_queries=["c"]),
        diagnostics={"c": dict(query_failed=True)})
    if change == "batch_status":
        report["batch_admitted"] = False
    batch.write_text(json.dumps(report))
    log.write_text("[blastall] WARNING: [test] c: SetUpBlastSearch failed.\n")
    if change == "failure":
        log.write_text("")
    candidate = b+a if change == "order" else a+b
    if change == "bytes":
        candidate = candidate.replace("\t50\n", "\t51\n")
    if change == "truncate":
        candidate = candidate.rstrip("\n")
    blast.write_text(candidate)
    return blast, fasta, log, prefix, [batch], output


def test_integrated_content_audit_never_admits_execution(tmp_path):
    args = inputs(tmp_path)
    result = audit(*args)
    assert result["query_blocks"] == 2
    assert result["content"]["failed_queries_with_subject_hits"] == 1
    assert result["dispositions"]["queries_without_hits_and_without_logged_failure"] == 1
    assert result["search_admitted"] is result["downstream_execution_authorized"] is False
    assert json.loads(args[-1].read_text()) == result
    with pytest.raises(FileExistsError):
        audit(*args)


@pytest.mark.parametrize("change", ["cutoff", "coordinates", "boundary", "post_boundary",
    "batch_status", "failure", "order", "bytes", "truncate"])
def test_corrupt_candidate_or_inventory_has_no_success_report(tmp_path, change):
    args = inputs(tmp_path, change)
    with pytest.raises(ValueError):
        audit(*args)
    assert not args[-1].exists()
