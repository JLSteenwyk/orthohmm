import pytest
import json

from benchmark_tools import admit_blast_recovery_batch as admission
from benchmark_tools.prepare_ob_candidate_neighborhood import record

from benchmark_tools.admit_blast_recovery_batch import completed_task, index_blocks, coverage
from benchmark_tools.audit_orthomcl_search_table import audit_table


def line(query="a", subject="subject"):
    return f"{query}\t{subject}\t100.00\t3\t0\t0\t1\t3\t1\t3\t1e-9\t42\n".encode()


@pytest.mark.parametrize("problem", [None, "running", "exit", "cpu", "node", "duplicate", "missing"])
def test_scheduler_gate(problem):
    row = ["22103_0", "22104", "COMPLETED", "0:0", "bizon", "180", "00:01:00"]
    if problem in {"running", "exit", "cpu", "node"}:
        index, value = {"running": (2, "RUNNING"), "exit": (3, "1:0"), "cpu": (5, "1"), "node": (4, "other")}[problem]
        row[index] = value
    text = "JobID|JobIDRaw|State|ExitCode|NodeList|AllocCPUS|Elapsed\n"
    if problem != "missing":
        text += "|".join(row) + "\n"
    if problem == "duplicate":
        text += "|".join(row) + "\n"
    if problem:
        with pytest.raises(ValueError, match="completed"):
            completed_task(text, 0)
    else:
        assert completed_task(text, 0)["JobIDRaw"] == "22104"


def test_block_ranges_reproduce_every_byte(tmp_path):
    path = tmp_path / "hits"
    content = line()*2 + line("c")
    path.write_bytes(content)
    blocks = index_blocks(path, ["a", "nohit", "c"])
    assert [r["query"] for r in blocks] == ["a", "c"]
    assert [r["rows"] for r in blocks] == [2, 1]
    assert blocks[0]["start"] == 0
    assert blocks[0]["end"] == blocks[1]["start"] == len(line())*2
    assert blocks[1]["end"] == len(content)
    import hashlib
    for b in blocks:
        assert hashlib.sha256(content[b["start"]:b["end"]]).hexdigest() == b["sha256"]


@pytest.mark.parametrize("content,genes", [(line("other"), ["a"]), (line("b")+line("a"), ["a", "b"]),
    (line("a")+line("b")+line("a"), ["a", "b"]), (line()[:-1], ["a"]),
    (line().replace(b"42", b"4\x00"), ["a"]), (b"a\tb\n", ["a"]), (line(), ["a", "a"])])
def test_malformed_or_wrong_query_blocks_rejected(tmp_path, content, genes):
    path = tmp_path / "hits"
    path.write_bytes(content)
    with pytest.raises(ValueError):
        index_blocks(path, genes)


def test_empty_batch_output_and_failure_are_distinct(tmp_path):
    path = tmp_path / "hits"
    path.write_bytes(b"")
    blocks = index_blocks(path, ["failed", "nohit"])
    assert blocks == []
    result = coverage(["failed", "nohit"], blocks, {"failed": {"query_failed": True}})
    assert result["queries_without_hits"] == ["failed", "nohit"]
    assert result["failed_queries"] == ["failed"]
    assert result["no_hits_without_logged_failure"] == ["nohit"]


def test_unknown_diagnostic_rejected():
    with pytest.raises(ValueError, match="outside"):
        coverage(["a"], [], {"b": {"query_failed": True}})


def test_full_database_subject_is_valid_outside_query_subset(tmp_path):
    blast, fasta, log = (tmp_path / name for name in ("hits", "full.fa", "log"))
    blast.write_bytes(line())
    fasta.write_text(">a\nAAA\n>subject\nAAA\n>unqueried\nAAA\n")
    log.write_text("")
    blocks = index_blocks(blast, ["a"])
    structural = audit_table(blast, fasta, log)
    assert structural["input_proteins"] == 3
    assert structural["hsp_rows"] == blocks[0]["rows"] == 1
    assert coverage(["a"], blocks, {})["queries_without_hits"] == []


@pytest.mark.parametrize("replacement", [False, True])
@pytest.mark.parametrize("problem", [None, "empty", "changed_output", "wrong_contract",
    "wrong_identity", "partial_file", "failed_status", "above_cutoff", "unknown_subject"])
def test_admission_contract_end_to_end(tmp_path, monkeypatch, problem, replacement):
    index = 14 if replacement else 0
    array = "22160" if replacement else "22103"
    recovery = dict(original_scheduler={"State": "TIMEOUT"}, records=[],
                    replacement_task="22160_14", partial_rows_reused=False)
    monkeypatch.setattr(admission, "interrupted_attempt", lambda *a: recovery)
    directory = tmp_path / "batch_00"
    directory.mkdir()
    fasta, query = tmp_path / "full.fa", tmp_path / "query.fa"
    fasta.write_text(">a\nAAA\n>subject\nAAA\n")
    query.write_text(">a\nAAA\n")
    blast, log, timing = [directory / n for n in ("hits.blast", "blast.log", "blast.time.txt")]
    content = line()
    if problem == "empty":
        content = b""
    elif problem == "above_cutoff":
        content = content.replace(b"1e-9", b"1e-3")
    elif problem == "unknown_subject":
        content = line(subject="unknown")
    blast.write_bytes(content)
    log.write_text("")
    timing.write_text("Exit status: 0\n")
    expected = dict(directory=str(directory), command=["blastall", "-d", str(fasta)],
        checked_records=[record(fasta), record(query)],
        batch=dict(input=record(query), queries=1, first_query="a", last_query="a"))
    identity = dict(index=index, job_id="22104", array_job_id=array, node="bizon", started_epoch=1)
    initial = dict(**identity, preflight=expected, status="starting")
    terminal = dict(**identity, preflight=expected, status="native_batch_completed_pending_admission",
        finished_epoch=2, exit_code=0, outputs=[record(p) for p in (log, timing, blast)],
        search_admitted=False, reuse_authorized=False, publication_ready=False)
    if problem == "changed_output":
        blast.write_bytes(line()*2)
    elif problem == "wrong_contract":
        terminal["preflight"] = dict(expected, command=["different"])
    elif problem == "wrong_identity":
        initial["job_id"] = "999"
    elif problem == "partial_file":
        (directory / "hits.blast.partial").write_bytes(b"")
    elif problem == "failed_status":
        terminal["status"] = "failed_or_interrupted"
    (directory / "preflight.json").write_text(json.dumps(initial))
    (directory / "status.json").write_text(json.dumps(terminal))

    def check_output(command, **kwargs):
        if command[0] == "sacct":
            return ("JobID|JobIDRaw|State|ExitCode|NodeList|AllocCPUS|Elapsed\n"
                    f"{array}_{index}|22104|COMPLETED|0:0|bizon|180|00:01:00\n")
        if command[0] == "git":
            return admission.EXECUTOR_COMMIT + "\n"
        assert command[:3] == ["/home/bizon/anaconda3/bin/python", "-B", "-c"]
        return json.dumps(expected)

    monkeypatch.setattr(admission.subprocess, "check_output", check_output)
    monkeypatch.setattr(admission.subprocess, "run", lambda *a, **k: None)
    output = tmp_path / "admission.json"
    if problem not in (None, "empty"):
        with pytest.raises(ValueError):
            admission.admit(tmp_path, index, output, replacement)
        assert not output.exists()
        return
    result = admission.admit(tmp_path, index, output, replacement)
    assert result.get("replacement") == (recovery if replacement else None)
    assert json.loads(output.read_text()) == result
    assert result["batch_admitted"] is True
    assert result["search_admitted"] is result["reuse_authorized"] is False
    assert result["coverage"]["queries_with_hits"] == (0 if problem == "empty" else 1)
    if problem == "empty":
        assert result["coverage"]["no_hits_without_logged_failure"] == ["a"]
    with pytest.raises(FileExistsError):
        admission.admit(tmp_path, index, output, replacement)


@pytest.mark.parametrize("problem", [None, "wrong_array", "wrong_index", "running", "duplicate", "old_only"])
def test_replacement_scheduler_identity(problem):
    header = "JobID|JobIDRaw|State|ExitCode|NodeList|AllocCPUS|Elapsed\n"
    old = "22103_14|22143|TIMEOUT|0:0|bizon|180|1-10:51:21\n"
    new = "22160_14|22160|COMPLETED|0:0|bizon|180|00:01:00\n"
    if problem == "wrong_array":
        new = new.replace("22160_14", "99999_14")
    elif problem == "running":
        new = new.replace("COMPLETED", "RUNNING")
    elif problem == "duplicate":
        new *= 2
    elif problem == "old_only":
        new = ""
    index = 15 if problem == "wrong_index" else 14
    if problem:
        with pytest.raises(ValueError):
            completed_task(header + old + new, index, True)
    else:
        assert completed_task(header + old + new, 14, True)["JobIDRaw"] == "22160"
    with pytest.raises(ValueError):
        completed_task(header + old + new, 14)


@pytest.mark.parametrize("problem", [None, "bytes", "extra", "missing", "state", "raw", "duplicate"])
def test_interrupted_evidence_preserved(tmp_path, monkeypatch, problem):
    import hashlib
    directory = tmp_path / "benchmarks/results/qfo_blast_recovery_v1/batch_14_interrupted_22103_14_20260925"
    directory.mkdir(parents=True)
    for name in admission.INTERRUPTED_HASHES:
        (directory / name).write_bytes(b"fixture")
    hashes = {name: hashlib.sha256(b"fixture").hexdigest() for name in admission.INTERRUPTED_HASHES}
    monkeypatch.setattr(admission, "INTERRUPTED_HASHES", hashes)
    if problem == "bytes":
        (directory / "status.json").write_bytes(b"changed")
    elif problem == "extra":
        (directory / "hits.blast").touch()
    elif problem == "missing":
        (directory / "status.json").unlink()
    row = "22103_14|22143|TIMEOUT|0:0|bizon|180\n"
    if problem == "state":
        row = row.replace("TIMEOUT", "RUNNING")
    elif problem == "raw":
        row = row.replace("22143", "99999")
    elif problem == "duplicate":
        row *= 2
    accounting = "JobID|JobIDRaw|State|ExitCode|NodeList|AllocCPUS\n" + row
    if problem:
        with pytest.raises(ValueError):
            admission.interrupted_attempt(tmp_path, accounting)
    else:
        result = admission.interrupted_attempt(tmp_path, accounting)
        assert len(result["records"]) == 5
        assert result["partial_rows_reused"] is False
