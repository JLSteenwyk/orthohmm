import hashlib
import json

import pytest

import benchmark_tools.run_blast_recovery_merge as module


def accounting():
    rows = ["JobID|State|ExitCode|NodeList|AllocCPUS|Elapsed", "22148|COMPLETED|0:0|bizon|2|00:01:25"]
    rows += [f"{job}_{i}|COMPLETED|0:0|bizon|{cpu}|00:01:00"
             for job, cpu in ((22103, 180), (22105, 2)) for i in range(20)]
    return "\n".join(rows)+"\n"


def test_complete_prerequisites():
    module.require_completed(accounting())


@pytest.mark.parametrize("change", ["missing", "failed", "pending", "cpu", "duplicate", "node"])
def test_bad_prerequisites(change):
    text = accounting()
    target = "22103_19|COMPLETED|0:0|bizon|180|00:01:00"
    replacement = {"missing": "", "failed": target.replace("0:0", "1:0"),
        "pending": target.replace("COMPLETED", "PENDING"), "cpu": target.replace("180", "2"),
        "duplicate": target+"\n"+target, "node": target.replace("bizon", "other")}[change]
    with pytest.raises(ValueError):
        module.require_completed(text.replace(target, replacement))


def diagnostic(failed=True, line=1):
    return dict(query_failed=failed, messages=[dict(level="WARNING", category="setup_failure" if failed else "selenocysteine_replacement",
        message="SetUpBlastSearch failed." if failed else "replacement", line=line)])


def test_diagnostics_ignore_line_numbers_not_multiplicity():
    original, replay = {"a": diagnostic(), "b": diagnostic(False)}, {"a": diagnostic(line=7)}
    assert set(module.select_diagnostics(original, replay, ["a"])) == {"a", "b"}
    replay["a"]["messages"] *= 2
    with pytest.raises(ValueError):
        module.select_diagnostics(original, replay, ["a"])


@pytest.mark.parametrize("original,replay,ids", [({"a": diagnostic()}, {}, ["a"]),
    ({"a": diagnostic()}, {}, []), ({}, {"a": diagnostic()}, [])])
def test_missing_or_misplaced_failure(original, replay, ids):
    with pytest.raises(ValueError):
        module.select_diagnostics(original, replay, ids)


def test_actual_incomplete_panel_fails_before_copy(tmp_path, monkeypatch):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: accounting().replace("22105_19|COMPLETED", "22105_19|PENDING"))
    with pytest.raises(ValueError, match="22105_19"):
        module.prepare(tmp_path, tmp_path/"new")
    assert not (tmp_path/"new").exists()


def test_changed_review_rejected_before_copy(tmp_path, monkeypatch):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: accounting())
    path = tmp_path/"benchmark_tools/results/QFO_BLAST_PREFIX_REUSE_REVIEW_20260923.md"
    path.parent.mkdir(parents=True)
    path.write_text("changed review\n")
    with pytest.raises(ValueError, match="Changed merge gate"):
        module.prepare(tmp_path, tmp_path/"new")
    assert not (tmp_path/"new").exists()


def test_transaction_records_prerequisite_failure(tmp_path, monkeypatch):
    monkeypatch.setenv("SLURM_JOB_ID", "test")
    monkeypatch.setenv("SLURM_JOB_NODELIST", "bizon")
    def fail(*args):
        raise ValueError("Prerequisite not ready")
    monkeypatch.setattr(module, "prepare", fail)
    with pytest.raises(ValueError):
        module.run(tmp_path, tmp_path/"output")
    status = json.loads((tmp_path/"output/status.json").read_text())
    assert status["status"] == "merge_failed_preserved" and status["search_admitted"] is False
    assert not (tmp_path/"output/table").exists()


def test_candidate_transaction_never_admits_search(tmp_path, monkeypatch):
    monkeypatch.setenv("SLURM_JOB_ID", "test")
    monkeypatch.setenv("SLURM_JOB_NODELIST", "bizon")
    def row(q):
        return (q+"\ts\t100\t10\t0\t0\t1\t10\t1\t10\t1e-8\t50\n").encode()
    prefix = tmp_path/"prefix"
    prefix.write_bytes(row("a"))
    replay = tmp_path/"replay"
    replay.write_bytes(row("b"))
    def block(q, source):
        return dict(query=q, path=source, start=0, end=len(row(q)), rows=1,
                    sha256=hashlib.sha256(row(q)).hexdigest(), final_observed_query=False)
    blocks = tmp_path/"blocks.jsonl"
    blocks.write_text(json.dumps(block("a", "prefix"))+"\n")
    review = tmp_path/"benchmark_tools/results/QFO_BLAST_PREFIX_REUSE_REVIEW_20260923.md"
    review.parent.mkdir(parents=True)
    review.write_text("fixture only\n")
    prepared = dict(genes=["a", "b"], replay_ids=["b"], blocks=[block("b", "new")],
        sources={"prefix": prefix, "new": replay}, recheck={"result": dict(retained_end=len(row("a")), retained_rows=1)},
        prefix_audit={"blocks": {"path": str(blocks)}}, panel={"coverage": {"hsp_rows": 1}},
        raw_diagnostics={}, selected_diagnostics={}, records=[], accounting="fixture")
    monkeypatch.setattr(module, "prepare", lambda *args: prepared)
    status = module.run(tmp_path, tmp_path/"output")
    assert status["status"] == "merged_candidate_pending_full_table_admission"
    assert status["search_admitted"] is status["reuse_authorized"] is False
    assert (tmp_path/"output/table/all.blast.candidate").read_bytes() == row("a")+row("b")
    with pytest.raises(FileExistsError):
        module.run(tmp_path, tmp_path/"output")
