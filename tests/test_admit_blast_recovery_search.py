from copy import deepcopy
from pathlib import Path

import pytest

import benchmark_tools.admit_blast_recovery_search as module


def accounting():
    return "JobID|State|ExitCode|NodeList|AllocCPUS|ReqMem|Elapsed\n22150|COMPLETED|0:0|bizon|2|64G|00:10:00\n"


def test_completed_merge_allocation():
    assert module.completed_merge(accounting())["JobID"] == "22150"


@pytest.mark.parametrize("old,new", [("COMPLETED", "RUNNING"), ("0:0", "1:0"),
    ("bizon", "other"), ("|2|", "|1|"), ("64G", "8G"), ("22150", "21713")])
def test_wrong_merge_execution(old, new):
    with pytest.raises(ValueError):
        module.completed_merge(accounting().replace(old, new))


def test_duplicate_scheduler_record():
    text = accounting()
    with pytest.raises(ValueError):
        module.completed_merge(text + text.splitlines()[1] + "\n")


def status(directory):
    source = dict(path="frozen/driver.py", bytes=1, sha256="fixture")
    candidate = dict(status="merged_candidate_requires_full_admission",
        path=str(directory / "table/all.blast.candidate"), bytes=10, sha256="fixture",
        search_admitted=False, reuse_authorized=False, publication_ready=False)
    return dict(status="merged_candidate_pending_full_table_admission", job_id="22150",
        source=source, candidate=candidate, selected_log={"path": str(directory / "selected.blast.log")},
        search_admitted=False, reuse_authorized=False, publication_ready=False)


def test_merge_output_binding():
    directory = Path("/fixture")
    report = status(directory)
    assert module.validate_merge(report, report["source"], directory)["bytes"] == 10


@pytest.mark.parametrize("change", ["job", "source", "status", "admitted", "candidate_path",
    "candidate_status", "candidate_admitted", "log"])
def test_wrong_merge_evidence(change):
    directory = Path("/fixture")
    report = status(directory)
    source = deepcopy(report["source"])
    if change == "job":
        report["job_id"] = "21713"
    elif change == "source":
        report["source"]["sha256"] = "changed"
    elif change == "status":
        report["status"] = "merge_failed_preserved"
    elif change == "admitted":
        report["search_admitted"] = True
    elif change == "candidate_path":
        report["candidate"]["path"] = "/other/all.blast"
    elif change == "candidate_status":
        report["candidate"]["status"] = "partial"
    elif change == "candidate_admitted":
        report["candidate"]["reuse_authorized"] = True
    elif change == "log":
        report["selected_log"]["path"] = "/other/blast.log"
    with pytest.raises(ValueError):
        module.validate_merge(report, source, directory)


def test_pending_merge_rejected_before_artifacts(tmp_path, monkeypatch):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
        accounting().replace("COMPLETED", "PENDING"))
    with pytest.raises(ValueError, match="completed frozen"):
        module.admit(tmp_path, tmp_path / "admission")
    assert not (tmp_path / "admission").exists()


def test_existing_output_preserved(tmp_path):
    output = tmp_path / "admission"
    output.mkdir()
    with pytest.raises(FileExistsError):
        module.admit(tmp_path, output)
