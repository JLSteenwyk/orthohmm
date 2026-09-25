from pathlib import Path

import pytest

from benchmark_tools.admit_blast_recovery_search import completed_merge, validate_merge


@pytest.mark.parametrize("replacement", [False, True])
@pytest.mark.parametrize("problem", [None, "job", "state", "exit", "cpu", "memory", "duplicate"])
def test_merge_scheduler(replacement, problem):
    job = "22162" if replacement else "22150"
    row = [job, "COMPLETED", "0:0", "bizon", "2", "64G"]
    if problem in {"job", "state", "exit", "cpu", "memory"}:
        index, value = {"job": (0, "99999"), "state": (1, "RUNNING"),
                        "exit": (2, "1:0"), "cpu": (4, "180"), "memory": (5, "1G")}[problem]
        row[index] = value
    text = "JobID|State|ExitCode|NodeList|AllocCPUS|ReqMem\n" + ("|".join(row) + "\n") * (2 if problem == "duplicate" else 1)
    if problem:
        with pytest.raises(ValueError):
            completed_merge(text, replacement)
    else:
        assert completed_merge(text, replacement)["JobID"] == job
        with pytest.raises(ValueError):
            completed_merge(text, not replacement)


@pytest.mark.parametrize("replacement", [False, True])
@pytest.mark.parametrize("problem", [None, "job", "source", "path", "promoted", "candidate"])
def test_merge_content_identity(replacement, problem):
    directory = Path("/fresh")
    source = dict(path="/frozen/run_blast_recovery_merge.py", sha256="source", bytes=1)
    flags = dict(search_admitted=False, reuse_authorized=False, publication_ready=False)
    status = dict(status="merged_candidate_pending_full_table_admission",
                  job_id="22162" if replacement else "22150", source=source,
                  selected_log=dict(path="/fresh/selected.blast.log"), **flags)
    status["candidate"] = dict(status="merged_candidate_requires_full_admission",
                               path="/fresh/table/all.blast.candidate", bytes=10, sha256="table", **flags)
    if problem == "job":
        status["job_id"] = "99999"
    elif problem == "source":
        status["source"] = dict(source, sha256="changed")
    elif problem == "path":
        status["candidate"]["path"] = "/old/table"
    elif problem == "promoted":
        status["search_admitted"] = True
    elif problem == "candidate":
        status["candidate"]["search_admitted"] = True
    if problem:
        with pytest.raises(ValueError):
            validate_merge(status, source, directory, replacement)
    else:
        assert validate_merge(status, source, directory, replacement)["bytes"] == 10
