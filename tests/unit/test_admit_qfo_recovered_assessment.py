import pytest

from benchmark_tools import admit_qfo_recovered_assessment as audit


def trace():
    names = ["validate_input_file", "convertPredictions", "consolidate", "vgnc_benchmark (1)",
             "ec_benchmark (1)", "go_benchmark (1)", "fas_benchmark (1)",
             "reference_genetrees_benchmark (SwissTrees)", "reference_genetrees_benchmark (TreeFam-A)"]
    names.extend(f"scheduleMetrics ({i})" for i in range(1, 7))
    return [[str(i), n, "COMPLETED", "-" if n.startswith("scheduleMetrics (") else "0"] for i, n in enumerate(names)]


@pytest.mark.parametrize("problem", [None, "missing", "duplicate", "failed", "cached", "exit", "task_id"])
def test_native_trace_contract(problem):
    rows = trace()
    if problem == "missing":
        rows.pop()
    elif problem == "duplicate":
        rows[-1] = rows[0]
    elif problem == "failed":
        rows[0][2] = "FAILED"
    elif problem == "cached":
        rows[0][2] = "CACHED"
    elif problem == "exit":
        rows[0][3] = "1"
    elif problem == "task_id":
        rows[0][0] = rows[1][0]
    text = "task_id\tname\tstatus\texit\n" + "\n".join("\t".join(r) for r in rows)
    if problem:
        with pytest.raises(ValueError):
            audit.validate_trace(text)
    else:
        assert len(audit.validate_trace(text)) == 15


def test_nonterminal_array_rejected_before_output_creation(tmp_path, monkeypatch):
    monkeypatch.setattr(audit.subprocess, "check_output", lambda *a, **k:
                        "JobID|JobIDRaw|State|ExitCode|Elapsed\n21548_0|21548|RUNNING|0:0|00:01:00\n")
    with pytest.raises(ValueError, match="terminal"):
        audit.admit(tmp_path, tmp_path / "output")
    assert not (tmp_path / "output").exists()


def test_existing_output_refused(tmp_path):
    with pytest.raises(FileExistsError):
        audit.admit(tmp_path, tmp_path)
