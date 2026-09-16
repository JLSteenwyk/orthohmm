from pathlib import Path

import pytest

from benchmark_tools import assemble_orthobench_factorial as assembler


def accounting(state="FAILED"):
    return "JobID|State|ExitCode\n" + "".join(f"21248_{i}|{state}|1:0\n" for i in range(4))


def test_accepts_complete_terminal_task_set_without_relabeling():
    rows = assembler.require_terminal_array(accounting())
    assert len(rows) == 4
    assert all(r["State"] == "FAILED" for r in rows)


@pytest.mark.parametrize("state", ["RUNNING", "PENDING", "COMPLETING", "UNKNOWN"])
def test_rejects_nonterminal_tasks(state):
    with pytest.raises(ValueError):
        assembler.require_terminal_array(accounting(state))


@pytest.mark.parametrize("change", ["duplicate", "missing"])
def test_rejects_ambiguous_task_set(change):
    text = accounting()
    text = text + "21248_0|FAILED|1:0\n" if change == "duplicate" else text.replace("21248_3|FAILED|1:0\n", "")
    with pytest.raises(ValueError):
        assembler.require_terminal_array(text)


def test_running_array_never_reaches_validation_or_scoring(tmp_path, monkeypatch):
    monkeypatch.setattr(assembler.subprocess, "check_output", lambda *a, **k: accounting("RUNNING"))
    def forbidden(*args, **kwargs):
        pytest.fail("Read inference or reference outcomes before terminal gate")
    monkeypatch.setattr(assembler, "validate", forbidden)
    monkeypatch.setattr(assembler, "read_frozen", forbidden)
    monkeypatch.setattr(assembler, "score_partition", forbidden)
    output = tmp_path / "results"
    with pytest.raises(ValueError, match="no partial scoring"):
        assembler.assemble(Path.cwd(), output)
    assert not output.exists()


def test_official_rounding_is_percentage_not_fraction():
    exact = {"f_score": 74.106, "precision": 81.770, "recall": 67.755, "exact_refogs": 20}
    printed = {"f_score": 74.1, "precision": 81.8, "recall": 67.8, "exact_refogs": 20}
    assembler.compare_official(exact, printed)


@pytest.mark.parametrize("change", ["f_score", "precision", "recall", "exact_refogs", "missing", "nan"])
def test_rejects_scorer_disagreement(change):
    exact = {"f_score": 70.0, "precision": 70.0, "recall": 70.0, "exact_refogs": 20}
    printed = exact.copy()
    if change == "missing":
        del printed["precision"]
    elif change == "nan":
        printed["recall"] = float("nan")
    else:
        printed[change] += 1
    with pytest.raises(ValueError):
        assembler.compare_official(exact, printed)
