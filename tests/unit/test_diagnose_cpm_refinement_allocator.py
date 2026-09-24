import copy

import pytest

from benchmark_tools.diagnose_cpm_refinement_allocator import environment, record, run, validate_child


def test_diagnostic_environment(monkeypatch, tmp_path):
    for key in ("PYTHONHOME", "LD_PRELOAD", "LD_LIBRARY_PATH"):
        monkeypatch.setenv(key, "bad")
    monkeypatch.setenv("OMP_NUM_THREADS", "32")
    env = environment(tmp_path)
    assert not set(("PYTHONHOME", "LD_PRELOAD", "LD_LIBRARY_PATH")) & env.keys()
    assert env["PYTHONMALLOC"] == "debug"
    assert env["PYTHONFAULTHANDLER"] == "1"
    assert env["OMP_NUM_THREADS"] == env["OPENBLAS_NUM_THREADS"] == env["MKL_NUM_THREADS"] == "1"
    assert env["PYTHONPATH"] == str(tmp_path)


@pytest.mark.parametrize("problem", [None, "metadata", "bytes", "record", "duplicate"])
def test_child_check(tmp_path, problem):
    path = tmp_path / "output.txt"
    path.write_text("a a\n" if problem == "duplicate" else "a b\n")
    reference = dict(groups=1, genes=2, modules=["frozen"], accuracy_evaluated=False, output=record(path))
    reference["output"]["path"] = "original"
    child = copy.deepcopy(reference)
    child["output"] = record(path)
    if problem == "metadata":
        child["modules"] = ["changed"]
    elif problem == "bytes":
        reference["output"]["sha256"] = "wrong"
    elif problem == "record":
        child["output"]["path"] = "wrong"
    if problem:
        with pytest.raises(ValueError):
            validate_child(child, reference, path, ["a", "b"])
    else:
        assert validate_child(child, reference, path, ["a", "b"]) == dict(genes=2, groups=1)


def test_unscheduled_run_cannot_create_output(tmp_path, monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    with pytest.raises(ValueError, match="allocation"):
        run(tmp_path, "fixture")
    assert not list(tmp_path.iterdir())
