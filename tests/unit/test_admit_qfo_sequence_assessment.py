import pytest

from benchmark_tools import admit_qfo_sequence_assessment as module


@pytest.mark.parametrize("problem", [None, "empty", "changed", "extra", "missing_trace", "extra_trace"])
def test_native_inventory(tmp_path, problem):
    results = tmp_path / "results"
    (results / "stats").mkdir(parents=True)
    trace = results / "stats/trace_1.txt"
    trace.write_text("fixture\n")
    metric = results / "metric.json"
    metric.write_text("{}\n")
    expected = [module.record(p) for p in sorted(results.rglob("*")) if p.is_file()]
    if problem == "empty":
        trace.unlink()
        metric.unlink()
        expected = []
    elif problem == "changed":
        metric.write_text("changed\n")
    elif problem == "extra":
        (results / "unrecorded").write_text("x")
    elif problem == "missing_trace":
        trace.unlink()
        expected = [module.record(metric)]
    elif problem == "extra_trace":
        (results / "stats/trace_2.txt").write_text("other\n")
        expected = [module.record(p) for p in sorted(results.rglob("*")) if p.is_file()]
    if problem:
        with pytest.raises(ValueError):
            module.outputs(results, expected)
    else:
        assert module.outputs(results, expected) == (expected, trace)


def test_live_job_blocks_admission(tmp_path, monkeypatch):
    def fail(*args):
        raise ValueError("still running")
    monkeypatch.setattr(module, "completed", fail)
    output = tmp_path / "admission.json"
    with pytest.raises(ValueError, match="still running"):
        module.admit(tmp_path, "all_hits", "job", "conversion", "sha", output)
    assert not output.exists()


def test_no_overwrite(tmp_path):
    output = tmp_path / "admission.json"
    output.write_text("existing")
    with pytest.raises(FileExistsError):
        module.admit(tmp_path, "all_hits", "job", "conversion", "sha", output)
    assert output.read_text() == "existing"
