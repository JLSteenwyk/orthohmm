import json

import pytest

from benchmark_tools import readback_large_copy_components as reader
from benchmark_tools import run_large_copy_components as runner
from tests.unit.test_run_large_copy_components import worker_fixture


def completed(tmp_path, monkeypatch, dataset="OrthoBench"):
    path, calls = worker_fixture(tmp_path, monkeypatch, dataset)
    construct = runner.construct
    monkeypatch.setattr(runner, "construct", lambda *a: (construct(*a)[0], {
        "component_counts": {}, "ownership": {"genes": runner.EXPECTED_GENES[dataset]}}))
    monkeypatch.setenv("SLURM_JOB_ID", "999")
    producer = runner.run(path)
    assert producer["status"] == "independent_content_readback_verified"
    monkeypatch.setattr(reader, "replay_native", runner.replay_native)
    monkeypatch.setattr(reader, "readback_scores", runner.readback_scores)
    monkeypatch.setattr(reader.subprocess, "check_output", lambda *a, **k:
        "JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n999|COMPLETED|0:0|00:01:00|bizon|16\n")
    return path, calls


@pytest.mark.parametrize("dataset", list(runner.EXPECTED_GENES))
def test_independent_reader_reconstructs_candidates_and_rereads_reused_scores(tmp_path, monkeypatch, dataset):
    path, calls = completed(tmp_path, monkeypatch, dataset)
    calls.clear()
    result = reader.readback(path, tmp_path / "independent")
    assert result["status"] == "terminal_component_candidates_and_raw_scores_independently_verified"
    assert result["dataset"] == dataset and calls == ["native_readback", "score_readback"]
    assert "native_execution" not in result["resources"]
    with pytest.raises(ValueError, match="untouched"):
        reader.readback(path, tmp_path / "independent")


@pytest.mark.parametrize("mutation", ["status", "postflight", "coverage", "scores", "products", "trace"])
def test_reader_rejects_altered_claims_or_pinned_outputs(tmp_path, monkeypatch, mutation):
    path, calls = completed(tmp_path, monkeypatch)
    producer_path = tmp_path / "results.json"
    producer = json.loads(producer_path.read_text())
    if mutation == "status":
        producer["status"] = "native_verified_scoring_pending"
    elif mutation == "postflight":
        producer["inputs_reverified_unchanged"] = False
    elif mutation == "coverage":
        producer["coverage"]["genes"] -= 1
    elif mutation == "scores":
        producer["scores"][next(iter(producer["scores"]))] = .6
    elif mutation == "products":
        producer["native_products"] = {"unexpected": {}}
    else:
        producer["structural_controls"]["components500"]["copy_split_trace"]["sha256"] = "bad"
    producer_path.write_text(json.dumps(producer))
    with pytest.raises(ValueError):
        reader.readback(path, tmp_path / "independent")
    assert not (tmp_path / "independent/report.json").exists()


def test_reader_rejects_a_live_scheduler_before_reconstruction(tmp_path, monkeypatch):
    path, calls = completed(tmp_path, monkeypatch)
    monkeypatch.setattr(reader.subprocess, "check_output", lambda *a, **k:
        "JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n999|RUNNING|0:0|00:01:00|bizon|16\n")
    monkeypatch.setattr(runner, "construct", lambda *a: pytest.fail("No reconstruction before terminal gate"))
    with pytest.raises(ValueError, match="terminal"):
        reader.readback(path, tmp_path / "independent")


def test_reader_rejects_changed_component_reconstruction(tmp_path, monkeypatch):
    path, calls = completed(tmp_path, monkeypatch)
    construct = runner.construct
    monkeypatch.setattr(runner, "construct", lambda *a: (construct(*a)[0], {
        "component_counts": {"families": 1}, "ownership": {"genes": runner.EXPECTED_GENES["OrthoBench"]}}))
    with pytest.raises(ValueError, match="activation or ownership"):
        reader.readback(path, tmp_path / "independent")
