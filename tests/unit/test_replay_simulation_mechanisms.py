"""Actual checked-count integration and portable replay guards."""

import json
from pathlib import Path

import pytest

from benchmark_tools import replay_simulation_mechanisms as replay


ROOT = Path(__file__).resolve().parents[2]


def reports():
    return {role: json.loads((ROOT / "benchmark_tools/results" / name).read_text())
            for role, (name, _) in replay.INPUTS.items()}


def test_actual_three_diagnostics_have_consistent_counts_and_complete_scope():
    summary = replay.summarize(**reports())
    assert summary["screened_cells"] == 70
    assert summary["screened_candidates"] == 10125
    assert summary["eligible_candidates"] == 1377
    assert summary["upstream_true_pair_rows"] == 163527
    assert summary["residual"]["counts"] == {"tp": 380, "fp": 66, "fn": 62}
    divergent = next(r for r in summary["conditions"] if r["condition"] == "divergent")
    assert divergent["oracle_fn"] == 12766
    assert divergent["cross_candidate_fn"] == 12736
    assert divergent["upstream_totals"]["different_graph_components"] == 11442
    assert divergent["upstream_totals"]["connected_but_separated"] == 1294


def test_sparse_counter_omitted_zero_native_errors_are_not_missing_measurements():
    data = reports()
    zeros = [c for c in data["upstream"]["cells"] if "native_fn" not in c["totals"]]
    assert zeros
    original = {c["label"]: c for c in data["oracle"]["cells"]}
    assert all(original[c["label"]]["arms"]["inferred"]["fn"] == 0 for c in zeros)
    assert replay.summarize(**data)["screened_cells"] == 70


@pytest.mark.parametrize("change", ["missing_cell", "duplicate_cell", "metric", "decomposition",
    "upstream_partition", "condition_mean", "residual_missing", "residual_duplicate", "residual_class", "not_admitted"])
def test_arithmetic_replay_rejects_incoherent_cross_diagnostic_evidence(change):
    data = reports()
    if change == "missing_cell":
        data["oracle"]["cells"].pop()
    elif change == "duplicate_cell":
        data["upstream"]["cells"][-1] = data["upstream"]["cells"][0]
    elif change == "metric":
        data["oracle"]["cells"][0]["arms"]["inferred"]["f1"] += 0.1
    elif change == "decomposition":
        data["oracle"]["cells"][0]["residual_by_arm"]["inferred"]["oracle_eligible"]["fn"] += 1
    elif change == "upstream_partition":
        data["upstream"]["cells"][0]["totals"]["connected_but_separated"] = 1
    elif change == "condition_mean":
        data["oracle"]["summary"][0]["mean_metrics"]["inferred"]["f1"] += 0.1
    elif change == "residual_missing":
        data["residual"]["candidates"].pop()
    elif change == "residual_duplicate":
        data["residual"]["candidates"].append(data["residual"]["candidates"][0])
    elif change == "residual_class":
        data["residual"]["candidates"][0]["error_classes"]["unsupported_satellite_constraint"] += 1
    else:
        data["residual"]["status"] = "failed"
    with pytest.raises(ValueError):
        replay.summarize(**data)


def test_render_preserves_all_conditions_and_noncausal_limits():
    text = replay.render(replay.summarize(**reports()))
    assert all(f"| {condition} |" in text for condition in replay.CONDITIONS)
    assert "unsupported_satellite_constraint | 62" in text
    assert "163,527" in text
    assert "not raw admission, native inference" in text
    assert "do not isolate search-stage causes" in text


def test_copied_checked_reports_run_without_historical_metadata_reads(tmp_path):
    directory = tmp_path / "inputs"
    directory.mkdir()
    for name, _ in replay.INPUTS.values():
        (directory / name).write_bytes((ROOT / "benchmark_tools/results" / name).read_bytes())
    result = replay.run(directory, tmp_path / "result")
    assert result["historical_native_paths_accessed"] is False
    assert (tmp_path / "result/summary.md").read_text() == replay.render(result["summary"])
    assert result["publication_ready"] is False
    with pytest.raises(ValueError, match="existing output"):
        replay.run(directory, tmp_path / "result")


def test_changed_input_rejected_before_creating_output(tmp_path):
    directory = tmp_path / "inputs"
    directory.mkdir()
    for name, _ in replay.INPUTS.values():
        data = (ROOT / "benchmark_tools/results" / name).read_bytes()
        (directory / name).write_bytes(data + b" ")
    with pytest.raises(ValueError, match="Changed checked input"):
        replay.run(directory, tmp_path / "output")
    assert not (tmp_path / "output").exists()
