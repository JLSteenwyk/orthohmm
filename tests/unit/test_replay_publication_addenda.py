"""Portable reporting replay keeps identities, unavailable costs and bounded scope."""

from copy import deepcopy
import hashlib
import json
from pathlib import Path

import pytest

from benchmark_tools import replay_publication_addenda as replay


ROOT = Path(__file__).resolve().parents[2]
PATHS = {
    "inventory": "development_family_inventory_20261004/inventory.json",
    "trace": "candidate_trace_variation_20261004/diagnostic.json",
    "native": "factorial_native_resource_linkage_20261004/linkage.json",
    "qfo": "qfo_orthohmm_stage_metadata_20261004/register.json",
    "factorial": "factorial_retained_resources_20261004/resources.json",
    "all_tools": "all_benchmark_metadata_integrated_20261004_v3/register.json",
}


def reports():
    return {key: json.loads((ROOT / "benchmark_tools/results" / path).read_bytes()) for key, path in PATHS.items()}


def test_actual_reporting_projection_matches_all_retained_claims():
    docs = reports()
    before = deepcopy(docs)
    result = replay.summarize(**docs)
    assert docs == before
    assert result["exposure"]["files"] == 1899
    assert result["exposure"]["families"] == 88
    assert result["exposure"]["associations"] == 9156
    assert result["exposure"]["summary"] == docs["inventory"]["summary"]
    assert [r["unattached_ids"] for r in result["candidate_trace"]["fixture_cases"]] == [[8], [0], [7]]
    assert all(r["unattached_count"] == 1 for r in result["candidate_trace"]["fixture_cases"])
    assert result["candidate_trace"]["retained_trace_common_merges"] == [8428, 8439, 8426]
    assert [r["resources"]["wall_seconds"]["median"] for r in result["native_costs"]["summaries"]] == [2720.077391309, 4392.871580988]
    assert result["register"]["metric_positions"] == 64
    assert result["register"]["secondary_mean_positions"] == 8
    assert result["register"]["stage_associations"] == 5
    assert result["register"]["distinct_stage_observations"] == 4
    assert result["register"]["original_full_cost_cells_unavailable"] == 16


@pytest.mark.parametrize("damage", ["file_duplicate", "byte_count", "block_duplicate", "unknown_family", "family_projection", "causal", "fixture_id", "fixture_partition", "native_repeat", "native_median", "native_invalid", "score_changed", "full_cost", "factorial_cost", "secondary_mean"])
def test_corruption_cannot_become_an_accepted_summary(damage):
    d = reports()
    if damage == "file_duplicate": d["inventory"]["files"].append(d["inventory"]["files"][0])
    elif damage == "byte_count": d["inventory"]["snapshot_json_bytes"] += 1
    elif damage == "block_duplicate": d["inventory"]["evaluations"].append(d["inventory"]["evaluations"][0])
    elif damage == "unknown_family": d["inventory"]["evaluations"][0]["families"].append("unknown")
    elif damage == "family_projection": d["inventory"]["family_rows"][0]["scored_evidence_blocks"] += 1
    elif damage == "causal": d["inventory"]["family_rows"][0]["causal_tuning_influence_established"] = True
    elif damage == "fixture_id": d["trace"]["fixture"]["cases"][0]["unattached_satellites"] = [7]
    elif damage == "fixture_partition": d["trace"]["fixture"]["cases"][0]["partition"][0].append(8)
    elif damage == "native_repeat": d["native"]["points"][0]["repeat"] = 1
    elif damage == "native_median": d["native"]["summaries"][0]["resources"]["wall_seconds"]["median"] += 1
    elif damage == "native_invalid": d["native"]["points"][0]["resources"]["wall_seconds"] = float("nan")
    elif damage == "score_changed": d["qfo"]["rows"][0]["scores"]["weighted_refog_F1"] = 0
    elif damage == "factorial_cost": d["factorial"]["rows"][0]["full_pipeline_wall_s"] = 1
    elif damage == "secondary_mean":
        for report in (d["all_tools"], d["qfo"]):
            next(r for r in report["rows"] if r["dataset"] == "QfO")["secondary_mean"] += .1
    else:
        next(r for r in d["qfo"]["rows"] if "qfo_stage_provenance" in r)["qfo_stage_provenance"]["full_pipeline_wall_s"] = 1
    with pytest.raises(ValueError): replay.summarize(**d)


def test_actual_input_hashes_match_the_delivered_metadata():
    for key, path in PATHS.items():
        assert hashlib.sha256((ROOT / "benchmark_tools/results" / path).read_bytes()).hexdigest() == replay.INPUTS[key][1]


def test_copied_inputs_ignore_metadata_paths_and_refuse_changed_bytes(tmp_path, monkeypatch):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    pins = {}
    for key, document in reports().items():
        document["nonexistent_historical_path"] = "/not-an-input/old-scientific-run"
        name = replay.INPUTS[key][0]
        data = json.dumps(document).encode()
        (inputs / name).write_bytes(data)
        pins[key] = name, hashlib.sha256(data).hexdigest()
    monkeypatch.setattr(replay, "INPUTS", pins)
    result = replay.run(inputs, tmp_path / "out")
    assert result["raw_or_historical_paths_accessed"] is False
    expected = ROOT / "benchmark_tools/results/development_family_inventory_20261004/family_evidence.tsv"
    assert (tmp_path / "out/family_evidence.tsv").read_bytes() == expected.read_bytes()
    with pytest.raises(FileExistsError): replay.run(inputs, tmp_path / "out")
    (inputs / pins["inventory"][0]).write_text("changed")
    with pytest.raises(ValueError): replay.run(inputs, tmp_path / "changed")
    assert not (tmp_path / "changed").exists()


def test_broken_symlink_output_is_occupied(tmp_path):
    target = tmp_path / "output"
    target.symlink_to(tmp_path / "absent")
    with pytest.raises(FileExistsError): replay.run(tmp_path, target)
