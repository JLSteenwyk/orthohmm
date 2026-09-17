from pathlib import Path

import pytest

from benchmark_tools.run_simulation_tree_mode_control import (
    compare_pairs, compare_inventory, native_tree, artifact_inventory, run, baseline_records,
)


def test_pairs_compare_membership_not_orientation_or_duplicates():
    result = compare_pairs([("a", "b"), ("b", "a")], [("b", "a")])
    assert result == {"identical": True, "before_pairs": 1, "after_pairs": 1, "added": [], "removed": []}


def test_pair_changes_are_retained():
    result = compare_pairs([("a", "b")], [("c", "a")])
    assert not result["identical"]
    assert result["removed"] == [["a", "b"]]
    assert result["added"] == [["a", "c"]]


def test_artifacts_compare_bytes_and_hashes_not_absolute_paths():
    before = {"x": {"path": "/old/x", "bytes": 2, "sha256": "aa"}}
    after = {"x": {"path": "/new/x", "bytes": 2, "sha256": "aa"}}
    assert compare_inventory(before, after)["identical"]
    after["x"]["sha256"] = "bb"
    assert compare_inventory(before, after)["changed"] == ["x"]
    assert compare_inventory(before, {})["missing"] == ["x"]
    assert compare_inventory({}, after)["extra"] == ["x"]


def test_missing_or_duplicate_of_tree_rejected(tmp_path):
    with pytest.raises(ValueError, match="exactly one"):
        native_tree("orthofinder_full", tmp_path)
    for name in ("a", "b"):
        path = tmp_path / name / "Species_Tree/SpeciesTree_rooted.txt"
        path.parent.mkdir(parents=True)
        path.write_text("((a,b),(c,d));")
    with pytest.raises(ValueError, match="exactly one"):
        native_tree("orthofinder_full", tmp_path)


def test_absent_mandatory_artifacts_rejected(tmp_path):
    with pytest.raises(ValueError, match="Incomplete"):
        artifact_inventory("orthohmm_satellite_v2", tmp_path)


def test_no_ambiguous_families_is_valid_inventory(tmp_path):
    path = tmp_path / "orthohmm_working_res/phylogeny_candidate_superfamilies.txt"
    path.parent.mkdir()
    path.write_text("a b c d\n")
    assert len(artifact_inventory("orthohmm_satellite_v2", tmp_path)) == 1


def test_existing_output_refused(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(FileExistsError):
        run(tmp_path, "baseline_20261101", tmp_path)


@pytest.mark.parametrize("records", [[], [{"condition": "baseline", "seed": 1,
    "method": "orthohmm_satellite_v2", "status": "failed"}]])
def test_missing_or_failed_baseline_cannot_be_control(records):
    with pytest.raises(ValueError, match="completed inferred baseline"):
        baseline_records({"condition": "baseline", "seed": 1}, {"records": records}, {}, {})
