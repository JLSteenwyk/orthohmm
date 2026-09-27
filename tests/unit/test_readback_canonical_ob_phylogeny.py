import pytest

from benchmark_tools import readback_canonical_ob_phylogeny as reader
from benchmark_tools.score_orthobench_partition import score_partition


def test_equal_scores_do_not_imply_identical_partitions():
    original = [frozenset("ab"), frozenset("cd")]
    current = [frozenset("ab"), frozenset("c"), frozenset("d")]
    references = {"R1": set("ab")}
    before = score_partition(original, references, {})
    after = score_partition(current, references, {})
    result = reader.comparison(original, current, before, after)
    assert result["score_objects_equal"]
    assert result["changed_reference_families"] == []
    assert not result["partitions"]["label_invariant_equal"]
    assert result["partitions"]["genes_in_changed_groups"] == 2


def test_equal_aggregate_retains_opposing_family_changes():
    references = {"R1": set("ab"), "R2": set("cd")}
    original = [frozenset("ab"), frozenset("c"), frozenset("d")]
    current = [frozenset("a"), frozenset("b"), frozenset("cd")]
    result = reader.comparison(original, current,
        score_partition(original, references, {}), score_partition(current, references, {}))
    assert all(value == 0 for value in result["score_differences_percentage_points"].values())
    assert [r["refog"] for r in result["changed_reference_families"]] == ["R1", "R2"]
    assert not result["score_objects_equal"]


def test_label_order_does_not_change_partition_comparison():
    groups = [frozenset("ab"), frozenset("cd")]
    score = score_partition(groups, {"R1": set("ab")}, {})
    result = reader.comparison(groups, groups[::-1], score, score)
    assert result["partitions"]["label_invariant_equal"]


@pytest.mark.parametrize("duplicate", [False, True])
def test_invalid_family_inventory_rejected(duplicate):
    groups = [frozenset("ab")]
    score = score_partition(groups, {"R1": set("ab")}, {})
    changed = dict(score, refog_records=score["refog_records"] * (2 if duplicate else 0))
    with pytest.raises(ValueError, match="inventories"):
        reader.comparison(groups, groups, score, changed)


def test_live_or_failed_native_job_blocks_workflow_before_output(tmp_path, monkeypatch):
    def reject(*args):
        raise ValueError("Native job not successfully terminal")
    monkeypatch.setattr(reader.admission, "audit", reject)
    output = tmp_path / "readback"
    with pytest.raises(ValueError, match="terminal"):
        reader.readback(tmp_path, tmp_path, 123, output)
    assert not output.exists()


def test_existing_readback_is_never_overwritten(tmp_path):
    with pytest.raises(FileExistsError):
        reader.readback(tmp_path, tmp_path, 123, tmp_path)


def test_save_never_overwrites_partial_failure_evidence(tmp_path):
    path = tmp_path / "stage.json"
    reader.save(path, {"first_attempt": True})
    before = path.read_bytes()
    with pytest.raises(FileExistsError):
        reader.save(path, {"replacement": True})
    assert path.read_bytes() == before
