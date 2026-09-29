import pytest

from benchmark_tools import admit_restored_archive_ob as module


def test_existing_output_refused(tmp_path):
    with pytest.raises(FileExistsError):
        module.audit(tmp_path / "unused", tmp_path)


def test_unfinished_job_creates_no_output(tmp_path, monkeypatch):
    def reject(_):
        raise ValueError("Not completed")
    monkeypatch.setattr(module, "verify", reject)
    output = tmp_path / "audit"
    with pytest.raises(ValueError, match="Not completed"):
        module.audit(tmp_path, output)
    assert not output.exists()


def test_changed_baseline_rejected(tmp_path, monkeypatch):
    monkeypatch.setattr(module, "RESULTS", tmp_path)
    (tmp_path / "reconstructed_full_ob_result_22376.json").write_text("{}")
    with pytest.raises(ValueError, match="Changed admitted baseline"):
        module.baseline()


@pytest.mark.parametrize("partition_equal,score_equal,native_equal,expected", [
    (True, True, True, True), (False, True, True, False),
    (True, False, True, False), (True, True, False, False),
])
def test_equality_requires_all_dimensions(partition_equal, score_equal, native_equal, expected):
    compared = dict(partitions=dict(label_invariant_equal=partition_equal), score_objects_equal=score_equal)
    native = {name: dict(byte_equal=True) for name in module.NATIVE_FILES}
    native[module.NATIVE_FILES[0]]["byte_equal"] = native_equal
    assert module.reproduction_equal(compared, native) is expected


def test_missing_native_comparisons_rejected():
    with pytest.raises(ValueError, match="four"):
        module.reproduction_equal({}, {})


def test_four_wrong_native_names_rejected():
    with pytest.raises(ValueError, match="four"):
        module.reproduction_equal({}, {str(i): dict(byte_equal=True) for i in range(4)})


def test_equal_refog_scores_do_not_hide_unscored_partition_change():
    original = [frozenset({"a", "b"}), frozenset({"x", "y"})]
    current = [frozenset({"a", "b"}), frozenset({"x"}), frozenset({"y"})]
    refs = {"ref1": {"a", "b"}}
    old_score = module.score_partition(original, refs, {})
    new_score = module.score_partition(current, refs, {})
    assert old_score == new_score
    compared = module.comparison(original, current, old_score, new_score)
    native = {name: dict(byte_equal=True) for name in module.NATIVE_FILES}
    assert not module.reproduction_equal(compared, native)
    assert compared["partitions"]["genes_in_changed_groups"] == 2


def test_group_order_does_not_change_reproduction_equality():
    original = [frozenset({"a", "b"}), frozenset({"x", "y"})]
    current = list(reversed(original))
    score = module.score_partition(original, {"ref1": {"a", "b"}}, {})
    compared = module.comparison(original, current, score, score)
    native = {name: dict(byte_equal=True) for name in module.NATIVE_FILES}
    assert module.reproduction_equal(compared, native)
