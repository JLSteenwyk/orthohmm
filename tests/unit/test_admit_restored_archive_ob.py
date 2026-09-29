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
    native = {str(i): dict(byte_equal=True) for i in range(4)}
    native["0"]["byte_equal"] = native_equal
    assert module.reproduction_equal(compared, native) is expected


def test_missing_native_comparisons_rejected():
    with pytest.raises(ValueError, match="four"):
        module.reproduction_equal({}, {})
