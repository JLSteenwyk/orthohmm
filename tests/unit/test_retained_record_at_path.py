from copy import deepcopy

import pytest

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def test_only_path_changes_and_original_is_preserved(tmp_path, retained_record_at_path):
    target = tmp_path / "explicit retained file $literal's.txt"
    target.write_bytes(b"retained bytes\n")
    original = {**record(target), "path": "/historical/checkout/file.txt"}
    saved = deepcopy(original)
    bound = retained_record_at_path(original, target)
    assert bound == record(target)
    assert original == saved
    assert bound is not original
    check(bound)


@pytest.mark.parametrize("problem", ["changed", "missing"])
def test_binding_does_not_repin_or_hide_bad_payload(tmp_path, problem, retained_record_at_path):
    target = tmp_path / "file.txt"
    target.write_bytes(b"retained bytes\n")
    original = {**record(target), "path": "/historical/checkout/file.txt"}
    if problem == "changed":
        target.write_bytes(b"different bytes\n")
    else:
        target.unlink()
    bound = retained_record_at_path(original, target)
    assert {key: value for key, value in bound.items() if key != "path"} == {
        key: value for key, value in original.items() if key != "path"}
    with pytest.raises(ValueError if problem == "changed" else FileNotFoundError):
        check(bound)
