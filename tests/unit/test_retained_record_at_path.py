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


@pytest.mark.parametrize("problem", [
    "changed", "missing", "same_size", "truncated", "wrong_size", "wrong_sha",
])
def test_binding_does_not_repin_or_hide_bad_payload(tmp_path, problem, retained_record_at_path):
    target = tmp_path / "file.txt"
    target.write_bytes(b"retained bytes\n")
    original = {**record(target), "path": "/historical/checkout/file.txt"}
    if problem == "changed":
        target.write_bytes(b"different bytes\n")
    elif problem == "missing":
        target.unlink()
    elif problem == "same_size":
        target.write_bytes(b"X" + target.read_bytes()[1:])
        assert target.stat().st_size == original["bytes"]
    elif problem == "truncated":
        target.write_bytes(target.read_bytes()[:-1])
    elif problem == "wrong_size":
        original["bytes"] += 1
    elif problem == "wrong_sha":
        original["sha256"] = "0" * 64
    bound = retained_record_at_path(original, target)
    assert {key: value for key, value in bound.items() if key != "path"} == {
        key: value for key, value in original.items() if key != "path"}
    with pytest.raises(FileNotFoundError if problem == "missing" else ValueError):
        check(bound)
