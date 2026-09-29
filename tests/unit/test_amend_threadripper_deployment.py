from copy import deepcopy

import pytest

from benchmark_tools.amend_threadripper_deployment import amend_plan


def panel():
    methods = ("orthohmm_high_sensitivity", "orthohmm_satellite_v2", "orthofinder_full")
    rows = []
    for repeat in range(3):
        for size_index, size in enumerate((4, 8, 12)):
            for offset in range(3):
                method = methods[(repeat + size_index + offset) % 3]
                rows.append(dict(index=len(rows), proteomes=size, repeat=repeat, native_method=method,
                                 native_argv=["old", "-m", method], configuration=dict(argv=["old", "adapter"]),
                                 original_native_argv=["historical"], input_checksum="unchanged"))
    return dict(runs=rows, resources={"cpu": 32}, other=["unchanged"])


def test_only_two_interpreter_fields_change_per_orthohmm_row():
    original = panel()
    saved = deepcopy(original)
    amended, indices = amend_plan(original, "old", "new")
    assert original == saved
    assert len(indices) == 18
    for i in indices:
        assert amended["runs"][i]["native_argv"][0] == "new"
        assert amended["runs"][i]["configuration"]["argv"][0] == "new"
        amended["runs"][i]["native_argv"][0] = "old"
        amended["runs"][i]["configuration"]["argv"][0] = "old"
    assert amended == original


@pytest.mark.parametrize("field,value", [("index", 9), ("repeat", 5), ("proteomes", 16),
                                         ("native_method", "unknown")])
def test_changed_identity_rejected(field, value):
    original = panel()
    original["runs"][0][field] = value
    with pytest.raises(ValueError, match="identities"):
        amend_plan(original, "old", "new")


def test_wrong_interpreter_rejected():
    with pytest.raises(ValueError, match="interpreter"):
        amend_plan(panel(), "different", "new")
