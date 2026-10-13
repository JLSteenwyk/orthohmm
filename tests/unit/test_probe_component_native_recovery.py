import copy

import pytest

from benchmark_tools.probe_component_native_recovery import validate_native_receipt


def receipt():
    return {"status": "complete", "command": ["python", "native"], "git": {"commit": "frozen", "dirty": False}}


def test_native_receipt_requires_matching_clean_source_and_planned_command():
    item = receipt()
    before = copy.deepcopy(item)
    validate_native_receipt(item, ["python", "native"], "frozen", "\tExit status: 0\n")
    assert item == before


@pytest.mark.parametrize("mutation", ["status", "command", "revision", "dirty"])
def test_mismatched_native_receipt_is_rejected(mutation):
    item = receipt()
    if mutation in ("status", "command"):
        item[mutation] = "wrong"
    elif mutation == "revision":
        item["git"]["commit"] = "wrong"
    else:
        item["git"]["dirty"] = True
    with pytest.raises(ValueError):
        validate_native_receipt(item, ["python", "native"], "frozen", "Exit status: 0\n")


@pytest.mark.parametrize("text", ["", "Exit status: 1\n", "Exit status: 0\nExit status: 0\n"])
def test_missing_nonzero_or_duplicate_exit_status_is_rejected(text):
    with pytest.raises(ValueError, match="unambiguous"):
        validate_native_receipt(receipt(), ["python", "native"], "frozen", text)
