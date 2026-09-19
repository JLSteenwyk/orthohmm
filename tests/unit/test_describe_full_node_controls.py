import copy

import pytest

from benchmark_tools import describe_full_node_controls as module


def rows():
    return [dict(index=i, block=block, mode=mode, status="validated",
        common=dict(distributions=dict(metric=dict(median={"steady": 5, "churn": 2, "contended": 8}[mode]))))
        for i, (block, mode) in enumerate((b, m) for b, modes in enumerate(module.ORDER) for m in modes)]


def test_inventory_and_signed_block_differences():
    data = rows()
    module.inventory_check(dict(trials=data, validated_trials=9))
    result = module.block_differences(data)
    for block in result:
        assert block["common_interval_median_differences"] == {
            "churn_minus_steady": {"metric": -3}, "contended_minus_steady": {"metric": 3}}


@pytest.mark.parametrize("fault", ["missing", "order", "invalid", "duplicate"])
def test_wrong_inventory_rejected(fault):
    data = rows()
    if fault == "missing":
        data.pop()
    elif fault == "order":
        data.reverse()
    elif fault == "invalid":
        data[0]["status"] = "failed"
    else:
        data[1] = copy.deepcopy(data[0])
    with pytest.raises(ValueError):
        module.inventory_check(dict(trials=data, validated_trials=9))


@pytest.mark.parametrize("fault", ["missing", "duplicate", "metrics"])
def test_incomplete_contrast_rejected(fault):
    data = rows()
    if fault == "missing":
        data.pop()
    elif fault == "duplicate":
        data.append(copy.deepcopy(data[0]))
    else:
        data[0]["common"]["distributions"] = {}
    with pytest.raises(ValueError):
        module.block_differences(data)


def test_wrong_audit_pin_rejected(tmp_path):
    path = tmp_path / "audit.gz"
    path.write_bytes(b"wrong")
    with pytest.raises(ValueError, match="pinned"):
        module.report(path)
