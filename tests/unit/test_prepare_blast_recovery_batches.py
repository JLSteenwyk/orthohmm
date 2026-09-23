import pytest

from benchmark_tools.prepare_blast_recovery_batches import batch_members


@pytest.mark.parametrize("count,sizes", [(1, [1]), (5000, [5000]), (5001, [5000, 1]),
                                        (98913, [5000]*19 + [3913])])
def test_batches_preserve_exhaustive_order(count, sizes):
    ids = [f"query_{i}" for i in range(count)]
    batches = batch_members(ids)
    assert [len(b) for b in batches] == sizes
    assert [q for b in batches for q in b] == ids


@pytest.mark.parametrize("ids", [[], ["a", "a"]])
def test_invalid_inventory(ids):
    with pytest.raises(ValueError):
        batch_members(ids)
