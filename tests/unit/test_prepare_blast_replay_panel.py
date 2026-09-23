import pytest

from benchmark_tools.prepare_blast_replay_panel import panel_ids, FIRST, FAILED, SELENO, BOUNDARY


def test_five_roles_keep_input_order():
    assert panel_ids([FIRST, "unused", FAILED, SELENO, "previous", BOUNDARY, "later"]) == [
        FIRST, FAILED, SELENO, "previous", BOUNDARY]


@pytest.mark.parametrize("ids", [[], ["wrong", FIRST, FAILED, SELENO, "previous", BOUNDARY],
    [FIRST, FAILED, SELENO, "previous", BOUNDARY, FIRST],
    [FIRST, FAILED, "previous", BOUNDARY], [FIRST, FAILED, SELENO, BOUNDARY]])
def test_invalid_panel_rejected(ids):
    with pytest.raises(ValueError):
        panel_ids(ids)
