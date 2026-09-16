import pytest
from benchmark_tools.audit_qfo_replay_inputs import verify_species_partition


def test_species_mapping_is_a_partition_not_order_assumption():
    assert verify_species_partition(["b", "a"], [9, 4], [("first", ["a"]), ("second", ["b"])]) == {"first": 4, "second": 9}


@pytest.mark.parametrize("records", [[("first", ["a"])], [("first", ["a", "x"])],
                                    [("first", ["a", "b"])], [("first", ["a"]), ("second", ["a", "b"])],
                                    [("first", []), ("second", ["a", "b"])]])
def test_rejects_wrong_fasta_to_checkpoint_mapping(records):
    with pytest.raises(ValueError):
        verify_species_partition(["a", "b"], [0, 1], records)
