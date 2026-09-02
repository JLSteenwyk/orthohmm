from benchmark_tools.create_orthobench_partition import partition_names


def test_partition_names_is_balanced_deterministic_and_disjoint():
    names = [f"RefOG{index:03d}.txt" for index in range(1, 8)]

    development, validation = partition_names(names)
    repeated = partition_names(list(reversed(names)) + [names[0]])

    assert (development, validation) == repeated
    assert len(development) == 4
    assert len(validation) == 3
    assert set(development).isdisjoint(validation)
    assert set(development) | set(validation) == set(names)
