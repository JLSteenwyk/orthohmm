import pytest

from benchmark_tools.merge_candidate_superfamilies import (
    validate_partition,
    write_clusters,
)


def test_validate_partition_requires_each_gene_once():
    validate_partition([[0, 2], [1]], 3)

    with pytest.raises(ValueError, match="exactly once"):
        validate_partition([[0, 1], [1]], 3)


def test_write_clusters_is_deterministic(tmp_path):
    output = tmp_path / "clusters.txt"

    write_clusters(output, [[2, 0], [3], [1]], ["a", "b", "c", "d"])

    assert output.read_text().splitlines() == ["a c", "b", "d"]
