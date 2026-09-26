from copy import deepcopy
import json

import pytest

from benchmark_tools.matched_graph_worker import run, validate_numeric, write_partition


DATA = dict(gene_names=["a", "b", "c"], gene_to_species=[0, 1, 1],
            hit_queries=[0, 1], hit_targets=[1, 0], hit_scores=[2.5, 3.5])


def test_valid_numeric_and_empty_hits():
    validate_numeric(DATA)
    value = deepcopy(DATA)
    for key in ("hit_queries", "hit_targets", "hit_scores"):
        value[key] = []
    validate_numeric(value)


@pytest.mark.parametrize("key,value", [
    ("gene_names", ["a", "a", "c"]), ("gene_names", ["b", "a", "c"]),
    ("gene_names", ["a a", "b", "c"]), ("gene_to_species", [0, 2, 2]),
    ("gene_to_species", [0, 1]), ("gene_to_species", [0, True, 1]),
    ("hit_queries", [-1, 1]), ("hit_targets", [3, 0]),
    ("hit_queries", [0.0, 1]), ("hit_scores", [float("nan"), 1]),
    ("hit_scores", [0, 1]), ("hit_scores", [True, 1]), ("hit_scores", [1]),
])
def test_invalid_numeric(key, value):
    data = deepcopy(DATA)
    data[key] = value
    with pytest.raises(ValueError):
        validate_numeric(data)


def test_duplicate_hits_and_unexpected_labels():
    data = deepcopy(DATA)
    data["hit_queries"] = [0, 0]
    data["hit_targets"] = [1, 1]
    with pytest.raises(ValueError):
        validate_numeric(data)
    data = dict(DATA, truth={})
    with pytest.raises(ValueError):
        validate_numeric(data)


def test_partition_has_complete_gene_universe(tmp_path):
    path = tmp_path / "partition.tsv"
    write_partition(path, [[2], [1, 0]], DATA["gene_names"])
    assert path.read_text() == "a\tb\nc\n"
    with pytest.raises(FileExistsError):
        write_partition(path, [[0], [1], [2]], DATA["gene_names"])


@pytest.mark.parametrize("groups", [[[0, 1]], [[0, 1], [1, 2]], [[], [0, 1, 2]], [[0, 1, 3]]])
def test_invalid_partition(tmp_path, groups):
    with pytest.raises(ValueError):
        write_partition(tmp_path / "partition", groups, DATA["gene_names"])


def test_nonisolated_rejected(tmp_path):
    with pytest.raises(RuntimeError, match="-I"):
        run(tmp_path / "numeric", tmp_path / "output")
