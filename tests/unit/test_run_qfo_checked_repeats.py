import copy
import hashlib

import numpy as np
import pytest

from benchmark_tools.run_qfo_checked_repeats import check_gate, constructor_digest, compare_partition, run


def fixture_gate():
    saved = {"edges": 3, "vertices": 4, "directed": False, "endpoint_sha256": "endpoints", "weight_sha256": "weights"}
    arguments = {"initial_membership": None, "weights": "weight", "n_iterations": 2,
                 "max_comm_size": 0, "seed": 4, "kwargs": {"resolution_parameter": .1},
                 "partition_type": "leidenalg.VertexPartition.CPMVertexPartition"}
    boundary = {"accuracy_evaluated": False, "calls": [{"status": "optimizer_returned",
                "arguments": arguments, "before": copy.deepcopy(saved), "after": copy.deepcopy(saved), "saved": copy.deepcopy(saved)}]}
    adapter = {"format": "python_pairs", "accuracy_evaluated": False, "calls": [{"status": "constructor_returned",
               "dtype": "int32", "shape": [3, 2], "n": 4, "directed": False, "ordered_input_bytes_sha256": "input"}]}
    return boundary, adapter, saved


def test_gate_accepts_exact_observation():
    check_gate(*fixture_gate(), "input")


@pytest.mark.parametrize("change", ["before", "after", "saved", "seed", "resolution", "format", "dtype", "digest", "count", "accuracy"])
def test_gate_rejects_changes(change):
    boundary, adapter, saved = fixture_gate()
    if change in ("before", "after", "saved"):
        boundary["calls"][0][change]["endpoint_sha256"] = "altered"
    elif change == "seed":
        boundary["calls"][0]["arguments"]["seed"] = 5
    elif change == "resolution":
        boundary["calls"][0]["arguments"]["kwargs"]["resolution_parameter"] = .2
    elif change == "format":
        adapter["format"] = "numpy"
    elif change == "dtype":
        adapter["calls"][0]["dtype"] = "int64"
    elif change == "digest":
        adapter["calls"][0]["ordered_input_bytes_sha256"] = "changed"
    elif change == "count":
        adapter["calls"].append(copy.deepcopy(adapter["calls"][0]))
    else:
        boundary["accuracy_evaluated"] = True
    with pytest.raises(ValueError):
        check_gate(boundary, adapter, saved, "input")


def test_digest_preserves_orientation_and_chunk_boundaries(tmp_path):
    sources = np.arange(100003, dtype=np.int32)
    targets = sources[::-1].copy()
    np.save(tmp_path / "sources.npy", sources)
    np.save(tmp_path / "targets.npy", targets)
    assert constructor_digest(tmp_path) == hashlib.sha256(np.column_stack((sources, targets)).tobytes()).hexdigest()


@pytest.mark.parametrize("content", ["a\n", "a a\nb\n", "a\na b\n", "a b c\n"])
def test_self_comparison_requires_full_valid_partition(tmp_path, content):
    partition = tmp_path / "partition.txt"
    partition.write_text(content)
    with pytest.raises(ValueError):
        compare_partition(partition, partition, {"a", "b"})


def test_existing_output_never_reused(tmp_path):
    with pytest.raises(FileExistsError):
        run(tmp_path, tmp_path)
