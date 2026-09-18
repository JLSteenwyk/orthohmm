import sys

import numpy as np
import pytest

from benchmark_tools.estimate_rbnh_array_payload import bounds, estimate
from benchmark_tools import estimate_rbnh_array_payload as model
from orthohmm.accuracy import build_rbnh_edges


@pytest.mark.parametrize("seed", range(8))
def test_bounds_contain_actual_named_allocations(seed):
    rng = np.random.default_rng(seed)
    genes, hits = 23, 200
    species = (np.arange(genes) % 5).astype(np.int32)
    q = rng.integers(0, genes, size=hits, dtype=np.int32)
    t = rng.integers(0, genes, size=hits, dtype=np.int32)
    scores = rng.integers(1, 4, size=hits).astype(np.float64)
    result = bounds(genes, 5, hits, int(np.count_nonzero(q == t)))
    fields = set(result["fixed_named_array_bytes"]) | {"winning_rows", "winning_keys", "first", "selected_rows"}
    snapshots = []

    def trace(frame, event, arg):
        if (event == "line" and frame.f_code is build_rbnh_edges.__code__
                and not snapshots and fields <= frame.f_locals.keys()):
            snapshots.append({k: frame.f_locals[k].nbytes for k in fields})
        return trace

    previous = sys.gettrace()
    try:
        sys.settrace(trace)
        build_rbnh_edges([str(i) for i in range(genes)], species, q, t, scores)
    finally:
        sys.settrace(previous)
    assert len(snapshots) == 1
    actual = snapshots[0]
    for name, value in result["fixed_named_array_bytes"].items():
        assert actual[name] == value
    assert result["named_array_payload_lower_bytes"] <= sum(actual.values())
    assert sum(actual.values()) <= result["named_array_payload_upper_bytes_at_this_snapshot_only"]
    assert result["input_array_logical_bytes_separate"] == sum(a.nbytes for a in (species, q, t, scores))


@pytest.mark.parametrize("hits", [0, 10])
def test_no_nonself_hits_do_not_reach_snapshot(hits):
    result = bounds(10, 2, hits, hits)
    assert not result["snapshot_reached"]
    assert result["named_array_payload_lower_bytes"] == 0
    assert result["named_array_payload_upper_bytes_at_this_snapshot_only"] == 0
    assert not result["total_peak_ram_bound_available"]
    assert not result["graph_feasibility_admitted"]


def test_large_dimensions_use_exact_integer_arithmetic():
    result = bounds(984137, 78, 10**10, 984137)
    assert result["gene_species_slots"] == 76762686
    assert result["fixed_named_array_bytes"]["keys"] == 8 * (10**10 - 984137)
    assert result["named_array_payload_lower_bytes"] > 192 * 1024**3


@pytest.mark.parametrize("args", [(0, 2, 0, 0), (10, 0, 0, 0), (10, 2, -1, 0),
                                  (10, 2, 10, 11), (10, 2, 10, -1), (True, 2, 0, 0),
                                  (10, 2, 10.0, 0)])
def test_invalid_dimensions(args):
    with pytest.raises(ValueError):
        bounds(*args)


def test_wrong_core_rejected_before_checkpoint_access(tmp_path):
    core = tmp_path / "accuracy.py"
    core.write_text("# Not frozen source\n")
    with pytest.raises(ValueError, match="frozen accuracy.py"):
        estimate(tmp_path / "missing", "0" * 64, core)


def test_checkpoint_changed_after_audit_is_rejected(tmp_path, monkeypatch):
    core = tmp_path / "accuracy.py"
    core.write_text("# Fixture\n")
    monkeypatch.setattr(model, "CORE_SHA", model.record(core)["sha256"])
    data = tmp_path / "gene_to_species.npy"
    np.save(data, np.array([0, 1], dtype=np.int32))
    expected = model.record(data)
    np.save(data, np.array([0, 2], dtype=np.int32))
    monkeypatch.setattr(model, "audit", lambda *args: {})
    monkeypatch.setattr(model, "read_frozen", lambda *args: {"files": {data.name: expected}})
    with pytest.raises(ValueError, match="changed after numeric audit"):
        estimate(tmp_path, "0" * 64, core)
