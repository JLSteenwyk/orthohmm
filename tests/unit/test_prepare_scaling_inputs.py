from collections import Counter

import pytest

from benchmark_tools.prepare_scaling_inputs import ordered_inputs, planned_runs, prepare, METHODS, SIZES


def inputs(prefix="/first"):
    return [{"path": f"{prefix}/species_{i}.fa", "sha256": str(i), "bytes": i + 1} for i in range(12)]


def test_order_is_input_order_and_location_independent():
    first = ordered_inputs(inputs())
    assert ordered_inputs(list(reversed(inputs()))) == first
    assert [r["sha256"] for r in ordered_inputs(inputs("/second"))] == [r["sha256"] for r in first]
    assert set(r["path"] for r in first[:4]) < set(r["path"] for r in first[:8]) < set(r["path"] for r in first)


@pytest.mark.parametrize("problem", ["missing", "duplicate_path", "duplicate_basename"])
def test_reject_ambiguous_inputs(problem):
    rows = inputs()
    if problem == "missing":
        rows.pop()
    elif problem == "duplicate_path":
        rows[0] = rows[1]
    else:
        rows[0]["path"] = "/elsewhere/species_1.fa"
    with pytest.raises(ValueError):
        ordered_inputs(rows)


def test_complete_balanced_run_plan():
    runs = planned_runs()
    assert [r["index"] for r in runs] == list(range(27))
    assert Counter((r["proteomes"], r["method"]) for r in runs) == Counter({(s, m): 3 for s in SIZES for m in METHODS})
    for size in SIZES:
        positions = Counter()
        for repeat in range(3):
            block = [r for r in runs if r["proteomes"] == size and r["repeat"] == repeat]
            assert {r["method"] for r in block} == set(METHODS)
            positions.update((r["method"], i) for i, r in enumerate(block))
        assert positions == Counter({(m, i): 1 for m in METHODS for i in range(3)})


def test_refuses_existing_output(tmp_path):
    with pytest.raises(FileExistsError):
        prepare(tmp_path, tmp_path)
