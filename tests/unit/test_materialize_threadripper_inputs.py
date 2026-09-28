from copy import deepcopy

import pytest

from benchmark_tools import materialize_threadripper_inputs as module


def fixture(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    rows = []
    for i in range(12):
        path = raw / f"species_{i}.fa"
        path.write_text(f">gene_{i}\nACDEFG\n")
        rows.append(module.record(path))
    inputs = {"datasets": [dict(proteomes=n, inputs=rows[:n]) for n in (4, 8, 12)]}
    reference = {"datasets": [dict(proteomes=n,
        native_order=[f"species_{i}.fa" for i in reversed(range(n))],
        inputs_in_native_order=list(reversed(rows[:n]))) for n in (4, 8, 12)]}
    return inputs, reference


@pytest.mark.parametrize("matching", [True, False])
@pytest.mark.parametrize("reverse_creation", [True, False])
def test_copy_and_actual_order_check(tmp_path, monkeypatch, matching, reverse_creation):
    inputs, reference = fixture(tmp_path)
    def snapshot(core, copied, runtime):
        for dataset in copied["datasets"]:
            assert all(module.record(r["path"]) == r for r in dataset["inputs"])
            n = dataset["proteomes"]
            indices = range(n) if reverse_creation else reversed(range(n))
            assert [module.Path(r["path"]).name for r in dataset["inputs"]] == [f"species_{i}.fa" for i in indices]
        datasets = deepcopy(reference["datasets"])
        if not matching:
            datasets[0]["native_order"].reverse()
        return {"datasets": datasets}
    monkeypatch.setattr(module, "snapshot", snapshot)
    result = module.materialize(inputs, reference, {"core_root": "unused", "core_sources": []},
                                tmp_path / "out", reverse_creation=reverse_creation)
    assert result["status"] == ("native_order_matched" if matching else "native_order_mismatched")
    assert result["scientific_execution_authorized"] is False
    assert (tmp_path / "out" / "manifest.json").is_file()


@pytest.mark.parametrize("change", ["bytes", "duplicate", "order", "sizes", "existing"])
def test_reject_invalid_inputs_before_copy(tmp_path, change):
    inputs, reference = fixture(tmp_path)
    if change == "bytes":
        inputs["datasets"][0]["inputs"][0]["sha256"] = "bad"
    elif change == "duplicate":
        reference["datasets"][0]["native_order"][0] = reference["datasets"][0]["native_order"][1]
    elif change == "order":
        reference["datasets"][0]["inputs_in_native_order"].reverse()
    elif change == "sizes":
        inputs["datasets"].pop()
    else:
        (tmp_path / "out").mkdir()
    with pytest.raises((ValueError, FileExistsError)):
        module.materialize(inputs, reference, {}, tmp_path / "out")
