import pytest

from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.prepare_scaling_inputs import planned_runs
from benchmark_tools.prepare_scaling_transfer import portable_inputs, copy_verified, prepare


def panel():
    rows = [{"input": {"path": "/original/p%d.fa" % n, "bytes": 10, "sha256": "a" * 64},
             "proteins": 1, "sequence_characters": 3} for n in range(12)]
    return {"ordered_proteomes": rows, "planned_runs": planned_runs(), "datasets": [
        {"proteomes": n, "inputs": [r["input"] for r in rows[:n]], "proteins": n,
         "sequence_characters": n * 3} for n in (4, 8, 12)]}


def test_portable_nested_inputs_preserve_inventory():
    files, datasets = portable_inputs(panel())
    assert len(files) == 12
    assert all(r["path"].startswith("inputs/") and not r["path"].startswith("/") for r in files)
    assert [len(d["inputs"]) for d in datasets] == [4, 8, 12]
    assert datasets[0]["inputs"] == datasets[2]["inputs"][:4]


@pytest.mark.parametrize("problem", ["order", "duplicate", "nested", "counts", "size"])
def test_changed_transfer_plan_rejected(problem):
    data = panel()
    if problem == "order":
        data["planned_runs"].reverse()
    elif problem == "duplicate":
        data["ordered_proteomes"][-1] = data["ordered_proteomes"][0]
    elif problem == "nested":
        data["datasets"][0]["inputs"].pop()
    elif problem == "counts":
        data["datasets"][0]["proteins"] += 1
    else:
        data["datasets"].pop()
    with pytest.raises(ValueError):
        portable_inputs(data)


def test_exact_copy_and_no_overwrite(tmp_path):
    source, target = tmp_path / "source", tmp_path / "nested/target"
    source.write_bytes(b"fixture")
    expected = record(source)
    copy_verified(source, target, expected)
    assert target.read_bytes() == b"fixture" and not target.is_symlink()
    with pytest.raises(FileExistsError):
        copy_verified(source, target, expected)


def test_changed_source_rejected(tmp_path):
    source, target = tmp_path / "source", tmp_path / "target"
    source.write_bytes(b"before")
    expected = record(source)
    source.write_bytes(b"after")
    with pytest.raises(ValueError):
        copy_verified(source, target, expected)
    assert not target.exists()


def test_existing_bundle_not_overwritten(tmp_path):
    with pytest.raises(FileExistsError):
        prepare(tmp_path, tmp_path)
