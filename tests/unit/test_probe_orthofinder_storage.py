import pytest

from benchmark_tools.probe_orthofinder_storage import record, validate


def fixture(tmp_path):
    inputs, output = tmp_path / "input", tmp_path / "out"
    inputs.mkdir()
    (inputs / "a.fa").write_text(">a\nACDEFG\n")
    before = [record(inputs / "a.fa")]
    working = output / "Results_test" / "WorkingDirectory"
    working.mkdir(parents=True)
    (working / "SpeciesIDs.txt").write_text("0: a.fa\n")
    (working / "Species0.fa").write_text(">0_0\nACDEFG\n")
    return inputs, output, before, working


def test_valid(tmp_path):
    inputs, output, before, _ = fixture(tmp_path)
    assert validate(inputs, output, before)["native_order"] == ["a.fa"]


@pytest.mark.parametrize("change", ["input_write", "mapping", "missing", "symlink"])
def test_rejects_bad_placement(tmp_path, change):
    inputs, output, before, working = fixture(tmp_path)
    if change == "input_write":
        (inputs / "generated").write_text("unexpected")
    elif change == "mapping":
        (working / "SpeciesIDs.txt").write_text("0: b.fa\n")
    elif change == "missing":
        (working / "Species0.fa").unlink()
    else:
        (output / "link").symlink_to(inputs)
    with pytest.raises(ValueError):
        validate(inputs, output, before)
