from pathlib import Path

import pytest

from benchmark_tools.validate_threadripper_outputs import orthofinder_mapping


def fixture(tmp_path):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    (inputs / "a.fa").write_text(">a1 annotation\nACDE\n>a2\nFGHI\n")
    (inputs / "b.fa").write_text(">b1\nACDE\n")
    output = tmp_path / "output"
    output.mkdir()
    (output / "SpeciesIDs.txt").write_text("0: a.fa\n1: b.fa\n")
    (output / "SequenceIDs.txt").write_text("0_0: a1 annotation\n0_1: a2\n1_0: b1\n")
    return dict(configuration=dict(output=str(output)), prepared_input_directory=str(inputs),
                expected_native_order=["a.fa", "b.fa"],
                dataset=dict(inputs=[dict(path=str(inputs / n)) for n in ("a.fa", "b.fa")]))


def test_exact_native_numbering(tmp_path):
    result = orthofinder_mapping(fixture(tmp_path))
    assert result["species"] == 2 and result["sequences"] == 3
    assert len(result["checked_files"]) == 2


@pytest.mark.parametrize("change", ["species", "sequence", "missing", "duplicate", "empty", "extra", "symlink"])
def test_rejects_native_mapping_drift(tmp_path, change):
    run = fixture(tmp_path)
    output = Path(run["configuration"]["output"])
    path = output / "SequenceIDs.txt"
    if change == "species":
        (output / "SpeciesIDs.txt").write_text("0: b.fa\n1: a.fa\n")
    elif change == "sequence":
        path.write_text("0_0: a2\n0_1: a1\n1_0: b1\n")
    elif change == "missing":
        path.write_text("0_0: a1\n1_0: b1\n")
    elif change == "duplicate":
        path.write_text(path.read_text() + "0_0: a1\n")
    elif change == "empty":
        path.write_text("0_0: \n")
    elif change == "extra":
        path.write_text(path.read_text() + "1_1: extra\n")
    else:
        saved = tmp_path / "saved"
        path.rename(saved)
        path.symlink_to(saved)
    with pytest.raises(ValueError):
        orthofinder_mapping(run)
