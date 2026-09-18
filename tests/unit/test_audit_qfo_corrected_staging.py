import pytest

from benchmark_tools.audit_qfo_corrected_staging import inventory, record


def setup(tmp_path):
    (tmp_path / "staging_manifest.json").write_text("{}")
    (tmp_path / "a.fasta").write_text(">sp|A|A_A\nMA\n>tr|B|B_A\nMC\n")
    (tmp_path / "b.fasta").write_text(">sp|C|C_B\nMUZ\n")
    mapping = {"species": ["A", "B"], "Goff": [0, 2, 3], "mapping": {"A": 1, "B": 2, "C": 3}}
    return mapping


def run(tmp_path, mapping):
    return inventory(tmp_path, [record(tmp_path / name) for name in ("a.fasta", "b.fasta")], mapping)


def test_complete_species_intervals_and_no_normalization(tmp_path):
    mapping = setup(tmp_path)
    result = run(tmp_path, mapping)
    assert result["total_sequences"] == 3 and result["total_residues"] == 7
    assert [row["reference_species"] for row in result["files"]] == ["A", "B"]
    assert "MUZ" in (tmp_path / "b.fasta").read_text()


@pytest.mark.parametrize("text", ["", ">sp|C|C_B\n", ">bad\nMA\n",
    ">sp|X|X_B\nMA\n", ">sp|A|A_A\nMA\n", ">sp|C|C_B\nMA\n>sp|C|C_B\nMA\n"])
def test_bad_fasta(tmp_path, text):
    mapping = setup(tmp_path)
    (tmp_path / "b.fasta").write_text(text)
    with pytest.raises(ValueError):
        run(tmp_path, mapping)


@pytest.mark.parametrize("offsets", [[0, 0, 3], [1, 2, 3], [0, 2], [0, True, 3]])
def test_invalid_offsets(tmp_path, offsets):
    mapping = setup(tmp_path)
    mapping["Goff"] = offsets
    with pytest.raises(ValueError):
        run(tmp_path, mapping)


def test_cross_species_mixture(tmp_path):
    mapping = setup(tmp_path)
    mapping["mapping"].update(B=3, C=2)
    with pytest.raises(ValueError, match="mixed reference species"):
        run(tmp_path, mapping)


def test_missing_numeric_id(tmp_path):
    mapping = setup(tmp_path)
    (tmp_path / "a.fasta").write_text(">sp|A|A_A\nMA\n")
    with pytest.raises(ValueError, match="incomplete proteome"):
        run(tmp_path, mapping)


def test_duplicate_numeric_alias(tmp_path):
    mapping = setup(tmp_path)
    mapping["mapping"]["B"] = 1
    with pytest.raises(ValueError, match="invalid staged identity"):
        run(tmp_path, mapping)


def test_extra_file(tmp_path):
    mapping = setup(tmp_path)
    (tmp_path / "extra.fa").write_text(">extra\nMA\n")
    with pytest.raises(ValueError, match="directory entries"):
        run(tmp_path, mapping)


def test_changed_checksum(tmp_path):
    mapping = setup(tmp_path)
    inputs = [record(tmp_path / name) for name in ("a.fasta", "b.fasta")]
    (tmp_path / "b.fasta").write_text(">sp|C|C_B\nMA\n")
    with pytest.raises(ValueError):
        inventory(tmp_path, inputs, mapping)


def test_symlink(tmp_path):
    mapping = setup(tmp_path)
    inputs = [record(tmp_path / name) for name in ("a.fasta", "b.fasta")]
    target = tmp_path / "b.fasta"
    target.unlink()
    target.symlink_to(tmp_path / "a.fasta")
    with pytest.raises(ValueError, match="regular file"):
        inventory(tmp_path, inputs, mapping)
