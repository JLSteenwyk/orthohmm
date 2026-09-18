import pytest

from benchmark_tools.validate_sonicparanoid_tables import validate_tables


HEADER = "Size\tRelations\tOrthoA\tOrthoB\n"
ROW = "2\t1\ta 1\tb 0.5\n"


def make(tmp_path, text=HEADER + ROW):
    inputs = []
    for name, gene in (("s1", "a"), ("s2", "b"), ("s3", "c")):
        p = tmp_path / (name + ".fasta")
        p.write_text(f">{gene}\nACDE\n")
        inputs.append(p)
    directory = tmp_path / "tables"
    for left, right in (("s1", "s2"), ("s1", "s3"), ("s2", "s3")):
        p = directory / (left + ".fasta") / f"{left}.fasta-{right}.fasta"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text if right == "s2" else HEADER)
    return directory, inputs


def test_empty_tables_and_repeated_relations(tmp_path):
    result = validate_tables(*make(tmp_path, HEADER + ROW + ROW))
    assert result == {"input_species": 3, "input_accessions": 3,
                      "species_pair_tables": 3, "native_rows": 2,
                      "raw_relations": 2, "distinct_pairs": 1, "duplicate_relations": 1}


@pytest.mark.parametrize("row", [
    "2\t1\tc 1\tb 1\n", "2\t1\tunknown 1\tb 1\n",
    "2\t1\ta nan\tb 1\n", "2\t1\ta inf\tb 1\n",
    "3\t1\ta 1\tb 1\n", "2\t2\ta 1\tb 1\n",
    "3\t2\tsp|a|x 1 tr|a|y 1\tb 1\n", "bad\n",
])
def test_reject_invalid_members_and_counts(tmp_path, row):
    with pytest.raises(ValueError):
        validate_tables(*make(tmp_path, HEADER + row))


def test_reject_entire_missing_species(tmp_path):
    directory, inputs = make(tmp_path)
    for path in directory.rglob("*s3.fasta"):
        if path.is_file():
            path.unlink()
    with pytest.raises(ValueError, match="Incomplete"):
        validate_tables(directory, inputs)


def test_reject_reversed_duplicate_table(tmp_path):
    directory, inputs = make(tmp_path)
    (directory / "s2.fasta" / "s2.fasta-s1.fasta").write_text(HEADER)
    with pytest.raises(ValueError, match="duplicate"):
        validate_tables(directory, inputs)
