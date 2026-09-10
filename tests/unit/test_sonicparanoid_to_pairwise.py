from io import StringIO

import pytest

from benchmark_tools.sonicparanoid_to_pairwise import iter_pairs, write_pairs


def test_converts_native_sonicparanoid_relations(tmp_path):
    species = tmp_path / "species_a"
    species.mkdir()
    table = species / "species_a-species_b"
    table.write_text(
        "Size\tRelations\tOrthoA\tOrthoB\n"
        "3\t2\tsp|B|GENE_B 1\ttr|A|GENE_A 1 tr|C|GENE_C 0.5\n"
    )
    output = StringIO()

    assert write_pairs(tmp_path, output) == (1, 2, 0)
    assert output.getvalue() == "A\tB\nB\tC\n"


def test_rejects_incorrect_relation_count(tmp_path):
    table = tmp_path / "pairs"
    table.write_text(
        "Size\tRelations\tOrthoA\tOrthoB\n"
        "3\t1\ta 1\tb 1 c 0.5\n"
    )

    with pytest.raises(ValueError, match="relation mismatch"):
        list(iter_pairs(table))


def test_rejects_duplicate_species_pair_tables(tmp_path):
    for left, right in (("species_a", "species_b"), ("species_b", "species_a")):
        directory = tmp_path / left
        directory.mkdir()
        (directory / f"{left}-{right}").write_text(
            "Size\tRelations\tOrthoA\tOrthoB\n"
            "2\t1\ta 1\tb 1\n"
        )

    with pytest.raises(ValueError, match="Duplicate SonicParanoid species-pair"):
        write_pairs(tmp_path, StringIO())


def test_rejects_incomplete_species_pair_matrix(tmp_path):
    for left, right in (("species_a", "species_b"), ("species_a", "species_c")):
        directory = tmp_path / left
        directory.mkdir(exist_ok=True)
        (directory / f"{left}-{right}").write_text(
            "Size\tRelations\tOrthoA\tOrthoB\n"
            "2\t1\ta 1\tb 1\n"
        )

    with pytest.raises(ValueError, match="Incomplete SonicParanoid"):
        write_pairs(tmp_path, StringIO())


def test_deduplicates_relations_within_species_pair(tmp_path):
    table = tmp_path / "pairs"
    table.write_text(
        "Size\tRelations\tOrthoA\tOrthoB\n"
        "2\t1\ta 1\tb 1\n"
        "2\t1\ta 0.5\tb 0.5\n"
    )
    stats = {"duplicates": 0}

    assert list(iter_pairs(table, stats)) == [("a", "b")]
    assert stats == {"duplicates": 1}
