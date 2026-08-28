from io import StringIO

import pytest

from benchmark_tools.orthofinder_to_pairwise import iter_pairs, write_pairs


def _write_pair_table(path, species_a, species_b, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"Orthogroup\t{species_a}\t{species_b}"]
    lines.extend("\t".join(row) for row in rows)
    path.write_text("\n".join(lines) + "\n")


def test_iter_pairs_uses_native_relationships_once(tmp_path):
    results = tmp_path / "Results_Test"
    orthologues = results / "Orthologues"
    _write_pair_table(
        orthologues / "Orthologues_A" / "A__v__B.tsv",
        "A",
        "B",
        [("OG1", "sp|A1|A_ONE, tr|A2|A_TWO", "sp|B1|B_ONE")],
    )
    _write_pair_table(
        orthologues / "Orthologues_B" / "B__v__A.tsv",
        "B",
        "A",
        [("OG1", "sp|B1|B_ONE", "sp|A1|A_ONE, tr|A2|A_TWO")],
    )

    assert list(iter_pairs(results)) == [("A1", "B1"), ("A2", "B1")]

    output = StringIO()
    assert write_pairs(results, output) == 2
    assert output.getvalue() == "A1\tB1\nA2\tB1\n"


def test_iter_pairs_rejects_missing_tables(tmp_path):
    with pytest.raises(FileNotFoundError, match="Orthologues directory"):
        list(iter_pairs(tmp_path / "Results_Test"))


def test_iter_pairs_rejects_unexpected_header(tmp_path):
    results = tmp_path / "Results_Test"
    table = results / "Orthologues" / "Orthologues_A" / "A__v__B.tsv"
    table.parent.mkdir(parents=True)
    table.write_text("bad\theader\n")

    with pytest.raises(ValueError, match="Unexpected OrthoFinder table header"):
        list(iter_pairs(results))
