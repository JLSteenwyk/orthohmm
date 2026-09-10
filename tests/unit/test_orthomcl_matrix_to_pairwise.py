from io import StringIO

import pytest

from benchmark_tools.orthomcl_matrix_to_pairwise import write_pairs


def _write_fixture(tmp_path, matrix):
    gg = tmp_path / "all.gg"
    gg.write_text("A: sp|A1|ONE sp|A2|TWO\nB: sp|B1|THREE\n")
    index = tmp_path / "all_ortho.idx"
    index.write_text("0\tsp|A1|ONE\n1\tsp|A2|TWO\n2\tsp|B1|THREE\n")
    matrix_path = tmp_path / "all_ortho.mtx"
    matrix_path.write_text(matrix)
    return matrix_path, index, gg


def test_write_pairs_uses_cross_species_edges_and_validates_symmetry(tmp_path):
    matrix, index, gg = _write_fixture(
        tmp_path,
        """(mclheader
mcltype matrix
dimensions 3x3
)

(mclmatrix
begin

0    1:0.400 2:1.500 $
1    0:0.400 2:0.700 $
2    0:1.500 1:0.700 $
)
""",
    )
    output = StringIO()

    stats = write_pairs(matrix, index, gg, output)

    assert output.getvalue() == "A1\tB1\nA2\tB1\n"
    assert stats == {
        "matrix_genes": 3,
        "directed_edges": 6,
        "same_species_directed_edges": 2,
        "cross_species_pairs": 2,
    }


def test_write_pairs_rejects_nonreciprocal_edges(tmp_path):
    matrix, index, gg = _write_fixture(
        tmp_path,
        """(mclheader
dimensions 3x3
)
(mclmatrix
begin
0 2:1.500 $
1 $
2 $
)
""",
    )

    with pytest.raises(ValueError, match="not reciprocal"):
        write_pairs(matrix, index, gg, StringIO())


def test_write_pairs_rejects_score_asymmetry(tmp_path):
    matrix, index, gg = _write_fixture(
        tmp_path,
        """(mclheader
dimensions 3x3
)
(mclmatrix
begin
0 2:1.500 $
1 $
2 0:1.400 $
)
""",
    )

    with pytest.raises(ValueError, match="not reciprocal"):
        write_pairs(matrix, index, gg, StringIO())


def test_write_pairs_rejects_duplicate_targets(tmp_path):
    matrix, index, gg = _write_fixture(
        tmp_path,
        """(mclheader
dimensions 3x3
)
(mclmatrix
begin
0 2:1.500 2:1.500 $
1 $
2 0:1.500 $
)
""",
    )

    with pytest.raises(ValueError, match="Duplicate OrthoMCL edge"):
        write_pairs(matrix, index, gg, StringIO())


def test_write_pairs_rejects_index_gene_missing_from_gg(tmp_path):
    matrix, index, gg = _write_fixture(tmp_path, "")
    index.write_text("0\tmissing\n")

    with pytest.raises(ValueError, match="absent from GG"):
        write_pairs(matrix, index, gg, StringIO())
