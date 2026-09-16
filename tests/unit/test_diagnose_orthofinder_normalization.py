import numpy as np
import pytest

from benchmark_tools.diagnose_orthofinder_normalization import design_summary, species_lengths


def test_lengths_follow_native_indices_not_file_order(tmp_path):
    path = tmp_path / "Species2.fa"
    path.write_text(">2_1\nACDE\n>2_0\nAAA\n")
    assert species_lengths(path, 2).tolist() == [3, 4]


@pytest.mark.parametrize("text", [">1_0\nAAA\n", ">2_1\nAAA\n", ">2_0\n", ">2_0\nAAA\n>2_0\nCCC\n"])
def test_invalid_native_fasta_rejected(tmp_path, text):
    path = tmp_path / "Species2.fa"
    path.write_text(text)
    with pytest.raises(ValueError):
        species_lengths(path, 2)


def test_rank_deficiency_can_occur_after_selecting_from_varied_lengths():
    assert design_summary([100, 200, 300])["design_rank"] == 2
    assert design_summary([100, 100])["design_rank"] == 1
    assert design_summary([])["design_rank"] == 0
    assert design_summary([100, 100])["observations"] == 2


@pytest.mark.parametrize("products", [[0], [-1], [np.nan], [np.inf]])
def test_invalid_products_rejected(products):
    with pytest.raises(ValueError):
        design_summary(products)
