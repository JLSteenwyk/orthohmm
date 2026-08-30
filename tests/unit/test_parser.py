import pytest

from orthohmm.parser import create_parser


@pytest.fixture
def parser():
    return create_parser()


class TestParser(object):
    def test_required_only(self, parser):
        fasta_directory = "./tests/samples/"
        parsed = parser.parse_args([fasta_directory])
        assert parsed.fasta_directory == fasta_directory
        assert parsed.clustering == "leiden"
        assert parsed.cpm_resolution == "0.1"
        assert parsed.threads_per_worker == 8
        assert parsed.phylogeny == "off"
        assert parsed.species_tree_mode == "supplied"
        assert parsed.species_tree is None
        assert parsed.aligner == "mafft"
        assert parsed.tree_builder == "FastTree"
