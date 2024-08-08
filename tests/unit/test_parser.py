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
