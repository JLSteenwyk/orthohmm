import pytest

from orthohmm.parser import create_parser


@pytest.fixture
def parser():
    return create_parser()


class TestParser(object):
    def test_required_only(self, parser):
        input_directory = "./tests/samples/"
        parsed = parser.parse_args([input_directory])
        assert parsed.input == input_directory
