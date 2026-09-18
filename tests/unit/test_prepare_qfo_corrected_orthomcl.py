import pytest

from benchmark_tools.prepare_qfo_corrected_orthomcl import commands, require_parity


def test_legacy_search_settings_preserved():
    actual = commands("/fresh", "/blastall", "/formatdb")
    assert actual["formatdb"] == ["/formatdb", "-i", "/fresh/all.fa", "-p", "t"]
    assert actual["blast"] == ["/blastall", "-p", "blastp", "-i", "/fresh/all.fa", "-d", "/fresh/all.fa",
                               "-e", "1e-5", "-o", "/fresh/all.blast.partial", "-m", "8", "-a", "180",
                               "-v", "1000", "-b", "1000"]


@pytest.mark.parametrize("threads", [0, -1, True, 1.5])
def test_invalid_threads(threads):
    with pytest.raises(ValueError):
        commands("/fresh", "/blastall", "/formatdb", threads)


@pytest.mark.parametrize("key,value", [("proteomes", 77), ("total_sequences", 976504),
                                      ("identical_sequences", 984136), ("sequence_difference_count", 1),
                                      ("differences", ["changed"]), ("genome_map_complete_and_correct", False)])
def test_requires_exact_corrected_parity(key, value):
    data = {"proteomes": 78, "total_sequences": 984137, "identical_sequences": 984137,
            "sequence_difference_count": 0, "differences": [], "genome_map_complete_and_correct": True}
    require_parity(data)
    data[key] = value
    with pytest.raises(ValueError, match="parity failed"):
        require_parity(data)
