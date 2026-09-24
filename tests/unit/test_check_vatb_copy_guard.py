from collections import Counter
from pathlib import Path

import pytest

from benchmark_tools.check_vatb_copy_guard import frozen_predicate, singleton_members


@pytest.mark.parametrize("size,copies,species,expected", [(150, 10, 50, True), (149, 10, 50, False),
    (150, 9, 50, False), (150, 10, 49, False), (276, 10, 78, True)])
def test_frozen_predicate_boundaries(size, copies, species, expected):
    source = Path(__file__).resolve().parents[2] / "orthohmm/refinement.py"
    defaults, predicate = frozen_predicate(source.read_text())
    assert defaults["DEFAULT_COPY_SPLIT_MIN_SIZE"] == 150
    assert predicate(size, Counter({0: copies}), species, 150, 10, 50) is expected


@pytest.mark.parametrize("text,result", [("a\nb\nc d\n", 2), ("a b c\n", 0), ("a\na b c\n", None), ("a\nb\n", None)])
def test_full_group_singleton_check(tmp_path, text, result):
    path = tmp_path / "partition"
    path.write_text(text)
    if result is None:
        with pytest.raises(ValueError):
            singleton_members(path, {"a", "b", "c"})
    else:
        assert singleton_members(path, {"a", "b", "c"}) == result
