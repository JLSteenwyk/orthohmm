from pathlib import Path

import pytest

from benchmark_tools.verify_frozen_scoring import _definitions, verify


def test_decorators_and_import_side_effects_are_not_executed():
    env = {}
    _definitions("import nonexistent\n@nonexistent\ndef score(x):\n return x + 1\n",
                 {"score"}, env)
    assert env["score"](2) == 3


def test_missing_definition_rejected():
    with pytest.raises(ValueError, match="absent"):
        _definitions("def other():\n return 0\n", {"score"}, {})


def test_exact_frozen_fixtures_are_bounded():
    result = verify(Path(__file__).resolve().parents[2])
    assert len(result["source"]) == 9
    assert result["fixtures"]["ACDE_self_raw_score"] == 24
    assert result["fixtures"]["half_gap_column_mask"] == [True, False]
    for key in ("native_backend_validation", "calibration_established",
                "biological_validation", "frozen_method_changed", "publication_ready"):
        assert result[key] is False
