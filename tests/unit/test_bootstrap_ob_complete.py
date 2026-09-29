from fractions import Fraction

import pytest

from benchmark_tools.bootstrap_ob_complete import run
from benchmark_tools.check_ob_complete_uncertainty import ratios


def test_independent_ratio_arithmetic():
    result = ratios([Fraction(2), Fraction(1), Fraction(1, 2)])
    assert result == dict(f_score=Fraction(800, 11), precision=Fraction(200, 3), recall=Fraction(80))


def test_no_overwrite(tmp_path):
    output = tmp_path / "output.json"
    output.write_text("retain")
    with pytest.raises(FileExistsError):
        run(tmp_path, output, "unused")
    assert output.read_text() == "retain"


def test_protocol_mismatch_prevents_analysis(tmp_path):
    (tmp_path / "OB_COMPLETE_UNCERTAINTY_PROTOCOL_20260928.md").write_text("changed")
    with pytest.raises(ValueError, match="Protocol"):
        run(tmp_path, tmp_path / "output.json", "0" * 64)
    assert not (tmp_path / "output.json").exists()
