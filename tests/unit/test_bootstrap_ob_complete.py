from fractions import Fraction
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from benchmark_tools.bootstrap_ob_complete import run
from benchmark_tools.check_ob_complete_uncertainty import (
    check_family_outcomes, family_outcomes, numerical_error, ratios,
    verify_portable,
)


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


@pytest.mark.parametrize("left,right", [
    (1., float("nan")), (float("nan"), 1.),
    (float("inf"), float("inf")), ([1., 2.], [1., float("nan")]),
    ([1., 2.], 1.), ([1., 2.], [[1., 2.]]), (1., 1.01),
])
def test_invalid_numerical_comparisons_rejected(left, right):
    with pytest.raises(ValueError):
        numerical_error(left, right)


def test_numerical_tolerance():
    assert numerical_error([1., 2.], [1., 2. + 1e-12]) < 1e-10


def family(name, tp, fp=0):
    return dict(refog=name, genes=3, true_positive=tp, false_positive=fp, false_negative=3-tp)


def test_family_outcomes_rational_and_order_independent():
    method = [family("win", 3), family("tie", 2), family("loss", 0)]
    baseline = [family("loss", 1), family("tie", 2), family("win", 2)]
    result = family_outcomes(method, baseline)
    assert result == dict(family_f1_wins=1, family_f1_ties=1, family_f1_losses=1)
    check_family_outcomes(result, result)


def test_fractional_counts_and_tie_tolerance():
    baseline = [family("a", 1.5)]
    assert family_outcomes([family("a", 1.5 + 1e-13)], baseline)["family_f1_ties"] == 1
    assert family_outcomes([family("a", 1.5 + 1e-9)], baseline)["family_f1_wins"] == 1


@pytest.mark.parametrize("records", [
    [], [family("a", 1), family("a", 1)], [family("b", 1)],
    [family("a", -1)], [dict(family("a", 1), genes=4)],
    [dict(family("a", 1), genes=True)],
])
def test_family_inventory_and_counts_rejected(records):
    with pytest.raises(ValueError):
        family_outcomes(records, [family("a", 1)])


@pytest.mark.parametrize("value", [0, True, 1., None])
def test_wrong_family_outcome_rejected(value):
    with pytest.raises(ValueError, match="outcome"):
        check_family_outcomes({"family_f1_wins": 1}, {"family_f1_wins": value})


def test_portable_checksum_failure_precedes_parsing(tmp_path):
    source = tmp_path / "data.json"
    source.write_text("invalid json")
    with pytest.raises(ValueError, match="checksum"):
        verify_portable(source, "0" * 64, tmp_path / "result")
    assert not (tmp_path / "result").exists()


def test_portable_unknown_schema(tmp_path):
    source = tmp_path / "data.json"
    source.write_text('{"schema": "unknown"}')
    with pytest.raises(ValueError, match="schema"):
        verify_portable(source, hashlib.sha256(source.read_bytes()).hexdigest(), tmp_path / "result")


def test_portable_no_overwrite(tmp_path):
    with pytest.raises(FileExistsError):
        verify_portable(tmp_path / "missing", "0" * 64, tmp_path)


def test_standalone_complete_panel(tmp_path):
    root = Path(__file__).resolve().parents[2]
    script = tmp_path / "check.py"
    source = tmp_path / "counts.json"
    shutil.copyfile(root / "benchmark_tools/check_ob_complete_uncertainty.py", script)
    shutil.copyfile(root / "benchmark_tools/results/ob_complete_portable_statistics_20260928.json", source)
    payload = json.loads(source.read_text())
    assert len(payload["scores"]) == 8
    for score in payload["scores"].values():
        assert len(score["refog_records"]) == 70
        assert all(set(row) == {"refog", "genes", "true_positive", "false_positive", "false_negative"}
                   for row in score["refog_records"])
    assert "/mnt/" not in source.read_text() and "/home/" not in source.read_text()
    run = subprocess.run([sys.executable, "-I", "-B", str(script), "--portable", str(source),
                          "--sha256", hashlib.sha256(source.read_bytes()).hexdigest(),
                          "--output", str(tmp_path / "result")], cwd=tmp_path,
                         text=True, capture_output=True, timeout=120)
    assert run.returncode == 0, run.stderr
    result = json.loads((tmp_path / "result/crosscheck.json").read_text())
    assert result["endpoints"] == 21 and len(result["family_outcomes"]) == 7
    assert result["maximum_absolute_error_percentage_points"] < 1e-10
    assert result["verification_scope"] == "portable_derived_counts_only_no_raw_input_verification"
