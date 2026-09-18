import math

import pytest

from benchmark_tools.inventory_swiss_sequences import collect, describe, summarize
from benchmark_tools.snapshot_orthohmm_input_order import record


def test_entropy_and_noncanonical_denominator():
    value = describe("AACDXX", "protein")
    assert value["canonical_entropy_bits"] == 1.5
    assert value["maximum_canonical_frequency"] == .5
    assert value["noncanonical_fraction"] == pytest.approx(1 / 3)
    assert value["noncanonical_counts"] == {"X": 2}
    assert describe("ACDEFGHIKLMNPQRSTVWY", "")["canonical_entropy_bits"] == pytest.approx(math.log2(20))


def test_unknown_only_not_zero_entropy():
    assert describe("XX", "")["canonical_entropy_bits"] is None
    assert describe("AAAA", "")["canonical_entropy_bits"] == 0


@pytest.mark.parametrize("sequence", ["", "AA A", "AA\nA"])
def test_invalid_sequence(sequence):
    with pytest.raises(ValueError):
        describe(sequence, "")


def test_fragment_literal_not_inference():
    assert describe("AA", "name (Fragments)")["explicit_fragment_description"]
    assert describe("AA", "name (fragment)")["explicit_fragment_description"]
    assert not describe("AA", "fragment-binding protein")["explicit_fragment_description"]


def test_collect_and_missing(tmp_path):
    path = tmp_path / "input.fa"
    path.write_text(">sp|P1|TEST name (Fragment)\nAACD\n>tr|P2|TEST name\nXXXX\n")
    result = collect({"A": ["P1", "P3"], "B": ["P2"]}, [record(path)])
    assert result["summary"]["matched_genes"] == 2
    assert result["summary"]["missing_genes"] == ["P3"]
    assert result["families"]["B"]["median_canonical_entropy_bits"] is None
    assert result["families"]["A"]["explicit_fragment_descriptions"] == 1
    assert summarize(["absent"], {})["median_length"] is None


@pytest.mark.parametrize("mode", ["shared", "changed", "duplicate_path", "duplicate_accession", "header"])
def test_reject_bad_sources(tmp_path, mode):
    path = tmp_path / "input.fa"
    path.write_text(">sp|P1|TEST name\nACD\n")
    families, inputs = {"A": ["P1"]}, [record(path)]
    if mode == "shared":
        families["B"] = ["P1"]
    elif mode == "changed":
        path.write_text(">sp|P1|TEST\nCCC\n")
    elif mode == "duplicate_path":
        inputs *= 2
    else:
        path.write_text(">P1\nCCC\n" if mode == "header" else ">sp|P1|ONE\nACD\n>tr|P1|TWO\nACD\n")
        inputs = [record(path)]
    with pytest.raises(ValueError):
        collect(families, inputs)
