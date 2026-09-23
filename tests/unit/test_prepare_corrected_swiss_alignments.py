import pytest

from benchmark_tools.prepare_corrected_swiss_alignments import extract_sequences, identity_bins
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def row(name, value):
    return dict(refog=name, status="alignment_validated", identity=dict(mean_pairwise_identity=value))


def test_median_split_ties_and_missing():
    result = identity_bins([row("a", .1), row("b", .5), row("c", .5), row("d", .9), row("e", None)])
    assert result == dict(median_family_identity=.5,
        strata=dict(lower_identity=["a", "b", "c"], higher_identity=["d"], missing_identity=["e"]))


def test_all_missing_not_zero():
    result = identity_bins([row("a", None)])
    assert result["median_family_identity"] is None
    assert result["strata"]["missing_identity"] == ["a"]


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -.1, 1.1, True])
def test_invalid_identity_rejected(value):
    with pytest.raises(ValueError, match="identity"):
        identity_bins([row("a", value)])


@pytest.mark.parametrize("rows", [[], [row("a", .5), row("a", .5)],
                                  [dict(row("a", .5), status="failed")]])
def test_partial_or_duplicate_panel_rejected(rows):
    with pytest.raises(ValueError, match="successful"):
        identity_bins(rows)


@pytest.mark.parametrize("problem", [None, "duplicate", "missing", "length", "alphabet", "hash"])
def test_reference_extraction(tmp_path, problem):
    fasta = tmp_path / "input.faa"
    fasta.write_text(">sp|A|name\nACDX\n>sp|B|name\nACD*\n>sp|OTHER|name\nAAA\n")
    if problem == "duplicate":
        fasta.write_text(fasta.read_text() + ">sp|A|again\nACDX\n")
    elif problem == "missing":
        fasta.write_text(">sp|A|name\nACDX\n")
    elif problem == "alphabet":
        fasta.write_text(">sp|A|name\nACDU\n>sp|B|name\nACD*\n")
    inventory = dict(family_memberships={"family": ["A", "B"]}, fasta_inputs=[record(fasta)],
                     genes={"A": {"length": 4}, "B": {"length": 4}})
    if problem == "length":
        inventory["genes"]["A"]["length"] = 5
    elif problem == "hash":
        fasta.write_text(fasta.read_text() + "\n")
    if problem:
        with pytest.raises(ValueError):
            extract_sequences(inventory)
    else:
        assert extract_sequences(inventory) == {"A": "ACDX", "B": "ACD*"}
