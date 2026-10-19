import json
from pathlib import Path

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import pytest

from orthohmm.search import selective_hmmer as hmmer
from orthohmm.search.matrices import get_background_freqs, get_matrix
from orthohmm.search.msa_profile import build_msa_profile
from tests.unit.test_selective_hmmer import tools


def write_alignment(path, sequences, rf):
    alignment = MultipleSeqAlignment([SeqRecord(Seq(sequence), id=name, description="") for name, sequence in sequences])
    alignment.column_annotations["reference_annotation"] = rf
    AlignIO.write(alignment, path, "stockholm")


def fixture(tmp_path):
    source, output = tmp_path / "source.sto", tmp_path / "annotated.sto"
    rows = [("g0", "--ACXDEFGH--"), ("g1", "--ACXDEFGH--"), ("g2", "---CXDEFG---")]
    annotated = [("g0", "~~ACXDEFGH~~"), ("g1", "..ACXDEFGH.."), ("g2", "~~~CXDEFG~~~")]
    rf = "..xxxxxxxx.."
    write_alignment(source, rows, rf)
    write_alignment(output, annotated, rf)
    return source, output, {"g0": "a", "g1": "b", "g2": "c", "query": "candidate"}, ["a", "b", "c"], 8


def test_fragment_missing_data_symbols_are_not_source_residues(tmp_path):
    args = fixture(tmp_path)
    result = hmmer.validate_annotated_alignment(*args)
    assert result["members"] == 3 and result["match_columns"] == 8
    assert result["missing_data_symbols"] == 10 and result["all_source_residues_and_owners_preserved"]


@pytest.mark.parametrize("mutation", ["substitution", "lost_X", "lost_residue", "wrong_owner", "duplicate_owner", "missing_owner", "unknown_symbol", "rf", "missing_rf", "expected_native"])
def test_genuine_source_or_match_mask_changes_still_fail(tmp_path, mutation):
    args = fixture(tmp_path)
    annotated = AlignIO.read(args[1], "stockholm")
    if mutation == "substitution": annotated[0].seq = Seq(str(annotated[0].seq).replace("A", "Y"))
    elif mutation == "lost_X": annotated[0].seq = Seq(str(annotated[0].seq).replace("X", "~"))
    elif mutation == "lost_residue": annotated[0].seq = Seq(str(annotated[0].seq).replace("A", "~"))
    elif mutation == "wrong_owner": annotated[0].id = "foreign"
    elif mutation == "duplicate_owner":
        args[2]["g2"] = "b"
        with pytest.raises(ValueError): hmmer.validate_annotated_alignment(*args)
        return
    elif mutation == "missing_owner": annotated = annotated[1:]
    elif mutation == "unknown_symbol": annotated[0].seq = Seq(str(annotated[0].seq).replace("A", "?"))
    elif mutation == "rf": annotated.column_annotations["reference_annotation"] = "...xxxxxxx.."
    elif mutation == "missing_rf": annotated.column_annotations.clear()
    else:
        with pytest.raises(ValueError): hmmer.validate_annotated_alignment(*args, expected_residues={"g0": "changed_native"})
        return
    AlignIO.write(annotated, args[1], "stockholm")
    with pytest.raises(ValueError): hmmer.validate_annotated_alignment(*args)


def test_empty_insertion_columns_can_be_removed_without_losing_source_residues(tmp_path):
    args = list(fixture(tmp_path))
    output = AlignIO.read(args[1], "stockholm")
    output = output[:, 2:-2]
    AlignIO.write(output, args[1], "stockholm")
    result = hmmer.validate_annotated_alignment(*args)
    assert result["source_columns"] == 12 and result["annotated_columns"] == 8
    assert result["match_columns"] == 8


def test_real_hmmer_fragment_source_is_scored_without_changing_build_controls(tmp_path, tools):
    rng = np.random.default_rng(20261019)
    sequence = "".join(rng.choice(list("ACDEFGHIKLMNPQRSTVWY"), 200))
    members = {"long_a": sequence, "long_b": sequence, "fragment": sequence[70:110]}
    profile = build_msa_profile(list(members.values()), list(members), get_matrix("BLOSUM62"),
        get_background_freqs("BLOSUM62"), projection_policy="residue_preserving")
    result = hmmer.score_family_hmmer(profile, members, {"query": sequence}, 10000, tmp_path / "actual", **tools)
    assert result["status"] == "hmmer_profile_scored"
    assert result["annotated_alignment_validation"]["missing_data_symbols"] > 0
    assert "--fragthresh" not in result["commands"]["build"]
    assert result["parameters"]["effective_sequence_number"] == 3
    assert "--nonull2" not in result["commands"]["search"]
    assert json.loads((tmp_path / "actual/measurements.json").read_text()) == result


def test_retained_twenty_profile_failure_is_format_only_not_changed_residues():
    root = Path(__file__).resolve().parents[2] / "benchmarks/results/projection_recovery_pilot_v1_20261019/projection_hmmer_score_only/graph/selective_profiles"
    paths = sorted(root.glob("*/measurements.json"))
    if not paths: pytest.skip("Requires retained terminal recovery receipts")
    assert len(paths) == 20, "Retained terminal recovery failure is required for this integration check"
    original_receipts = [path.read_bytes() for path in paths]
    statuses = []
    for path in paths:
        receipt = json.loads(path.read_text())
        statuses.append(receipt["status"])
        result = hmmer.validate_annotated_alignment(path.parent / "source.sto", path.parent / "hmmer_source.sto",
            receipt["token_to_gene"], receipt["members"], receipt["model_length"])
        assert result["all_source_residues_and_owners_preserved"]
        assert result["match_columns"] == receipt["model_length"]
    assert statuses.count("failed") == 19 and statuses.count("hmmer_profile_scored") == 1
    assert [path.read_bytes() for path in paths] == original_receipts
