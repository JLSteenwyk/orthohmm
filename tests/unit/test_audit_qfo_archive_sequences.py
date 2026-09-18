import hashlib
import io
import tarfile

import pytest

from benchmark_tools.compare_qfo_corrected_archive import compare
from benchmark_tools.audit_qfo_input_sequences import compare_database, sequence_identity
from benchmark_tools.audit_qfo_archive_sequences import classify_matches


def test_archive_to_native_exact_replacement_and_other(tmp_path):
    fasta = b">sp|A|ONE\nACD\n>sp|B|TWO\nAUO\n>tr|C|THREE\nABC\n"
    archive = tmp_path / "source.tar.gz"
    with tarfile.open(archive, "w:gz") as stream:
        member = tarfile.TarInfo("Eukaryota/test.fasta")
        member.size = len(fasta)
        stream.addfile(member, io.BytesIO(fasta))
    expected = {"test.fasta": {"bytes": len(fasta), "sha256": hashlib.sha256(fasta).hexdigest()}}
    sequences = {}
    report = compare(archive, expected, {"A": 1, "B": 2, "C": 3}, set(), sequence_records=sequences)
    assert report["changed_canonical_files"] == []
    assert sequences[2]["BOUZ_to_X_sequence"] == sequence_identity("AXX")
    database = tmp_path / "db"
    database.write_text("".join(f"<E><OS>S</OS><MAPIDS>{a}</MAPIDS><SEQ>{s}</SEQ></E>\n" for a, s in [("A", "ACD"), ("B", "AXX"), ("C", "AXD"), ("D", "AA")]))
    native = compare_database(database, sequences, 4)
    classes = classify_matches(native)
    assert classes["exact_sequence_matches"] == 1
    assert classes["BOUZ_to_X_only_numeric_ids"] == [2]
    assert classes["unexplained_sequence_difference_ids"] == [3]
    assert not classes["all_mapped_sequences_accounted_for_by_exact_or_BOUZ_to_X"]
    assert native["native_entries_without_input_by_species"] == {"S": 1}
    assert not (tmp_path / "test.fasta").exists()


def test_accounting_cannot_hide_difference():
    with pytest.raises(ValueError):
        classify_matches({"differences": [], "sequence_identical": 1, "mapped_input_sequences": 2})


def test_empty_differences_not_complete_reference_claim():
    classes = classify_matches({"differences": [], "sequence_identical": 1, "mapped_input_sequences": 1})
    assert classes["all_mapped_sequences_accounted_for_by_exact_or_BOUZ_to_X"]
    assert "complete_mapping_coverage" not in classes
