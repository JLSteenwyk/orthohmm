import hashlib
import io
import tarfile

import pytest

from benchmark_tools.compare_qfo_corrected_archive import compare, missing_by_species


def make_archive(tmp_path, entries):
    path = tmp_path / "corrected.tar.gz"
    with tarfile.open(path, "w:gz") as stream:
        for name, data in entries:
            member = tarfile.TarInfo("Eukaryota/" + name)
            member.size = len(data)
            stream.addfile(member, io.BytesIO(data))
    return path


def identity(data):
    return {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def test_species_boundaries_match_native_one_based_rule():
    assert missing_by_species([1, 2, 3, 5], ["A", "B"], [0, 2, 5]) == {"A": 2, "B": 2}


@pytest.mark.parametrize("numbers,species,offsets", [([0], ["A"], [0, 2]), ([3], ["A"], [0, 2]),
    ([1], ["A", "A"], [0, 1, 2]), ([1], ["A"], [1, 2]), ([1], ["A"], [0, 0])])
def test_invalid_species_intervals(numbers, species, offsets):
    with pytest.raises(ValueError):
        missing_by_species(numbers, species, offsets)


def test_changes_mapping_and_reference_recovery(tmp_path):
    old = b">sp|A|ONE\nACD\n"
    new = b">tr|B|TWO\nAAA\n>tr|C|THREE\nCCC\n"
    path = make_archive(tmp_path, [("a.fasta", old), ("b.fasta", new), ("a_additional.fasta", b"ignored")])
    result = compare(path, {"a.fasta": identity(old), "b.fasta": identity(old)}, {"A": 1, "B": 2, "D": 4}, {"B", "D"})
    assert result["changed_canonical_files"] == ["b.fasta"]
    assert result["sequences"] == 3
    assert result["unique_mapped_numeric_ids"] == 2
    assert result["unmapped_accessions"] == 1
    assert result["mapping_numeric_identity_count"] == 3
    assert result["mapping_numeric_ids_without_canonical_accession"] == [4]
    assert result["remaining_missing_accessions"] == ["D"]
    assert result["recovered_missing_reference_accessions"]["B"]["numeric_protein_id"] == 2


@pytest.mark.parametrize("kind", ["duplicate_file", "extra_file", "missing_file", "duplicate_accession", "numeric_alias", "bad_header"])
def test_reject_ambiguous_or_incomplete_inputs(tmp_path, kind):
    data = b">sp|A|ONE\nACD\n"
    expected = {"a.fasta": identity(data)}
    entries, mapping = [("a.fasta", data)], {"A": 1, "B": 1}
    if kind == "duplicate_file":
        entries *= 2
    elif kind == "extra_file":
        entries.append(("b.fasta", data))
    elif kind == "missing_file":
        expected["b.fasta"] = identity(data)
    elif kind == "duplicate_accession":
        entries = [("a.fasta", data + data)]
    elif kind == "numeric_alias":
        entries = [("a.fasta", data + b">sp|B|TWO\nAAA\n")]
    else:
        entries = [("a.fasta", b">A\nAAA\n")]
    with pytest.raises(ValueError):
        compare(make_archive(tmp_path, entries), expected, mapping, set())
