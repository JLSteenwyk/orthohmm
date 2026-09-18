import hashlib
import io
import tarfile

import pytest

from benchmark_tools.compare_qfo_corrected_archive import compare


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


def test_changes_mapping_and_reference_recovery(tmp_path):
    old = b">sp|A|ONE\nACD\n"
    new = b">tr|B|TWO\nAAA\n>tr|C|THREE\nCCC\n"
    path = make_archive(tmp_path, [("a.fasta", old), ("b.fasta", new), ("a_additional.fasta", b"ignored")])
    result = compare(path, {"a.fasta": identity(old), "b.fasta": identity(old)}, {"A": 1, "B": 2}, {"B", "D"})
    assert result["changed_canonical_files"] == ["b.fasta"]
    assert result["sequences"] == 3
    assert result["unique_mapped_numeric_ids"] == 2
    assert result["unmapped_accessions"] == 1
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
