from copy import deepcopy
import hashlib
import io
import tarfile

import pytest

from benchmark_tools.stage_qfo_corrected_inputs import (
    ALIASES_SHA, DATABASE_SHA, MAPPING_SHA, PREPARED_SHA, REFERENCE_COUNT,
    XENOPUS, extract_canonical, validate_reports,
)
from benchmark_tools.audit_qfo_archive_sequences import classify_matches


def reports():
    files = {f"species_{i}.fasta": {"identical_to_original": True, "sequences": 1,
             "mapped_accessions": 1, "unmapped_accessions": []} for i in range(77)}
    files[XENOPUS] = {"identical_to_original": False, "sequences": REFERENCE_COUNT - 77,
                      "mapped_accessions": REFERENCE_COUNT - 77, "unmapped_accessions": []}
    shared = [{"path": "archive", "bytes": 2648666198, "sha256": "archive-sha"},
              {"path": "prepared", "sha256": PREPARED_SHA}, {"path": "mapping", "sha256": MAPPING_SHA}]
    comparison = {"status": "qfo_canonical_archive_compared", "inputs": shared + [{"sha256": ALIASES_SHA}],
                  "canonical_files": files, "canonical_proteomes": 78, "changed_canonical_files": [XENOPUS],
                  "sequences": REFERENCE_COUNT, "unique_mapped_numeric_ids": REFERENCE_COUNT,
                  "unmapped_accessions": 0, "mapping_numeric_identity_count": REFERENCE_COUNT,
                  "mapping_numeric_ids_without_canonical_accession": [],
                  "missing_numeric_ids_by_species": {f"species_{i}": 0 for i in range(78)},
                  "gzip_read_to_eof": True, "remaining_missing_accessions": [],
                  "recovered_missing_reference_accessions": {str(i): {"numeric_protein_id": i + 1} for i in range(14)}}
    native = {"native_entries": REFERENCE_COUNT, "mapped_input_sequences": REFERENCE_COUNT,
              "sequence_identical": REFERENCE_COUNT, "sequence_different": 0, "differences": [],
              "native_entries_without_input_by_species": {}}
    sequences = {"status": "archive_native_sequence_comparison", "inputs": shared + [{"sha256": DATABASE_SHA}],
                 "archive_comparison": deepcopy(comparison), "inputs_modified": False,
                 "complete_mapping_coverage": True, "native_sequence_comparison": native,
                 "sequence_match_classes": classify_matches(native)}
    return comparison, sequences, set(files), set(comparison["recovered_missing_reference_accessions"])


def test_compatible_reports():
    c, s, names, missing = reports()
    assert validate_reports(c, s, names, missing) == c["canonical_files"]


@pytest.mark.parametrize("value", [1, -1, False, "0", None])
def test_invalid_per_species_missing_count(value):
    c, s, names, missing = reports()
    c["missing_numeric_ids_by_species"]["species_0"] = value
    s["archive_comparison"]["missing_numeric_ids_by_species"]["species_0"] = value
    with pytest.raises(ValueError, match="per-species missing"):
        validate_reports(c, s, names, missing)


def test_empty_missing_count_dictionary_is_not_full_inventory():
    c, s, names, missing = reports()
    c["missing_numeric_ids_by_species"] = {}
    s["archive_comparison"]["missing_numeric_ids_by_species"] = {}
    with pytest.raises(ValueError, match="per-species missing"):
        validate_reports(c, s, names, missing)


def test_representation_match_is_not_exact():
    c, s, names, missing = reports()
    native = s["native_sequence_comparison"]
    native.update(sequence_identical=REFERENCE_COUNT - 1, sequence_different=1,
                  differences=[{"numeric_id": 1, "BOUZ_to_X_sequence": "AX", "native_sequence": "AX"}])
    s["sequence_match_classes"] = classify_matches(native)
    validate_reports(c, s, names, missing)
    assert s["sequence_match_classes"]["exact_sequence_matches"] == REFERENCE_COUNT - 1
    native["differences"][0]["native_sequence"] = "AA"
    s["sequence_match_classes"] = classify_matches(native)
    with pytest.raises(ValueError, match="Unexplained"):
        validate_reports(c, s, names, missing)


@pytest.mark.parametrize("key,value", [("canonical_proteomes", 77), ("sequences", REFERENCE_COUNT - 1),
                                      ("unmapped_accessions", 1), ("gzip_read_to_eof", False),
                                      ("changed_canonical_files", [XENOPUS, "species_0.fasta"])])
def test_incompatible_archive_reports(key, value):
    c, s, names, missing = reports()
    c[key] = value
    s["archive_comparison"][key] = value
    with pytest.raises(ValueError):
        validate_reports(c, s, names, missing)


def test_cross_report_disagreement():
    c, s, names, missing = reports()
    s["archive_comparison"]["canonical_files"][XENOPUS]["sequences"] -= 1
    with pytest.raises(ValueError, match="disagree"):
        validate_reports(c, s, names, missing)


def test_missing_recovered_accession():
    c, s, names, missing = reports()
    del c["recovered_missing_reference_accessions"]["0"]
    with pytest.raises(ValueError, match="recovery"):
        validate_reports(c, s, names, missing)


def test_native_missing_entries():
    c, s, names, missing = reports()
    s["native_sequence_comparison"]["native_entries_without_input_by_species"] = {"XENTR": 1}
    with pytest.raises(ValueError, match="native database"):
        validate_reports(c, s, names, missing)


def test_changed_source_identity():
    c, s, names, missing = reports()
    s["inputs"] = deepcopy(s["inputs"])
    s["inputs"][0]["sha256"] = "different"
    with pytest.raises(ValueError, match="source identities"):
        validate_reports(c, s, names, missing)


def archive_fixture(tmp_path, member_name="Eukaryota/test.fasta", duplicate=False, symlink=False):
    data = b">sp|A|ONE description\nAUZ\n"
    archive = tmp_path / "archive.tar.gz"
    with tarfile.open(archive, "w:gz") as stream:
        for _ in range(2 if duplicate else 1):
            member = tarfile.TarInfo(member_name)
            member.size = len(data)
            if symlink:
                member.type = tarfile.SYMTYPE
                member.linkname = "/tmp/unwanted"
                member.size = 0
                stream.addfile(member)
            else:
                stream.addfile(member, io.BytesIO(data))
        ignored = tarfile.TarInfo("Eukaryota/test_DNA.fasta")
        ignored.size = 3
        stream.addfile(ignored, io.BytesIO(b"AAA"))
    files = {"test.fasta": {"archive_member": member_name, "bytes": len(data),
                             "sha256": hashlib.sha256(data).hexdigest()}}
    return archive, files, tmp_path / "stage", data


def test_flat_exclusive_unmodified_extraction(tmp_path):
    archive, files, destination, data = archive_fixture(tmp_path)
    result = extract_canonical(archive, files, destination)
    assert len(result) == 1
    assert (destination / "test.fasta").read_bytes() == data
    assert not (destination / "test_DNA.fasta").exists()
    with pytest.raises(FileExistsError):
        extract_canonical(archive, files, destination)


@pytest.mark.parametrize("options", [{"duplicate": True}, {"symlink": True},
                                     {"member_name": "../test.fasta"}, {"member_name": "/test.fasta"}])
def test_unsafe_archive(tmp_path, options):
    archive, files, destination, _ = archive_fixture(tmp_path, **options)
    with pytest.raises(ValueError):
        extract_canonical(archive, files, destination)
    assert not (destination / "staging_manifest.json").exists()


def test_checksum_mismatch_preserves_partial_directory(tmp_path):
    archive, files, destination, _ = archive_fixture(tmp_path)
    files["test.fasta"]["sha256"] = "wrong"
    with pytest.raises(ValueError, match="checksum"):
        extract_canonical(archive, files, destination)
    assert destination.exists() and not (destination / "staging_manifest.json").exists()


def test_missing_archive_member(tmp_path):
    archive, files, destination, _ = archive_fixture(tmp_path)
    files["missing.fasta"] = files["test.fasta"]
    with pytest.raises(ValueError, match="Missing canonical"):
        extract_canonical(archive, files, destination)


def test_unsafe_output_name(tmp_path):
    archive, files, destination, _ = archive_fixture(tmp_path)
    files["../test.fasta"] = files.pop("test.fasta")
    with pytest.raises(ValueError, match="Unsafe output"):
        extract_canonical(archive, files, destination)
    assert not destination.exists()


def test_truncated_gzip(tmp_path):
    archive, files, destination, _ = archive_fixture(tmp_path)
    archive.write_bytes(archive.read_bytes()[:-4])
    with pytest.raises(EOFError):
        extract_canonical(archive, files, destination)
    assert not (destination / "staging_manifest.json").exists()
