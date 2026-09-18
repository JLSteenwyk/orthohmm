import gzip
import hashlib

import pytest

from benchmark_tools.audit_three_kingdoms_sources import CODES, change_details, complete_hits, decompressed_identity, sequence_inventory, source_rows


def script():
    return "SPECIES=(\n" + "\n".join(f'"species_{code}|{code}|${{EBI}}/archive/{code}.fa.gz"' for code in sorted(CODES)) + "\n)"


def test_source_parser_lexes_without_executing():
    parsed = source_rows(script() + "\nexit 1\n")
    assert set(parsed) == CODES
    assert parsed["Arab"]["source_url"] == "https://ftp.ebi.ac.uk/archive/Arab.fa.gz"
    assert not parsed["Arab"]["moving_release_url"]


@pytest.mark.parametrize("mutation", ["missing", "duplicate", "unresolved"])
def test_invalid_source_inventory(mutation):
    text = script()
    if mutation == "missing":
        text = text.replace("|Arab|", "|Other|")
    elif mutation == "duplicate":
        text += '\n"species_Arab|Arab|https://example.invalid/a"'
    else:
        text = text.replace("${EBI}", "${UNKNOWN}")
    with pytest.raises(ValueError):
        source_rows(text)


def test_sequence_and_compressed_identity(tmp_path):
    path = tmp_path / "genes.fa"
    content = b">a description\nAbC\n>b\nDEF\n"
    path.write_bytes(content)
    genes = sequence_inventory(path)
    assert genes["a"] == {"length": 3, "sha256": hashlib.sha256(b"AbC").hexdigest()}
    compressed = tmp_path / "genes.gz"
    compressed.write_bytes(gzip.compress(content))
    assert decompressed_identity(compressed) == {"bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()}


@pytest.mark.parametrize("content", ["", ">a\n", ">a\nABC\n>a\nDEF\n"])
def test_bad_fasta_rejected(tmp_path, content):
    path = tmp_path / "genes.fa"
    path.write_text(content)
    with pytest.raises(ValueError):
        sequence_inventory(path)


HEADER = "# BUSCO version is: 5.8.2 \n# The lineage dataset is: eukaryota_odb10 (Creation date: 2024-01-08, number of genomes: 70, number of BUSCOs: 255)\n"


def test_complete_filter_preserves_native_reference_rule(tmp_path):
    path = tmp_path / "table.tsv"
    path.write_text(HEADER + "id1\tComplete\ta\n" + "id2\tMissing\n" + "id3\tDuplicated\tb\n")
    assert complete_hits(path, {"a"}) == {"id1": "a"}


@pytest.mark.parametrize("body", ["id1\tComplete\tunknown\n", "id1\tComplete\ta\nid1\tComplete\ta\n", "id1\tComplete\n", "id1\tUnexpected\ta\n"])
def test_bad_complete_rows_rejected(tmp_path, body):
    path = tmp_path / "table.tsv"
    path.write_text(HEADER + body)
    with pytest.raises(ValueError):
        complete_hits(path, {"a"})


def test_lineage_mismatch_rejected(tmp_path):
    path = tmp_path / "table.tsv"
    path.write_text(HEADER.replace("2024-01-08", "2020-01-01"))
    with pytest.raises(ValueError, match="lineage differs"):
        complete_hits(path, set())


def test_stop_removal_is_verified_not_assumed(tmp_path):
    raw, staged = tmp_path / "raw.fa", tmp_path / "staged.fa"
    raw.write_text(">a\nAB*C*\n>b\nABC\n")
    staged.write_text(">a\nABC\n>b\nABD\n")
    rows = change_details(raw, staged, {"a", "b"}, {"a"})
    assert rows[0]["exactly_explained_by_removing_stop_markers"] is True
    assert rows[0]["raw_stop_markers"] == 2
    assert rows[0]["in_scored_reference"] is True
    assert rows[1]["exactly_explained_by_removing_stop_markers"] is False
    assert rows[1]["in_scored_reference"] is False
    assert change_details(raw, staged, set(), set()) == []
