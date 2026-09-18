import hashlib
import gzip
import json

import pytest

from benchmark_tools.audit_swiss_additional_sequences import scan
from benchmark_tools import audit_swiss_additional_sequences as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def test_scan_records_literal_isoform_target(tmp_path):
    path = tmp_path / "input.fa"
    path.write_text(">tr|A|ONE Isoform of B, Protein\nACD\n>sp|C|TWO Protein mentioning Isoform of D,\nAAAA\n")
    count, records = scan(path)
    assert count == 2
    assert records["A"]["header_isoform_of"] == "B"
    assert records["C"]["header_isoform_of"] is None
    assert records["A"]["sequence_sha256"] == hashlib.sha256(b"ACD").hexdigest()
    assert scan(path, {"A"})[1] == {"A": records["A"]}
    assert scan(path, {"missing"}) == (2, {})


@pytest.mark.parametrize("text", [">A\nACD\n", ">sp|A|ONE\n", ">sp|A|ONE\nACD\n>tr|A|TWO\nAAA\n"])
def test_reject_invalid_or_duplicate_records(tmp_path, text):
    path = tmp_path / "input.fa"
    path.write_text(text)
    with pytest.raises(ValueError):
        scan(path)


def test_file_audit_records_different_numeric_target_without_remapping(tmp_path, monkeypatch):
    canonical, staged, additional = (tmp_path / name for name in ("canonical.fa", "staged.fa", "additional.fa"))
    canonical.write_text(">sp|P|ONE\nACD\n")
    staged.write_bytes(canonical.read_bytes())
    additional.write_text(">tr|M|TWO Isoform of P, Protein\nACDE\n")
    aliases = tmp_path / "aliases.json"
    aliases.write_text(json.dumps({"fasta_inputs": [record(staged)],
        "summary": {"missing_genes": ["M", "N"]},
        "identities": {"M": {"numeric_protein_id": 1}, "N": {"numeric_protein_id": 3}}}))
    mapping = tmp_path / "mapping.gz"
    with gzip.open(mapping, "wt") as stream:
        json.dump({"mapping": {"M": 1, "P": 2, "N": 3}}, stream)
    monkeypatch.setattr(module, "ALIAS_SHA", record(aliases)["sha256"])
    monkeypatch.setattr(module, "MAPPING_SHA", record(mapping)["sha256"])
    result = module.audit(aliases, canonical, additional, staged, mapping)
    assert result["found_in_additional"] == 1
    assert result["genes"]["M"]["numeric_protein_id"] == 1
    assert result["genes"]["M"]["header_target_numeric_id"] == 2
    assert result["genes"]["M"]["header_target_in_canonical"]
    assert not result["genes"]["N"]["found_in_additional"]
    canonical.write_text(">sp|P|ONE\nAAA\n")
    with pytest.raises(ValueError, match="FASTAs differ"):
        module.audit(aliases, canonical, additional, staged, mapping)
