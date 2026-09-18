import pytest
import gzip
import json

from benchmark_tools import resolve_swiss_sequence_aliases as module
from benchmark_tools.resolve_swiss_sequence_aliases import resolve
from benchmark_tools.inventory_swiss_sequences import collect
from benchmark_tools.snapshot_orthohmm_input_order import record


def test_exact_alias_missing_and_ambiguity():
    mapping = {"R1": 1, "R2": 2, "R3": 3, "R4": 4,
               "A2": 2, "A4": 4, "B4": 4}
    rows = resolve(["R1", "R2", "R3", "R4"], ["R1", "A2", "B4", "A4", "unmapped"], mapping)
    assert [rows[f"R{i}"]["status"] for i in range(1, 5)] == ["exact", "mapped_alias", "missing", "ambiguous"]
    assert rows["R4"]["input_accessions"] == ["A4", "B4"]
    assert rows["R2"]["numeric_protein_id"] == 2


def test_exact_does_not_hide_ambiguous_alias():
    assert resolve(["R"], ["R", "A"], {"R": 1, "A": 1})["R"]["status"] == "ambiguous"


@pytest.mark.parametrize("number", [None, True, "1", 0, -1, 1.5])
def test_invalid_reference_number(number):
    with pytest.raises(ValueError):
        resolve(["R"], [], {"R": number})


def test_shared_reference_identity():
    with pytest.raises(ValueError, match="Shared"):
        resolve(["R", "S"], [], {"R": 1, "S": 1})


def test_duplicate_input_identity():
    with pytest.raises(ValueError, match="Duplicate"):
        resolve(["R"], ["R", "R"], {"R": 1})


def test_invalid_input_number():
    with pytest.raises(ValueError, match="Invalid input"):
        resolve(["R"], ["A"], {"R": 1, "A": True})


def test_file_audit_preserves_exact_and_resolves_alias(tmp_path, monkeypatch):
    fasta = tmp_path / "input.fa"
    fasta.write_text(">sp|R1|ONE\nACD\n>tr|A2|TWO name (Fragment)\nAAAA\n")
    families = {"F": ["R1", "R2", "R3"]}
    inputs = [record(fasta)]
    counts = tmp_path / "counts.json"
    counts.write_text(json.dumps({"methods": [{"families": [{"family": "F", "represented_genes": families["F"]}]}]}))
    prepared = tmp_path / "prepared.json"
    prepared.write_text(json.dumps({"input_fastas": inputs}))
    inventory = tmp_path / "inventory.json"
    inventory.write_text(json.dumps({**collect(families, inputs), "inputs": [record(counts)], "fasta_inputs": inputs}))
    mapping = tmp_path / "mapping.json.gz"
    with gzip.open(mapping, "wt") as stream:
        json.dump({"mapping": {"R1": 1, "R2": 2, "A2": 2, "R3": 3}}, stream)
    for name, path in (("INVENTORY_SHA", inventory), ("PREPARED_SHA", prepared), ("MAPPING_SHA", mapping)):
        monkeypatch.setattr(module, name, record(path)["sha256"])
    result = module.audit(inventory, prepared, mapping)
    assert result["resolution_counts"] == {"exact": 1, "mapped_alias": 1, "missing": 1, "ambiguous": 0}
    assert result["genes"]["R2"]["resolved_input_accession"] == "A2"
    assert result["summary"]["explicit_fragment_descriptions"] == 1
    assert result["summary"]["missing_genes"] == ["R3"]
    assert not result["prediction_statistics_evaluated"]
    counts.write_text("{}")
    with pytest.raises(ValueError, match="Changed reference"):
        module.audit(inventory, prepared, mapping)
