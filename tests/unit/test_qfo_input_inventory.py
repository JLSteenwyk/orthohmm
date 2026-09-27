import json

import pytest

from benchmark_tools.qfo_input_inventory import record, verify_inventory


def fixture(directory):
    fasta = directory / "species.fasta"
    fasta.write_text(">gene\nMAA\n")
    fastas = [record(fasta)]
    metadata = directory / "staging_manifest.json"
    metadata.write_text(json.dumps(dict(input_fastas=fastas)))
    return fastas, record(metadata)


def test_explicit_metadata_is_admitted(tmp_path):
    fastas, metadata = fixture(tmp_path)
    assert verify_inventory(tmp_path, fastas, metadata) == [*fastas, metadata]


@pytest.mark.parametrize("mutation", ["extra_fasta", "extra_json", "directory",
    "missing_fasta", "missing_metadata", "changed_fasta", "changed_metadata",
    "duplicate", "wrong_metadata_inventory"])
def test_fail_closed(tmp_path, mutation):
    fastas, metadata = fixture(tmp_path)
    if mutation in {"extra_fasta", "extra_json"}:
        (tmp_path / mutation).write_text("unexpected")
    elif mutation == "directory":
        (tmp_path / "extra").mkdir()
    elif mutation.startswith("missing_"):
        (tmp_path / ("species.fasta" if mutation.endswith("fasta") else "staging_manifest.json")).unlink()
    elif mutation.startswith("changed_"):
        (tmp_path / ("species.fasta" if mutation.endswith("fasta") else "staging_manifest.json")).write_text("changed")
    elif mutation == "duplicate":
        fastas.append(fastas[0])
    else:
        path = tmp_path / "staging_manifest.json"
        path.write_text(json.dumps(dict(input_fastas=[])))
        metadata = record(path)
    with pytest.raises(ValueError):
        verify_inventory(tmp_path, fastas, metadata)
