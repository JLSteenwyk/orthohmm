import json

import pytest

from benchmark_tools.prepare_orthomcl_inputs import prepare_inputs


def test_prepare_inputs_matches_native_header_and_taxon_format(tmp_path):
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    (input_dir / "Beta.fasta").write_bytes(b">b1 description\r\nMKT\r\n")
    (input_dir / "Alpha.fasta").write_bytes(b">a1 first description\nMA\nTG\n>a2\nQQ\n")
    combined = tmp_path / "all.fa"
    gg = tmp_path / "all.gg"
    summary_path = tmp_path / "summary.json"

    summary = prepare_inputs(input_dir, combined, gg, summary_path)

    assert combined.read_bytes() == b">a1\nMA\nTG\n>a2\nQQ\n>b1\nMKT\n"
    assert gg.read_bytes() == b"Alpha: a1 a2\nBeta: b1\n"
    assert summary["taxon_count"] == 2
    assert summary["sequence_count"] == 3
    assert json.loads(summary_path.read_text()) == summary


def test_prepare_inputs_rejects_duplicate_ids_without_publishing(tmp_path):
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    (input_dir / "A.fasta").write_text(">same\nAA\n")
    (input_dir / "B.fasta").write_text(">same duplicate\nBB\n")
    combined = tmp_path / "all.fa"
    gg = tmp_path / "all.gg"

    with pytest.raises(ValueError, match="Duplicate FASTA ID: same"):
        prepare_inputs(input_dir, combined, gg)

    assert not combined.exists()
    assert not gg.exists()
