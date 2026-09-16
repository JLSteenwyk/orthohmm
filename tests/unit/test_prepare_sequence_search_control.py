from pathlib import Path

import pytest

from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.prepare_sequence_search_control import search_command, write_queries


def test_explicit_search_settings_and_raw_score_output():
    argv = search_command("diamond", "queries", "db", "hits")
    for key, value in {"--threads": "32", "--max-target-seqs": "0", "--max-hsps": "1",
                       "--matrix": "BLOSUM62", "--gapopen": "11", "--gapextend": "1",
                       "--comp-based-stats": "1", "--masking": "1", "--evalue": "0.0001"}.items():
        assert argv[argv.index(key) + 1] == value
    assert argv[argv.index("--outfmt") + 1:] == ["6", "qseqid", "sseqid", "qlen", "slen", "score", "bitscore", "evalue"]
    assert "--very-sensitive" in argv
    assert "--no-self-hits" not in argv


def test_query_preparation_preserves_ids_lengths_and_sequences(tmp_path):
    fasta = tmp_path / "species.fasta"
    fasta.write_text(">gene1 description\nACDE\n>gene2\nKLM\n")
    output = tmp_path / "query.fasta"
    owners, lengths = write_queries([file_provenance(fasta)], output)
    assert owners == {"gene1": "species.fasta", "gene2": "species.fasta"}
    assert lengths == {"gene1": 4, "gene2": 3}
    assert output.read_text() == fasta.read_text()
    with pytest.raises(FileExistsError):
        write_queries([file_provenance(fasta)], output)


def test_duplicate_ids_rejected(tmp_path):
    fasta = tmp_path / "species.fasta"
    fasta.write_text(">gene\nACDE\n>gene\nACDE\n")
    with pytest.raises(ValueError, match="Duplicate"):
        write_queries([file_provenance(fasta)], tmp_path / "output")


def test_changed_input_rejected(tmp_path):
    fasta = tmp_path / "species.fasta"
    fasta.write_text(">gene\nACDE\n")
    record = file_provenance(fasta)
    fasta.write_text(">gene\nKLMN\n")
    with pytest.raises(ValueError, match="Changed"):
        write_queries([record], tmp_path / "output")
