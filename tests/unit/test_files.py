"""Unit tests for orthohmm/files.py (FASTA discovery + output writers)."""
import os

import numpy as np
import pytest

from orthohmm.files import (
    fetch_fasta_files,
    write_clusters_file,
    write_copy_number_file,
    write_fasta_files_for_all_ogs,
    write_fasta_files_for_single_copy_orthologs,
    write_file_of_single_copy_ortholog_names,
)


class TestFetchFastaFiles:
    def test_finds_all_supported_extensions(self, tmp_path):
        for ext in (".fa", ".faa", ".fas", ".fasta", ".pep", ".prot"):
            (tmp_path / f"x{ext}").write_text(">a\nAA\n")
        files = fetch_fasta_files(str(tmp_path))
        assert sorted(files) == sorted(
            ["x.fa", "x.faa", "x.fas", "x.fasta", "x.pep", "x.prot"]
        )

    def test_skips_unsupported_extensions(self, tmp_path):
        (tmp_path / "x.fa").write_text(">a\nAA\n")
        (tmp_path / "x.txt").write_text("not a fasta")
        (tmp_path / "x.fastq").write_text(">a\nAA\n+\nIII\n")
        assert fetch_fasta_files(str(tmp_path)) == ["x.fa"]

    def test_empty_directory_returns_empty_list(self, tmp_path):
        assert fetch_fasta_files(str(tmp_path)) == []

    def test_returns_basenames_only(self, tmp_path):
        (tmp_path / "human.faa").write_text(">a\nAA\n")
        result = fetch_fasta_files(str(tmp_path))
        assert result == ["human.faa"]
        assert "/" not in result[0]


class TestWriteClustersFile:
    def test_one_cluster_one_line(self, tmp_path):
        write_clusters_file(str(tmp_path), [["OG0:", "gene_b", "gene_a", "gene_c"]])
        content = (tmp_path / "orthohmm_orthogroups.txt").read_text()
        # Genes after the OG label are sorted alphabetically
        assert content == "OG0: gene_a gene_b gene_c\n"

    def test_multiple_clusters_one_per_line(self, tmp_path):
        write_clusters_file(
            str(tmp_path),
            [["OG0:", "z", "a"], ["OG1:", "m"]],
        )
        lines = (tmp_path / "orthohmm_orthogroups.txt").read_text().splitlines()
        assert lines == ["OG0: a z", "OG1: m"]


class TestWriteCopyNumberFile:
    def test_header_row_plus_data_rows(self, tmp_path):
        write_copy_number_file(
            str(tmp_path),
            {"files:": ["spA", "spB"], "OG0:": ["1", "2"], "OG1:": ["0", "3"]},
        )
        lines = (tmp_path / "orthohmm_gene_count.txt").read_text().splitlines()
        assert lines[0] == "files: spA spB"
        assert "OG0: 1 2" in lines
        assert "OG1: 0 3" in lines


class TestWriteSingleCopyNames:
    def test_writes_ids_without_trailing_colon(self, tmp_path):
        # Input keys come in like "OG3:" (with trailing colon); writer strips it.
        write_file_of_single_copy_ortholog_names(
            str(tmp_path),
            {"files:": [], "OG3:": [], "OG7:": []},
        )
        content = (tmp_path / "orthohmm_single_copy_orthogroups.txt").read_text()
        assert content == "OG3\nOG7\n"


class TestWriteFastaFilesForAllOgs:
    def test_creates_one_fa_per_og(self, tmp_path):
        ogs = {
            "OG0": [">geneA", "MASTA", ">geneB", "MSTAR"],
            "OG1": [">geneC", "MMM"],
        }
        write_fasta_files_for_all_ogs(str(tmp_path), ogs)
        og_dir = tmp_path / "orthohmm_orthogroups"
        assert og_dir.is_dir()
        assert (og_dir / "OG0.fa").read_text() == ">geneA\nMASTA\n>geneB\nMSTAR\n"
        assert (og_dir / "OG1.fa").read_text() == ">geneC\nMMM\n"


class TestWriteFastaFilesForSingleCopyOrthologs:
    def test_header_gets_taxon_prefix_stripped_extension(self, tmp_path):
        gene_lengths = np.array(
            [
                ("human.fa",   "ENSP1", 100),
                ("mouse.fa",   "ENSM1", 110),
                ("unknown.fa", "X",      50),
            ],
            dtype=[("spp", object), ("name", object), ("length", int)],
        )
        ogs = {
            "OG0": [">ENSP1", "MASTA", ">ENSM1", "MSTAR"],
        }
        write_fasta_files_for_single_copy_orthologs(
            str(tmp_path), ogs, gene_lengths, ["OG0"], extensions=(".fa", ".faa"),
        )
        out = (tmp_path / "orthohmm_single_copy_orthogroups" / "OG0.fa").read_text()
        assert ">human|ENSP1" in out
        assert ">mouse|ENSM1" in out
        # Sequence lines come through untouched
        assert "MASTA" in out
        assert "MSTAR" in out

    def test_unknown_gene_falls_back_to_unknownspecies(self, tmp_path):
        gene_lengths = np.array(
            [("human.fa", "ENSP1", 100)],
            dtype=[("spp", object), ("name", object), ("length", int)],
        )
        ogs = {"OG0": [">UNSEEN_GENE", "AAA"]}
        write_fasta_files_for_single_copy_orthologs(
            str(tmp_path), ogs, gene_lengths, ["OG0"], extensions=(".fa",),
        )
        out = (tmp_path / "orthohmm_single_copy_orthogroups" / "OG0.fa").read_text()
        assert ">UnknownSpecies|UNSEEN_GENE" in out
