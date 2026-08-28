"""Unit tests for orthohmm/writer.py (banner + stats printers)."""
import re
import time

import pytest

from orthohmm.helpers import StartStep, StopStep, SubstitutionMatrix
from orthohmm.writer import write_output_stats, write_user_args


class TestWriteUserArgs:
    def test_prints_builtin_engine_by_default(self, capsys):
        write_user_args(
            fasta_directory="/in",
            output_directory="/out",
            phmmer="phmmer",
            mcl="mcl",
            cpu=8,
            single_copy_threshold=0.5,
            files=["a.fa", "b.fa"],
            start=None,
            stop=None,
            substitution_matrix=SubstitutionMatrix.blosum62,
            evalue_threshold=1e-4,
            inflation_value=1.5,
        )
        out = capsys.readouterr().out
        assert "built-in" in out
        assert "Accuracy profile: standard" in out
        assert "BLOSUM62" in out
        assert "Clustering: Leiden CPM (resolution=0.1)" in out
        assert "Path to mcl" not in out
        assert "Number of FASTA files: 2" in out
        assert "Step to start analysis: NA" in out
        assert "Step to stop analysis: NA" in out

    def test_phmmer_mode_shows_phmmer_path(self, capsys):
        write_user_args(
            fasta_directory="/in",
            output_directory="/out",
            phmmer="/path/to/phmmer",
            mcl="mcl",
            cpu=4,
            single_copy_threshold=0.5,
            files=[],
            start=None,
            stop=None,
            substitution_matrix=SubstitutionMatrix.blosum62,
            evalue_threshold=1e-4,
            inflation_value=1.5,
            search_mode="phmmer",
        )
        out = capsys.readouterr().out
        assert "phmmer (/path/to/phmmer)" in out

    def test_high_sensitivity_accuracy_profile_is_printed(self, capsys):
        write_user_args(
            fasta_directory="/in",
            output_directory="/out",
            phmmer="phmmer",
            mcl="mcl",
            cpu=4,
            single_copy_threshold=0.5,
            files=[],
            start=None,
            stop=None,
            substitution_matrix=SubstitutionMatrix.blosum62,
            evalue_threshold=1e-4,
            inflation_value=1.5,
            accuracy_profile="high_sensitivity",
        )
        out = capsys.readouterr().out
        assert "Accuracy profile: high_sensitivity" in out

    def test_mcl_mode_shows_mcl_settings(self, capsys):
        write_user_args(
            fasta_directory="/in",
            output_directory="/out",
            phmmer="phmmer",
            mcl="/path/to/mcl",
            cpu=4,
            single_copy_threshold=0.5,
            files=[],
            start=None,
            stop=None,
            substitution_matrix=SubstitutionMatrix.blosum62,
            evalue_threshold=1e-4,
            inflation_value=2.0,
            clustering="mcl",
        )
        out = capsys.readouterr().out
        assert "Clustering: MCL (/path/to/mcl, inflation=2.0)" in out

    def test_start_and_stop_step_names_are_printed(self, capsys):
        write_user_args(
            fasta_directory="/in",
            output_directory="/out",
            phmmer="phmmer",
            mcl="mcl",
            cpu=1,
            single_copy_threshold=0.7,
            files=["a"],
            start=StartStep.search_res,
            stop=StopStep.infer,
            substitution_matrix=SubstitutionMatrix.blosum62,
            evalue_threshold=1e-4,
            inflation_value=1.5,
        )
        out = capsys.readouterr().out
        assert "search_res" in out
        assert "infer" in out


class TestWriteOutputStats:
    def test_counts_match_inputs(self, capsys):
        write_output_stats(
            start_time=time.time() - 0.123,
            single_copy_ogs=["OG0", "OG1"],
            singletons=["g1", "g2", "g3"],
            ogs_dat={"OG0": [], "OG1": [], "OG2": [], "OG3": [], "OG4": []},
            edges={frozenset(("a", "b")): 0.1, frozenset(("c", "d")): 0.2},
            gene_lengths=list(range(17)),
        )
        out = capsys.readouterr().out
        assert "Number of genes processed: 17" in out
        assert "Number of orthogroups: 5" in out
        assert "Number of edges in network: 2" in out
        assert "Number of single-copy orthogroups: 2" in out
        assert "Number of singletons: 3" in out
        # Execution time is rendered to 3 decimal places
        assert re.search(r"Execution time: \d+\.\d{1,3}s", out)
