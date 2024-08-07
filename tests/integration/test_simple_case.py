import os
import pytest
from pathlib import Path

from orthohmm.orthohmm import execute


here = Path(__file__)


@pytest.mark.integration
class TestSimpleCase(object):
    @pytest.mark.slow
    def test_simple_case(self):
        """
        test simple case
        usage: orthohmm ./path_to_fasta_directory
        """
        kwargs = dict(
            fasta_directory=f"{here.parent.parent}/samples/",
            output_directory=f"{here.parent.parent}/samples/",
            phmmer="phmmer",
            cpu=8,
            single_copy_threshold=0.5,
            mcl="mcl",
            inflation_value=1.5,
        )

        execute(**kwargs)

        # read orthohmm_gene_count
        with open(f"{here.parent.parent}/expected/orthohmm_gene_count.txt", "r") as expected:
            expected_orthohmm_gene_count = expected.read()

        with open(f"{here.parent.parent}/samples/orthohmm_gene_count.txt", "r") as out_file:
            output_orthohmm_gene_count = out_file.read()

        # read orthogroups
        with open(f"{here.parent.parent}/expected/orthohmm_orthogroups.txt", "r") as expected:
            expected_orthohmm_orthogroups = expected.read()

        with open(f"{here.parent.parent}/samples/orthohmm_orthogroups.txt", "r") as out_file:
            output_orthohmm_orthogroups = out_file.read()

        # read single-copy orthogroups
        with open(f"{here.parent.parent}/expected/orthohmm_single_copy_orthogroups.txt", "r") as expected:
            expected_orthohmm_single_copy_orthogroups = expected.read()

        with open(f"{here.parent.parent}/samples/orthohmm_single_copy_orthogroups.txt", "r") as out_file:
            output_orthohmm_single_copy_orthogroups = out_file.read()

        # determine how many files were made in single-copy and all orthogroup directories
        files_orthohmm_single_copy_orthogroups = [
            f for f in os.listdir(f"{here.parent.parent}/samples/orthohmm_single_copy_orthogroups/") if os.path.isfile(os.path.join(f"{here.parent.parent}/samples/orthohmm_single_copy_orthogroups/", f))
        ]

        files_orthohmm_orthogroups = [
            f for f in os.listdir(f"{here.parent.parent}/samples/orthohmm_orthogroups/") if os.path.isfile(os.path.join(f"{here.parent.parent}/samples/orthohmm_orthogroups/", f))
        ]

        # read each individual ortholog fasta
        with open(f"{here.parent.parent}/expected/orthohmm_orthogroups/OG0.fa", "r") as expected:
            expected_OG0 = expected.read()

        with open(f"{here.parent.parent}/samples/orthohmm_orthogroups/OG0.fa", "r") as out_file:
            output_OG0 = out_file.read()

        with open(f"{here.parent.parent}/expected/orthohmm_orthogroups/OG1.fa", "r") as expected:
            expected_OG1 = expected.read()

        with open(f"{here.parent.parent}/samples/orthohmm_orthogroups/OG1.fa", "r") as out_file:
            output_OG1 = out_file.read()

        with open(f"{here.parent.parent}/expected/orthohmm_orthogroups/OG2.fa", "r") as expected:
            expected_OG2 = expected.read()

        with open(f"{here.parent.parent}/samples/orthohmm_orthogroups/OG2.fa", "r") as out_file:
            output_OG2 = out_file.read()

        with open(f"{here.parent.parent}/expected/orthohmm_orthogroups/OG3.fa", "r") as expected:
            expected_OG3 = expected.read()

        with open(f"{here.parent.parent}/samples/orthohmm_orthogroups/OG3.fa", "r") as out_file:
            output_OG3 = out_file.read()

        with open(f"{here.parent.parent}/expected/orthohmm_orthogroups/OG4.fa", "r") as expected:
            expected_OG4 = expected.read()

        with open(f"{here.parent.parent}/samples/orthohmm_orthogroups/OG4.fa", "r") as out_file:
            output_OG4 = out_file.read()

        with open(f"{here.parent.parent}/expected/orthohmm_single_copy_orthogroups/OG2.fa", "r") as expected:
            expected_OG2_sc = expected.read()

        with open(f"{here.parent.parent}/samples/orthohmm_single_copy_orthogroups/OG2.fa", "r") as out_file:
            output_OG2_sc = out_file.read()

        assert expected_orthohmm_gene_count == output_orthohmm_gene_count
        assert expected_orthohmm_orthogroups == output_orthohmm_orthogroups
        assert expected_orthohmm_single_copy_orthogroups == output_orthohmm_single_copy_orthogroups

        assert os.path.isdir(f"{here.parent.parent}/samples/orthohmm_single_copy_orthogroups/")
        assert os.path.isdir(f"{here.parent.parent}/samples/orthohmm_orthogroups/")

        assert len(files_orthohmm_single_copy_orthogroups) == 1
        assert len(files_orthohmm_orthogroups) == 5

        assert expected_OG0 == output_OG0
        assert expected_OG1 == output_OG1
        assert expected_OG2 == output_OG2
        assert expected_OG3 == output_OG3
        assert expected_OG4 == output_OG4
        assert expected_OG2_sc == output_OG2_sc
