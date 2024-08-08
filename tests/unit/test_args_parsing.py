from argparse import Namespace
import pytest
from pathlib import Path

from orthohmm.args_processing import process_args


here = Path(__file__)


@pytest.fixture
def args():
    kwargs = dict(
        fasta_directory=f"{here.parent.parent}/samples/",
        output_directory=f"{here.parent.parent}/samples/",
        phmmer="phmmer",
        cpu=8,
        single_copy_threshold=0.5,
        mcl="mcl",
        inflation_value=1.5,
    )
    return Namespace(**kwargs)


class TestArgsProcessing(object):
    def test_process_args_input_directory_dne(self, args):
        args.fasta_directory = "some/file/that/doesnt/exist"
        with pytest.raises(SystemExit):
            process_args(args)

    def test_process_args_output_directory_dne(self, args):
        args.output_directory = "some/file/that/doesnt/exist"
        with pytest.raises(SystemExit):
            process_args(args)

    def test_process_args_phmmer_not_installed(self, args):
        args.phmmer = "phmmer-that-dne"
        with pytest.raises(SystemExit):
            process_args(args)

    def test_process_args_mcl_not_installed(self, args):
        args.mcl = "mcl-that-dne"
        with pytest.raises(SystemExit):
            process_args(args)

    def test_process_args_default_single_copy_threshold(self, args):
        args.single_copy_threshold = None
        res = process_args(args)
        assert res["single_copy_threshold"] == 0.5

    def test_process_args_default_inflation_value(self, args):
        args.inflation_value = None
        res = process_args(args)
        assert res["inflation_value"] == 1.5

    def test_process_args_expected_keywords(self, args):
        res = process_args(args)
        expected_keys = [
            "fasta_directory",
            "output_directory",
            "phmmer",
            "cpu",
            "single_copy_threshold",
            "mcl",
            "inflation_value",
        ]
        assert sorted(res.keys()) == sorted(expected_keys)
