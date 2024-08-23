from argparse import Namespace
import pytest
from pathlib import Path

from orthohmm.args_processing import process_args
from orthohmm.helpers import (
    StartStep,
    StopStep,
    SubstitutionMatrix,
)


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
        start=None,
        stop=None,
        substitution_matrix=SubstitutionMatrix.blosum62
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

    def test_process_args_default_start(self, args):
        args.start = None
        res = process_args(args)
        assert res["start"] is None

    def test_process_args_search_res_start(self, args):
        args.start = "search_res"
        res = process_args(args)
        assert res["start"] is StartStep.search_res

    def test_process_args_substitution_matrix_default(self, args):
        args.substitution_matrix = "BLOSUM62"
        res = process_args(args)
        assert res["substitution_matrix"] is SubstitutionMatrix.blosum62

    def test_process_args_default_stop(self, args):
        args.stop = None
        res = process_args(args)
        assert res["stop"] is None

    def test_process_args_prepare_stop(self, args):
        args.stop = "prepare"
        res = process_args(args)
        assert res["stop"] == StopStep.prepare

    def test_process_args_infer_stop(self, args):
        args.stop = "infer"
        res = process_args(args)
        assert res["stop"] == StopStep.infer

    def test_process_args_write_stop(self, args):
        args.stop = "write"
        res = process_args(args)
        assert res["stop"] == StopStep.write

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
            "start",
            "stop",
            "substitution_matrix",
        ]
        assert sorted(res.keys()) == sorted(expected_keys)
