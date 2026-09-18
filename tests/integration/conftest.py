import os
from pathlib import Path
import shutil

import pytest

from orthohmm.helpers import SubstitutionMatrix
from orthohmm.orthohmm import execute


@pytest.fixture
def native_case(tmp_path):
    def run(suffix, reverse=False):
        base = Path(__file__).resolve().parents[1]
        source = base / "samples" / suffix
        inputs = tmp_path / "input"
        output = tmp_path / "output"
        inputs.mkdir()
        output.mkdir()
        files = sorted((p for p in source.iterdir() if p.is_file()
                        and p.suffix in (".fa", ".faa", ".fas", ".fasta", ".pep", ".prot")),
                       reverse=reverse)
        for path in files:
            shutil.copyfile(path, inputs / path.name)
        execute(fasta_directory=str(inputs), output_directory=str(output),
                phmmer=os.environ.get("ORTHOHMM_TEST_PHMMER", "phmmer"), cpu=8,
                single_copy_threshold=.5, mcl=os.environ.get("ORTHOHMM_TEST_MCL", "mcl"),
                inflation_value=1.5, start=None, stop=None,
                substitution_matrix=SubstitutionMatrix.blosum62, evalue_threshold=.0001,
                search_mode="phmmer", clustering="mcl")
        return output, base / "expected" / suffix, list(inputs.iterdir())
    return run
