from pathlib import Path
import shutil
import subprocess
import sys

from Bio import SeqIO
import pytest


@pytest.mark.integration
def test_module_cli_native_partition(tmp_path):
    root = Path(__file__).resolve().parents[2]
    source = root / "tests" / "samples"
    inputs = tmp_path / "input"
    outputs = tmp_path / "output"
    inputs.mkdir()
    outputs.mkdir()
    genes = []
    original_bytes = {}
    for path in sorted(source.iterdir()):
        if path.is_file() and path.suffix in (".fa", ".faa", ".fas", ".fasta", ".pep", ".prot"):
            original_bytes[path] = path.read_bytes()
            shutil.copyfile(path, inputs / path.name)
            genes.extend(record.id for record in SeqIO.parse(path, "fasta"))
    assert genes and len(genes) == len(set(genes))
    response = subprocess.run(
        [sys.executable, "-m", "orthohmm", str(inputs), "-o", str(outputs),
         "-c", "1", "--search_mode", "builtin", "--clustering", "leiden"],
        cwd=root, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, timeout=180, check=False,
    )
    assert response.returncode == 0, response.stdout
    assigned = []
    for line in (outputs / "orthohmm_orthogroups.txt").read_text().splitlines():
        group, members = line.split(":", 1)
        assert group and members.split()
        assigned.extend(members.split())
    assert sorted(assigned) == sorted(genes)
    for path, content in original_bytes.items():
        assert path.read_bytes() == content
        assert (inputs / path.name).read_bytes() == content
