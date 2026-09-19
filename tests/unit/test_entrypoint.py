import subprocess
import sys
from pathlib import Path
from unittest.mock import Mock

from orthohmm import orthohmm


def run_module(*args):
    return subprocess.run(
        [sys.executable, "-m", "orthohmm", *map(str, args)],
        cwd=Path(__file__).resolve().parents[2],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, timeout=60, check=False,
    )


class TestEntrypoint:
    def test_help(self):
        response = run_module("--help")
        assert response.returncode == 0
        assert "Usage: orthohmm" in response.stdout

    def test_run(self, tmp_path, monkeypatch):
        inputs = tmp_path / "input"
        outputs = tmp_path / "output"
        inputs.mkdir()
        outputs.mkdir()
        for species in ("species_a", "species_b"):
            (inputs / f"{species}.faa").write_text(f">{species}_gene\nMPEPTIDE\n")
        execute = Mock()
        monkeypatch.setattr(orthohmm, "execute", execute)
        monkeypatch.setattr(sys, "argv", [
            "orthohmm", str(inputs), "-o", str(outputs), "-c", "1",
            "--search_mode", "builtin", "--clustering", "leiden",
        ])
        orthohmm.main()
        execute.assert_called_once()
        args = execute.call_args.kwargs
        assert args["fasta_directory"] == str(inputs)
        assert args["output_directory"] == str(outputs)
        assert args["cpu"] == 1
        assert args["search_mode"] == "builtin"
        assert args["clustering"] == "leiden"
        assert list(outputs.iterdir()) == []

    def test_input_error(self, tmp_path):
        response = run_module(tmp_path / "missing")
        assert response.returncode == 0
        assert response.stdout == "Input directory does not exist\n"

    def test_run_no_args(self):
        response = run_module()
        assert response.returncode == 0
        assert "Usage: orthohmm" in response.stdout
