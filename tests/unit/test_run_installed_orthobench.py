from pathlib import Path

import pytest

from benchmark_tools.run_installed_orthobench import command, fasta_ids, run


def test_command_is_fresh_isolated_and_inferred():
    argv = command("python", "input", "output", "mafft", "FastTree")
    assert argv[:4] == ["python", "-I", "-m", "orthohmm"]
    assert argv[argv.index("--species_tree_mode") + 1] == "infer"
    assert argv[argv.index("--phylogeny_candidates") + 1] == "satellite_v2"
    assert argv[argv.index("-c") + 1] == "32"
    assert not any("checkpoint" in value or "resume" in value for value in argv)


def test_fasta_uniqueness(tmp_path):
    p = tmp_path / "input.fa"
    p.write_text(">a description\nAA\n>b\nBB\n")
    assert fasta_ids([p]) == {"a", "b"}
    with pytest.raises(ValueError):
        fasta_ids([p, p])
    p.write_text("> \nAA\n")
    with pytest.raises(ValueError):
        fasta_ids([p])


def test_run_requires_allocation(tmp_path, monkeypatch):
    p = tmp_path / "plan.json"
    p.write_text("{}")
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    with pytest.raises(ValueError, match="Slurm"):
        run(p)
