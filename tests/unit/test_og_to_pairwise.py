from pathlib import Path

import pytest

from qfo_benchmark import og_to_pairwise


def _write_input(input_dir: Path) -> None:
    input_dir.mkdir()
    (input_dir / "species_a.fasta").write_text(">sp|A|GENE_A\nAAAA\n>a2\nAAAA\n")
    (input_dir / "species_b.fasta").write_text(">tr|B|GENE_B\nAAAA\n")


def test_group_to_pairwise_uses_cross_species_pairs(tmp_path, capsys, monkeypatch):
    input_dir = tmp_path / "input"
    _write_input(input_dir)
    groups = tmp_path / "groups.txt"
    groups.write_text("sp|A|GENE_A a2 tr|B|GENE_B\n")
    monkeypatch.setattr(
        og_to_pairwise.sys, "argv", ["og_to_pairwise.py", str(groups), str(input_dir)]
    )

    og_to_pairwise.main()

    captured = capsys.readouterr()
    assert captured.out == "A\tB\nB\ta2\n"
    assert captured.err == "Emitted 2 pairs\n"


def test_group_to_pairwise_rejects_unknown_gene(tmp_path, monkeypatch):
    input_dir = tmp_path / "input"
    _write_input(input_dir)
    groups = tmp_path / "groups.txt"
    groups.write_text("sp|A|GENE_A missing\n")
    monkeypatch.setattr(
        og_to_pairwise.sys, "argv", ["og_to_pairwise.py", str(groups), str(input_dir)]
    )

    with pytest.raises(ValueError, match="absent from"):
        og_to_pairwise.main()


def test_group_to_pairwise_rejects_gene_in_multiple_groups(tmp_path, monkeypatch):
    input_dir = tmp_path / "input"
    _write_input(input_dir)
    groups = tmp_path / "groups.txt"
    groups.write_text("sp|A|GENE_A tr|B|GENE_B\nsp|A|GENE_A a2\n")
    monkeypatch.setattr(
        og_to_pairwise.sys, "argv", ["og_to_pairwise.py", str(groups), str(input_dir)]
    )

    with pytest.raises(ValueError, match="occurs in multiple groups"):
        og_to_pairwise.main()


def test_group_to_pairwise_validates_singletons(tmp_path, monkeypatch):
    input_dir = tmp_path / "input"
    _write_input(input_dir)
    groups = tmp_path / "groups.txt"
    groups.write_text("missing\n")
    monkeypatch.setattr(
        og_to_pairwise.sys, "argv", ["og_to_pairwise.py", str(groups), str(input_dir)]
    )

    with pytest.raises(ValueError, match="absent from"):
        og_to_pairwise.main()
