from io import StringIO

import pytest

from benchmark_tools.restore_orthofinder_ids import (
    build_restoration_map,
    restore_orthogroups,
)


def test_restores_orthofinder_sanitized_identifiers(tmp_path):
    (tmp_path / "species.fasta").write_text(
        ">gene:one description\nAAAA\n"
        ">gene,two\nAAAA\n"
        ">gene(three)\nAAAA\n"
        ">unchanged\nAAAA\n"
    )
    restoration = build_restoration_map(tmp_path)
    destination = StringIO()

    counts = restore_orthogroups(
        StringIO("OG0001: gene_one gene_two\nOG0002: gene_three_ unchanged\n"),
        destination,
        restoration,
    )

    assert counts == (2, 4, 3)
    assert destination.getvalue() == (
        "OG0001: gene:one gene,two\nOG0002: gene(three) unchanged\n"
    )


def test_rejects_ambiguous_sanitized_identifiers(tmp_path):
    (tmp_path / "species.faa").write_text(">gene:one\nAAAA\n>gene_one\nAAAA\n")

    with pytest.raises(ValueError, match="sanitization is ambiguous"):
        build_restoration_map(tmp_path)


def test_rejects_unknown_orthofinder_identifier(tmp_path):
    (tmp_path / "species.fa").write_text(">known\nAAAA\n")

    with pytest.raises(ValueError, match="Unknown OrthoFinder identifier"):
        restore_orthogroups(
            StringIO("OG0001: missing\n"),
            StringIO(),
            build_restoration_map(tmp_path),
        )
