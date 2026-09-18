import pytest

from tests.integration.output_checks import check_outputs


def fixture(tmp_path):
    inputs, expected, output = (tmp_path / name for name in ("input", "expected", "output"))
    for directory in (inputs, expected, output):
        directory.mkdir()
    (inputs / "a.fa").write_text(">a1\nAAA\n>a2\nCCC\n")
    (inputs / "b.faa").write_text(">b1\nBBB\n")
    (expected / "orthohmm_orthogroups.txt").write_text("OG000: a1 b1\nOG001: a2\n")
    (output / "orthohmm_orthogroups.txt").write_text("OG9: b1 a1\nOG3: a2\n")
    (output / "orthohmm_gene_count.txt").write_text("files: b.faa a.fa\nOG3: 0 1\nOG9: 1 1\n")
    (output / "orthohmm_single_copy_orthogroups.txt").write_text("OG9\n")
    all_groups = output / "orthohmm_orthogroups"
    single = output / "orthohmm_single_copy_orthogroups"
    all_groups.mkdir()
    single.mkdir()
    (all_groups / "OG9.fa").write_text(">b1\nBBB\n>a1\nAAA\n")
    (all_groups / "OG3.fa").write_text(">a2\nCCC\n")
    (single / "OG9.fa").write_text(">a|a1\nAAA\n>b|b1\nBBB\n")
    return output, expected, list(inputs.iterdir())


def test_species_group_and_record_order_do_not_change_meaning(tmp_path):
    assert check_outputs(*fixture(tmp_path)) == (2, 1)


@pytest.mark.parametrize("problem", ["count", "columns", "sequence", "duplicate_fasta", "prefix",
                                    "all_listed", "missing_single", "extra_file", "partition", "duplicate_gene"])
def test_actual_errors_are_not_normalized_away(tmp_path, problem):
    output, expected, inputs = fixture(tmp_path)
    if problem == "count":
        path = output / "orthohmm_gene_count.txt"
        path.write_text(path.read_text().replace("OG9: 1 1", "OG9: 2 1"))
    elif problem == "columns":
        path = output / "orthohmm_gene_count.txt"
        path.write_text(path.read_text().replace("b.faa a.fa", "a.fa a.fa"))
    elif problem in ("sequence", "duplicate_fasta"):
        path = output / "orthohmm_orthogroups/OG9.fa"
        path.write_text(path.read_text().replace("AAA", "AXA") if problem == "sequence"
                        else path.read_text() + ">a1\nAAA\n")
    elif problem == "prefix":
        path = output / "orthohmm_single_copy_orthogroups/OG9.fa"
        path.write_text(path.read_text().replace("a|a1", "wrong|a1"))
    elif problem in ("all_listed", "missing_single"):
        (output / "orthohmm_single_copy_orthogroups.txt").write_text(
            "OG9\nOG3\n" if problem == "all_listed" else "")
    elif problem == "extra_file":
        (output / "orthohmm_orthogroups/extra.fa").write_text(">x\nA\n")
    elif problem == "partition":
        (output / "orthohmm_orthogroups.txt").write_text("OG9: b1 a2\nOG3: a1\n")
    else:
        (output / "orthohmm_orthogroups.txt").write_text("OG9: b1 a1\nOG3: a1 a2\n")
    with pytest.raises(AssertionError):
        check_outputs(output, expected, inputs)
