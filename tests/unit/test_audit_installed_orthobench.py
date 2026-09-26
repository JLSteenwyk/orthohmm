import pytest

from benchmark_tools.audit_installed_orthobench import (
    compare_partitions, read_root_hogs, verify_input_inventory,
)
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def test_strict_partition_and_label_invariance(tmp_path):
    path = tmp_path / "groups.tsv"
    path.write_text("root_hog\tsource_family\tgenes\nH1\tF1\ta,b\nH2\tF2\tc\n")
    groups = read_root_hogs(path, {"a", "b", "c"})
    assert compare_partitions(groups, groups[::-1])["label_invariant_equal"]
    split = [frozenset("a"), frozenset("b"), frozenset("c")]
    result = compare_partitions(groups, split)
    assert result["genes_in_changed_groups"] == 2
    assert result["identical_groups"] == 1


@pytest.mark.parametrize("body", [
    "H1\tF1\ta,a\nH2\tF2\tb\n", "H1\tF1\ta,b\nH2\tF2\ta\n",
    "H1\tF1\ta\n", "H1\tF1\ta,b,x\n", "H1\tF1\ta\nH1\tF1\tb\n",
    "H1\t\ta,b\n", "H1\tF1\ta,b,\n", "H1\tF1\ta,b\textra\n", "H1\tF1\n",
])
def test_invalid_memberships_rejected(tmp_path, body):
    path = tmp_path / "groups.tsv"
    path.write_text("root_hog\tsource_family\tgenes\n" + body)
    with pytest.raises(ValueError):
        read_root_hogs(path, {"a", "b"})


@pytest.fixture
def inventory(tmp_path):
    directory = tmp_path / "input"
    directory.mkdir()
    path = directory / "one.fa"
    path.write_text(">a\nACDE\n>b\nFGHI\n")
    plan = dict(checked_records=[record(path)], expected_species=1, expected_genes=2)
    return directory, path, plan


def test_inventory_exact(inventory):
    directory, _, plan = inventory
    genes, records = verify_input_inventory(directory, plan)
    assert genes == {"a", "b"}
    assert records == plan["checked_records"]


@pytest.mark.parametrize("name", ["extra.fa", "extra.faa", "extra.fas", "extra.fasta",
                                  "extra.pep", "extra.prot", ".hidden.fa", "notes.txt"])
def test_inventory_extra_file(inventory, name):
    directory, _, plan = inventory
    (directory / name).write_text(">extra\nACDE\n")
    with pytest.raises(ValueError, match="inventory differs"):
        verify_input_inventory(directory, plan)


@pytest.mark.parametrize("change", ["missing", "content", "symlink", "directory",
                                   "duplicate_record", "wrong_gene_count", "wrong_species_count"])
def test_inventory_mutations(inventory, change):
    directory, path, plan = inventory
    if change == "missing":
        path.unlink()
    elif change == "content":
        path.write_text(">a\nAAAA\n>b\nFGHI\n")
    elif change == "symlink":
        target = directory.parent / "original.fa"
        path.rename(target)
        path.symlink_to(target)
    elif change == "directory":
        path.unlink()
        path.mkdir()
    elif change == "duplicate_record":
        plan["checked_records"] *= 2
        plan["expected_species"] = 2
    elif change == "wrong_gene_count":
        plan["expected_genes"] = 3
    else:
        plan["expected_species"] = 2
    with pytest.raises(ValueError):
        verify_input_inventory(directory, plan)
