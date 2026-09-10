from benchmark_tools.normalize_three_kingdoms_orthogroups import (
    iter_fastoma,
    iter_orthomcl,
    iter_root_hogs,
    iter_sonicparanoid,
    write_groups,
)


def test_normalize_sonicparanoid(tmp_path):
    source = tmp_path / "ortholog_groups.tsv"
    source.write_text(
        "group_id\tgroup_size\tsp_in_grp\tseed_ortholog_cnt\tsp1\tsp2\n"
        "1\t3\t2\t2\ta,b\tc\n"
        "2\t2\t1\t1\t*\td,e\n"
    )

    assert list(iter_sonicparanoid(source)) == [
        ("1", ("a", "b", "c")),
        ("2", ("d", "e")),
    ]


def test_normalize_sonicparanoid_rejects_bad_counts(tmp_path):
    source = tmp_path / "ortholog_groups.tsv"
    source.write_text(
        "group_id\tgroup_size\tsp_in_grp\tseed_ortholog_cnt\tsp1\tsp2\n"
        "1\t4\t2\t2\ta,b\tc\n"
    )

    import pytest

    with pytest.raises(ValueError, match="count mismatch"):
        list(iter_sonicparanoid(source))


def test_normalize_fastoma(tmp_path):
    source = tmp_path / "groups.tsv"
    source.write_text("Group\tProtein\nOG1\ta\nOG1\tb\nOG2\tc\n")

    assert list(iter_fastoma(source)) == [("OG1", ("a", "b")), ("OG2", ("c",))]


def test_normalize_root_hogs(tmp_path):
    source = tmp_path / "root_hogs.tsv"
    source.write_text(
        "root_hog\tsource_family\tgenes\n"
        "RootHOG1\tFamily1\ta,b\n"
        "RootHOG2\tFamily2\tc\n"
    )

    assert list(iter_root_hogs(source)) == [
        ("RootHOG1", ("a", "b")),
        ("RootHOG2", ("c",)),
    ]


def test_normalize_orthomcl_and_write(tmp_path):
    source = tmp_path / "all_orthomcl.out"
    source.write_text(
        "ORTHOMCL0(2 genes,2 taxa): a(Amph) b(Arab)\n"
        "ORTHOMCL1(1 genes,1 taxa): c(Scer)\n"
    )
    output = tmp_path / "orthogroups.txt"

    groups = iter_orthomcl(source)
    assert write_groups(groups, output) == (2, 3)
    assert output.read_text() == "a b\nc\n"
