import pytest

from benchmark_tools.audit_ob_fastoma_provenance import match_input, table_counts, task_records, tree_clades


def test_input_matches_ids_then_checks_residues():
    source = {"first.fa": dict(a="AC", b="DF"), "second.fa": dict(c="HI")}
    assert match_input(dict(b="DF", a="AC"), source) == ("first.fa", [])
    assert match_input(dict(a="AA", b="DF"), source) == ("first.fa", ["a"])
    with pytest.raises(ValueError, match="unique"):
        match_input(dict(a="AC"), source)
    with pytest.raises(ValueError, match="unique"):
        match_input(dict(c="HI"), dict(source, duplicate=dict(c="HI")))


@pytest.mark.parametrize("root_hogs", [False, True])
def test_counts_keep_final_and_root_formats_separate(tmp_path, root_hogs):
    path = tmp_path / "groups.tsv"
    text = "RootHOG\tProtein\tOMAmerRootHOG\ng1\ta\th1\ng1\tb\th1\n" if root_hogs else "Group\tProtein\ng1\ta\ng1\tb\n"
    path.write_text(text)
    assert table_counts(path, {"a", "b", "c"}, root_hogs) == dict(groups=1, assigned_genes=2)
    with pytest.raises(ValueError, match="header"):
        table_counts(path, {"a", "b", "c"}, not root_hogs)


@pytest.mark.parametrize("text", ["Group\tProtein\ng1\tx\n", "Group\tProtein\ng1\ta\ng2\ta\n",
                                  "Group\tProtein\n", "Group\tProtein\ng1\ta\textra\n"])
def test_bad_memberships_rejected(tmp_path, text):
    path = tmp_path / "groups.tsv"
    path.write_text(text)
    with pytest.raises(ValueError):
        table_counts(path, {"a", "b"})


def test_failed_tasks_are_retained(tmp_path):
    for suffix, code in (("first", "0"), ("second", "137")):
        directory = tmp_path / "aa" / suffix
        directory.mkdir(parents=True)
        for name in (".command.sh", ".command.run", ".command.log"):
            (directory / name).write_text("fixture\n")
        (directory / ".exitcode").write_text(code)
    rows = task_records(tmp_path)
    assert [r["exit_code"] for r in rows] == [0, 137]
    assert all(len(r["records"]) == 4 for r in rows)
    (tmp_path / "aa/second/.exitcode").write_text("missing")
    with pytest.raises(ValueError, match="Malformed"):
        task_records(tmp_path)


def test_tree_check_ignores_added_lengths_but_not_topology(tmp_path):
    a, b = tmp_path / "a.nwk", tmp_path / "b.nwk"
    a.write_text("((A,B),C);")
    b.write_text("((A:1,B:1):1,C:1):0;")
    assert tree_clades(a, ["A", "B", "C"]) == tree_clades(b, ["A", "B", "C"])
    b.write_text("((A,C),B);")
    assert tree_clades(a, ["A", "B", "C"]) != tree_clades(b, ["A", "B", "C"])
