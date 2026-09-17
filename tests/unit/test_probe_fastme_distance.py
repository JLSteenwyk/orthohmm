import pytest

from benchmark_tools.probe_fastme_distance import canonical_tree, example_matrices


def test_multiple_matrices():
    text = "4\na 0 1 1 1\nb 1 0 1 1\nc 1 1 0 1\nd 1 1 1 0\n"
    assert len(example_matrices(text + "\n" + text)) == 2
    with pytest.raises(ValueError):
        example_matrices(text + "\n5\n")


def test_root_and_child_order_invariance(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    a.write_text("((a:1,b:1):2,(c:1,d:1):3);")
    b.write_text("(d:1,c:1,(b:1,a:1):5);")
    assert canonical_tree(a, list("abcd")) == canonical_tree(b, list("abcd"))


@pytest.mark.parametrize("tree", ["(a:1,a:1,c:1);", "(a:1,b:1,z:1);", "(a,b:1,c:1);", "(a:1e999,b:1,c:1);"])
def test_bad_tree(tmp_path, tree):
    path = tmp_path / "tree"
    path.write_text(tree)
    with pytest.raises(ValueError):
        canonical_tree(path, list("abc"))
