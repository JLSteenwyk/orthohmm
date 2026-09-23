import pytest

from benchmark_tools.acquire_swisstree_models import model_path


def test_script_model_only_not_pool_or_raw():
    page = "<script>render('view','/ST/ST001/modeltree.nhx','/ST/ST001/pooltrees.nhx');</script>"
    assert model_path(page, "ST001") == "/ST/ST001/modeltree.nhx"


def test_onclick_tree_and_commented_alternative():
    page = ('<a onclick="openWin(\'/ST/ST013/ST013_treemodel.phyloxml\')">Reference</a>'
            '<!-- <a href="/ST/ST013/modeltree.nhx">Old</a> -->')
    assert model_path(page, "ST013") == "/ST/ST013/ST013_treemodel.phyloxml"


@pytest.mark.parametrize("page,identifier", [
    ("<a href='/ST/ST001/modeltree.nhx'>Wrong family</a>", "ST002"),
    ("<a href='/ST/ST001/modeltree.nhx.bad'>Wrong suffix</a>", "ST001"),
    ("<a href='/ST/ST001/modeltree.nhx'>a</a><a href='/ST/ST001/ST001_treemodel.phyloxml'>b</a>", "ST001"),
    ("No link", "ST001"), ("No link", "ST020")])
def test_missing_or_ambiguous_model_rejected(page, identifier):
    with pytest.raises(ValueError):
        model_path(page, identifier)
