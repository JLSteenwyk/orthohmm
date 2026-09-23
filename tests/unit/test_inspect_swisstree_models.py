import pytest

from benchmark_tools.inspect_swisstree_models import inspect_tree


def test_nhx_unknown_events_and_mapping_candidates(tmp_path):
    path = tmp_path / "model.nhx"
    path.write_text("(X_P12345,Y_P12345,Z_P23456)[&&NHX:D=Y];")
    result = inspect_tree(path, "nhx", ["P12345", "P23456", "MISSING"])
    assert result["explicit_duplications"] is None
    assert result["rooted_attribute"] is None
    assert result["ambiguous_label_candidates"] == {"P12345": ["X_P12345", "Y_P12345"]}
    assert result["exact_label_candidates"] == {"P23456": "Z_P23456"}
    assert result["missing_label_candidates"] == ["MISSING"]
    assert result["mapping_admitted"] is False


def test_xml_explicit_events_not_implicit_speciations(tmp_path):
    path = tmp_path / "model.xml"
    path.write_text('<phyloxml xmlns="http://www.phyloxml.org"><phylogeny rooted="false"><clade>'
        '<events><duplications>1</duplications></events><clade><name>A</name></clade>'
        '<clade><name>B</name></clade></clade></phylogeny></phyloxml>')
    result = inspect_tree(path, "phyloxml", ["A", "B"])
    assert result["explicit_duplications"] == result["nodes_with_explicit_duplication_count"] == 1
    assert result["rooted_attribute"] is False


def test_duplicate_leaves_rejected(tmp_path):
    path = tmp_path / "model.nhx"
    path.write_text("(A,A);")
    with pytest.raises(ValueError, match="duplicate"):
        inspect_tree(path, "nhx", ["A"])
