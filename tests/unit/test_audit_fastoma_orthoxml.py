import gzip
import xml.etree.ElementTree as ET

import pytest

from benchmark_tools.audit_fastoma_orthoxml import read_xml, check_root_table, check_pair_scope

XML = '''<orthoXML xmlns="http://orthoXML.org/2011/" origin="FastOMA 0.3.5" version="0.5">
<species name="s1" taxonId="1"><database><genes><gene id="1" protId="a"/>
<gene id="3" protId="c"/></genes></database></species>
<species name="s2" taxonId="2"><database><genes><gene id="2" protId="b"/></genes></database></species>
<taxonomy><taxon id="0" name="N0"><taxon id="1" name="s1"/><taxon id="2" name="s2"/></taxon></taxonomy>
<groups><orthologGroup id="HOG:1_0" taxonId="0"><score value="0.5"/>
<geneRef id="1"/><paralogGroup><orthologGroup id="HOG:1.1_2" taxonId="2"><geneRef id="2"/>
</orthologGroup></paralogGroup></orthologGroup></groups></orthoXML>'''
OWNERS = {"a": "s1", "b": "s2", "c": "s1", "d": "s2"}


def parse(tmp_path, text=XML):
    path = tmp_path / "data.xml"
    path.write_text(text)
    return read_xml(path, OWNERS)


def test_partial_coverage_reported_without_inventing_assignments(tmp_path):
    content, members = parse(tmp_path)
    assert members == {"a": "HOG:1", "b": "HOG:1"}
    assert content["input_proteins"] == 4
    assert content["declared_proteins"] == 3
    assert content["input_not_declared"] == 1
    assert content["declared_not_hog_referenced"] == 1
    assert content["hog_referenced_proteins"] == 2
    assert content["root_hogs"] == 1


@pytest.mark.parametrize("before,after", [
    ("FastOMA 0.3.5", "FastOMA 0.4.0"), ("http://orthoXML.org/2011/", "wrong"),
    ('protId="a"', 'protId="unknown"'), ('protId="b"', 'protId="a"'),
    ('gene id="2"', 'gene id="1"'), ('name="s2" taxonId="2"', 'name="s1" taxonId="2"'),
    ('geneRef id="1"', 'geneRef id="99"'), ('geneRef id="2"', 'geneRef id="1"'),
    ('taxonId="0"', 'taxonId="999"'), ('value="0.5"', 'value="nan"'),
    ('value="0.5"', 'value="inf"'), ('taxon id="2" name="s2"', 'taxon id="2" name="wrong"'),
    ('taxon id="2"', 'taxon id="1"'), ('HOG:1.1_2', 'HOG:1_0'),
    ('<geneRef id="1"/>', '<unknown/>'),
    ('</groups>', '<geneRef id="3"/></groups>'),
])
def test_invalid_xml_rejected(tmp_path, before, after):
    with pytest.raises(ValueError):
        parse(tmp_path, XML.replace(before, after))


def test_truncation_rejected(tmp_path):
    with pytest.raises(ET.ParseError):
        parse(tmp_path, XML[:-10])


@pytest.mark.parametrize("group", ['<orthologGroup id="HOG:2_0" taxonId="0"/>',
                                   '<orthologGroup id="HOG:1_1" taxonId="1"/>'])
def test_empty_or_ambiguous_root_group_rejected(tmp_path, group):
    with pytest.raises(ValueError):
        parse(tmp_path, XML.replace("</groups>", group + "</groups>"))


def test_wrong_root_table_header(tmp_path):
    path = tmp_path / "roots.tsv"
    path.write_text("Group\tProtein\n")
    with pytest.raises(ValueError, match="header"):
        check_root_table(path, {"a": "HOG:1"})


def test_empty_native_pairs_rejected(tmp_path):
    path = tmp_path / "pairs.gz"
    path.write_bytes(gzip.compress(b""))
    with pytest.raises(ValueError, match="Empty"):
        check_pair_scope(path, OWNERS, {"a": "HOG:1"})


def test_root_table_exact_membership(tmp_path):
    _, members = parse(tmp_path)
    path = tmp_path / "roots.tsv"
    path.write_text("RootHOG\tProtein\tOMAmerRootHOG\nHOG:1\tb\torigin\nHOG:1\ta\torigin\n")
    assert check_root_table(path, members) == 2


@pytest.mark.parametrize("rows", ["HOG:1\ta\tx\n", "HOG:2\ta\tx\nHOG:1\tb\tx\n",
                                  "HOG:1\ta\tx\nHOG:1\ta\tx\n", "HOG:1\tforeign\tx\n",
                                  "HOG:1\ta\n", "HOG:1\ta\t\n"])
def test_invalid_root_table(tmp_path, rows):
    path = tmp_path / "roots.tsv"
    path.write_text("RootHOG\tProtein\tOMAmerRootHOG\n" + rows)
    with pytest.raises(ValueError):
        check_root_table(path, {"a": "HOG:1", "b": "HOG:1"})


def test_pair_scope_and_absent_endpoints(tmp_path):
    path = tmp_path / "pairs.gz"
    path.write_bytes(gzip.compress(b"a\tb\n"))
    result = check_pair_scope(path, OWNERS, {"a": "r1", "b": "r1", "c": "r2"})
    assert result == {"native_pair_rows": 1, "pair_endpoint_proteins": 2, "hog_members_without_pairs": 1}


@pytest.mark.parametrize("members", [{"a": "r1", "b": "r2"}, {"a": "r1"}])
def test_pair_scope_rejects_cross_root_or_unreferenced_gene(tmp_path, members):
    path = tmp_path / "pairs.gz"
    path.write_bytes(gzip.compress(b"a\tb\n"))
    with pytest.raises(ValueError):
        check_pair_scope(path, OWNERS, members)
