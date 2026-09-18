import pytest

from benchmark_tools.correct_qfo_service_bylines import byline


def article(contributors, doi="test", extra=""):
    return (f'<article><front><article-meta><article-id pub-id-type="doi">{doi}</article-id>'
            f'<contrib-group>{contributors}</contrib-group></article-meta></front>{extra}</article>')


PERSON = '<contrib><name><surname>Smith</surname><given-names>A</given-names></name></contrib>'


def test_order_repetition_and_collective_preserved():
    authors = byline(article(PERSON + '<contrib><collab>the A team the B Consortium</collab></contrib>' + PERSON), "test")
    assert authors == [{"given": "A", "family": "Smith"},
                       {"literal": "the A team the B Consortium"},
                       {"given": "A", "family": "Smith"}]


def test_consortium_members_outside_byline_excluded():
    assert len(byline(article(PERSON, extra='<back><contrib-group>' + PERSON + '</contrib-group></back>'), "test")) == 1


def test_2020_literal_conversion():
    name = '<contrib><name><surname>for&#160;Orthologs&#160;Consortium</surname><given-names>the Quest</given-names></name></contrib>'
    assert byline(article(name, doi="10.1093/nar/gkaa308"), "10.1093/nar/gkaa308") == [
        {"literal": "the Quest for Orthologs Consortium"}]


@pytest.mark.parametrize("contributor", ["", '<contrib/>', '<contrib><collab/></contrib>',
    '<contrib><name><surname>X</surname></name></contrib>',
    '<contrib><collab>X<contrib/></collab></contrib>',
    '<contrib contrib-type="editor"><collab>X</collab></contrib>'])
def test_malformed_byline_rejected(contributor):
    with pytest.raises(ValueError):
        byline(article(contributor), "test")


def test_wrong_doi_rejected():
    with pytest.raises(ValueError, match="DOI"):
        byline(article(PERSON), "wrong")
