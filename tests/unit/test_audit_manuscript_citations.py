import json

import pytest

from benchmark_tools.audit_manuscript_citations import audit, inventory


def link(url):
    return {"t": "Link", "c": [["", [], []], [{"t": "Str", "c": "Reference"}], [url, ""]]}


def document(*nodes):
    return {"blocks": [{"t": "Para", "c": list(nodes)}]}


@pytest.mark.parametrize("url", ["https://doi.org/10.1000/ABC", "http://dx.doi.org/10.1000/abc", "https://doi.org/10.1000%2Fabc"])
def test_resolver_variants_match_same_doi(url):
    result = inventory(document(link(url)), [{"id": "a", "DOI": "10.1000/abc"}])
    assert result["matched_bibliography_ids"] == ["a"]
    assert result["unresolved_citations"] == []
    assert result["citation_adequacy_established"] is False


def test_unknown_doi_not_conflated_with_other_external_or_local_links():
    result = inventory(document(link("https://doi.org/10.1000/missing"),
        link("https://doi.org.example/10.1000/abc"), link("evidence.md")), [{"id": "a", "DOI": "10.1000/abc"}])
    assert result["unresolved_citations"][0]["status"] == "unmatched_doi"
    assert len(result["other_external_links"]) == 1
    assert result["internal_link_occurrences"] == 1
    assert result["bibliography_not_explicitly_matched"] == ["a"]


def test_exact_url_and_ambiguous_url_are_distinct():
    rows = [{"id": "a", "URL": "https://example.org"}]
    assert inventory(document(link("https://example.org")), rows)["matched_bibliography_ids"] == ["a"]
    rows.append({"id": "b", "URL": "https://example.org"})
    result = inventory(document(link("https://example.org")), rows)
    assert result["unresolved_citations"][0]["status"] == "ambiguous"
    assert result["matched_bibliography_ids"] == []


@pytest.mark.parametrize("rows", [[{"id": "a"}, {"id": "a"}],
    [{"id": "a", "DOI": "10.1/A"}, {"id": "b", "DOI": "10.1/a"}]])
def test_ambiguous_bibliography_rejected(rows):
    with pytest.raises(ValueError):
        inventory(document(), rows)


def test_citation_ids_and_section_context():
    parsed = document({"t": "Cite", "c": [[{"citationId": "a"}, {"citationId": "missing"}], []]})
    parsed["blocks"].insert(0, {"t": "Header", "c": [2, ["", [], []], [{"t": "Str", "c": "Methods"}]]})
    result = inventory(parsed, [{"id": "a"}])
    assert result["matched_bibliography_ids"] == ["a"]
    assert result["unresolved_citations"][0]["section"] == "Methods"
    assert result["unresolved_citations"][0]["status"] == "unknown_citation_id"


def test_real_parser_ignores_code_and_keeps_reference_style_links(tmp_path):
    manuscript = tmp_path / "draft.md"
    manuscript.write_text('## Methods\n\n[Reference][r] and [@b].\n\n[r]: https://doi.org/10.1000/A\n\n```text\n[not a citation](https://doi.org/10.1000/fake)\n```\n')
    bibliography = tmp_path / "refs.json"
    bibliography.write_text(json.dumps([{"id": "a", "DOI": "10.1000/a"}, {"id": "b"}]))
    result = audit(manuscript, bibliography)
    assert result["matched_bibliography_ids"] == ["a", "b"]
    assert len(result["occurrences"]) == 2
    assert result["unresolved_citations"] == []
