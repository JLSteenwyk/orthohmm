from copy import deepcopy
import json

import pytest

from benchmark_tools.export_publication_citations import citation, export, record


def fixture():
    return ({"DOI": "10.1234/test", "type": "journal-article", "title": "Study of <i>genes</i> &amp; trees",
             "author": [{"given": "A.", "family": "Author", "suffix": "III"}, {"name": "Consortium"}],
             "issued": {"date-parts": [[2020, 2, 1]]}, "container-title": "Journal", "volume": "2",
             "abstract": "Do not export", "reference": [{"DOI": "unrelated"}]},
            {"id": "test", "doi": "10.1234/test", "year": 2020})


def test_valid_citation_and_no_abstract():
    metadata, item = fixture()
    result = citation(metadata, item)
    assert result["type"] == "article-journal"
    assert result["title"] == "Study of genes & trees"
    assert result["author"][0]["suffix"] == "III"
    assert result["author"][1] == {"literal": "Consortium"}
    assert "abstract" not in result and "reference" not in result


def test_preprint_not_journal():
    metadata, item = fixture()
    metadata.update(type="posted-content", **{"container-title": []})
    result = citation(metadata, item)
    assert result["genre"] == "preprint" and result["type"] == "article"
    assert "container-title" not in result


@pytest.mark.parametrize("key,value", [("DOI", "10.9999/wrong"), ("type", "unsupported"), ("title", []),
                                      ("title", ""), ("author", []), ("author", [{"given": "Only"}]),
                                      ("issued", {"date-parts": [[2019]]}), ("volume", ["2"])])
def test_invalid_metadata(key, value):
    metadata, item = fixture()
    metadata[key] = deepcopy(value)
    with pytest.raises(ValueError):
        citation(metadata, item)


def test_article_identifier_is_not_pagination():
    metadata, item = fixture()
    metadata.update({"article-number": "eaaz5667", "page": "1-9"})
    result = citation(metadata, item)
    assert result["number"] == "eaaz5667"
    assert result["page"] == "1-9"
    metadata["article-number"] = ["invalid"]
    with pytest.raises(ValueError):
        citation(metadata, item)


@pytest.mark.parametrize("mutation", [None, "selection", "response", "inventory", "url"])
def test_offline_replay(tmp_path, monkeypatch, mutation):
    metadata, item = fixture()
    manifest = tmp_path / "selection.json"
    manifest.write_text(json.dumps([item]))
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir()
    raw = raw_dir / "test.json"
    raw.write_text(json.dumps(metadata))
    source = {"id": "test", "doi": item["doi"], "response": record(raw),
              "url": "https://api.crossref.org/works/10.1234%2Ftest/transform/application/vnd.citationstyles.csl+json",
              "retrieved_utc": "2026-09-18T00:00:00+00:00"}
    cached = {"selection": record(manifest), "records": [source]}
    if mutation == "selection":
        cached["selection"]["sha256"] = "wrong"
    elif mutation == "response":
        raw.write_text("{}")
    elif mutation == "inventory":
        source["id"] = "wrong"
    elif mutation == "url":
        source["url"] = "https://example.org"
    previous = tmp_path / "previous.json"
    previous.write_text(json.dumps(cached))

    def no_network(*args, **kwargs):
        pytest.fail("Offline replay attempted network access")

    monkeypatch.setattr("benchmark_tools.export_publication_citations.urlopen", no_network)
    output, provenance = tmp_path / "out.json", tmp_path / "provenance.json"
    if mutation:
        with pytest.raises(ValueError):
            export(manifest, raw_dir, output, provenance, previous)
        assert not output.exists() and not provenance.exists()
    else:
        result = export(manifest, raw_dir, output, provenance, previous)
        assert json.loads(output.read_text()) == [citation(metadata, item)]
        assert result["records"][0]["retrieved_utc"] == source["retrieved_utc"]
        with pytest.raises(FileExistsError):
            export(manifest, raw_dir, output, provenance, previous)
