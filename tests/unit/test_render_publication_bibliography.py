import copy
import hashlib
import json
from pathlib import Path
import shutil

import pytest

from benchmark_tools.render_publication_bibliography import audit, entries, render, text


ROOT = Path(__file__).resolve().parents[2]


def document(value, identity="ref-a"):
    return {"blocks": [{"t": "Div", "c": [[identity, ["csl-entry"], []],
        [{"t": "Para", "c": [{"t": "Str", "c": value}]}]]}]}


def test_visible_text_excludes_link_target_and_attributes():
    node = {"t": "Link", "c": [["secret", [], []],
        [{"t": "Emph", "c": [{"t": "Str", "c": "Visible"}]},
         {"t": "Space"}, {"t": "Str", "c": "label"}],
        ["https://hidden.example", "hidden title"]]}
    assert text(node) == "Visible label"
    assert entries(document("Cafe\u0301\n title")) == {"a": "Caf\u00e9 title"}


@pytest.mark.parametrize("problem", ["duplicate_render", "wrong_prefix", "missing", "extra", "duplicate_input"])
def test_inventory_mismatch_rejected(problem):
    records = [{"id": "a", "title": "Title"}]
    doc = document("Title")
    if problem == "duplicate_render":
        doc["blocks"] *= 2
    elif problem == "wrong_prefix":
        doc = document("Title", "a")
    elif problem == "missing":
        doc["blocks"] = []
    elif problem == "extra":
        doc["blocks"].extend(document("Other", "ref-b")["blocks"])
    else:
        records *= 2
    with pytest.raises(ValueError):
        audit(records, doc)


@pytest.mark.parametrize("field,fragment", [
    ("title", "Unique title"), ("genre", "preprint"),
    ("DOI", "10.123/example"), ("URL", "https://example.org"),
    ("container-title", "Journal name"), ("number", "Identifier: e123"),
    ("author_0_family", "Buida"), ("author_0_suffix", "III"),
    ("author_1_literal", "Test Consortium"), ("issued_year", "2024"),
    ("accessed", "Accessed: 2026-09-19"),
])
def test_missing_checked_field_is_reported(field, fragment):
    record = {"id": "a", "title": "Unique title", "genre": "preprint",
        "DOI": "10.123/example", "URL": "https://example.org",
        "container-title": "Journal name", "number": "e123",
        "author": [{"family": "Buida", "suffix": "III"}, {"literal": "Test Consortium"}],
        "issued": {"date-parts": [[2024]]}, "accessed": {"date-parts": [[2026, 9, 19]]}}
    value = "Unique title preprint 10.123/example https://example.org Journal name " \
        "Identifier: e123 Buida III Test Consortium 2024 Accessed: 2026-09-19"
    assert audit([record], document(value))["all_checked_fields_present"]
    result = audit([record], document(value.replace(fragment, "")))
    assert not result["all_checked_fields_present"]
    assert result["checks"][0]["missing"] == {field: fragment}


def test_existing_destination_and_missing_executable(tmp_path):
    with pytest.raises(FileExistsError):
        render(tmp_path / "unused", tmp_path / "unused", tmp_path)
    destination = tmp_path / "new"
    with pytest.raises(FileNotFoundError):
        render(tmp_path / "unused", tmp_path / "unused", destination,
               executable="nonexistent-orthohmm-citation-processor")
    assert not destination.exists()


@pytest.mark.skipif(shutil.which("pandoc") is None, reason="Pandoc required")
def test_actual_selected_bibliography_renders_after_relocation(tmp_path, monkeypatch):
    source = ROOT / "benchmark_tools/results/publication_bibliography_20260919_v3.csl.json"
    bibliography = tmp_path / "references.json"
    style = tmp_path / "review.csl"
    shutil.copyfile(source, bibliography)
    shutil.copyfile(ROOT / "benchmark_tools/publication-review.csl", style)
    before = bibliography.read_bytes()
    monkeypatch.chdir(tmp_path)
    result = render(bibliography, style, tmp_path / "rendered")
    assert result["field_audit"]["entries"] == 37
    assert result["field_audit"]["all_checked_fields_present"]
    assert result["stderr"] == ["", ""]
    assert not result["publication_ready"]
    assert bibliography.read_bytes() == before
    for item in result["outputs"]:
        assert hashlib.sha256(Path(item["path"]).read_bytes()).hexdigest() == item["sha256"]
    manifest = json.loads((tmp_path / "rendered/manifest.json").read_text())
    assert manifest == result
    rendered = {r["id"]: r["text"] for r in result["field_audit"]["checks"]}
    assert "III" in rendered["orthohmm2024preprint"]
    assert "preprint" in rendered["orthohmm2024preprint"]
    assert "Identifier: 238" in rendered["orthofinder2019"]
    assert "1695" in rendered["Csardi2006igraph"]
    assert "10." not in rendered["Csardi2006igraph"]
    doc = json.loads((tmp_path / "rendered/bibliography.pandoc.json").read_text())
    broken = copy.deepcopy(doc)
    def drop_suffix(node):
        if isinstance(node, dict):
            if node.get("t") == "Str" and "III" in node.get("c", ""):
                node["c"] = node["c"].replace("III", "")
            for value in node.values():
                drop_suffix(value)
        elif isinstance(node, list):
            for value in node:
                drop_suffix(value)
    drop_suffix(broken)
    assert not audit(json.loads(before), broken)["all_checked_fields_present"]
