import hashlib
import json
from pathlib import Path
import shutil

import pytest

from benchmark_tools.assemble_publication_citations import assemble


def setup_inputs(tmp_path, rows):
    source = tmp_path / "source.json"
    source.write_text(json.dumps(rows))
    selection = [{"file": source.name, "sha256": hashlib.sha256(source.read_bytes()).hexdigest()}]
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps(selection))
    return manifest, tmp_path / "out.json", tmp_path / "provenance.json"


def test_preserves_all_fields(tmp_path):
    rows = [{"id": "a", "type": "article", "title": "A", "genre": "preprint",
             "author": [{"literal": "Consortium"}], "number": "e123", "page": "1-5"},
            {"id": "web", "type": "webpage", "title": "Resource",
             "accessed": {"date-parts": [[2026, 9, 19]]}}]
    args = setup_inputs(tmp_path, rows)
    report = assemble(*args)
    assert json.loads(args[1].read_text()) == rows
    assert report["count"] == 2
    assert [r["id"] for r in report["records"]] == ["a", "web"]
    with pytest.raises(FileExistsError):
        assemble(*args)


@pytest.mark.parametrize("rows", [[], {}, [None], [{"id": "a"}],
    [{"id": "a", "type": "article", "title": " "}],
    [{"id": "a", "type": "article", "title": "A", "DOI": None}],
    [{"id": "a", "type": "article", "title": "A"}] * 2,
    [{"id": "a", "type": "article", "title": "A", "DOI": "10.1/ABC"},
     {"id": "b", "type": "article", "title": "B", "DOI": " 10.1/abc "}],
])
def test_invalid_rows_leave_no_outputs(tmp_path, rows):
    args = setup_inputs(tmp_path, rows)
    with pytest.raises(ValueError):
        assemble(*args)
    assert not args[1].exists() and not args[2].exists()


@pytest.mark.parametrize("change", ["hash", "duplicate", "escape", "empty", "same_output"])
def test_manifest_and_destinations(tmp_path, change):
    manifest, output, provenance = setup_inputs(tmp_path, [{"id": "a", "type": "article", "title": "A"}])
    selection = json.loads(manifest.read_text())
    if change == "hash":
        selection[0]["sha256"] = "wrong"
    elif change == "duplicate":
        selection *= 2
    elif change == "escape":
        selection[0]["file"] = "../source.json"
    elif change == "empty":
        selection = []
    else:
        provenance = output
    manifest.write_text(json.dumps(selection))
    with pytest.raises(ValueError):
        assemble(manifest, output, provenance)
    assert not output.exists() and not provenance.exists()


def test_relocated_real_selection(tmp_path, monkeypatch):
    root = Path(__file__).resolve().parents[2] / "benchmark_tools/results"
    name = "publication_bibliography_selection_20260919.json"
    shutil.copyfile(root / name, tmp_path / name)
    expected = []
    for source in json.loads((root / name).read_text()):
        shutil.copyfile(root / source["file"], tmp_path / source["file"])
        expected.extend(json.loads((root / source["file"]).read_text()))
    monkeypatch.chdir(tmp_path)
    output = tmp_path / "combined.json"
    result = assemble(tmp_path / name, output, tmp_path / "report.json")
    actual = json.loads(output.read_text())
    assert actual == expected and result["count"] == 36
    by_id = {r["id"]: r for r in actual}
    assert len(by_id["Li2006TreeFam"]["author"]) == 15
    assert len(by_id["qfo2022"]["author"]) == 31
    assert by_id["orthohmm2024preprint"]["genre"] == "preprint"
    assert "issued" not in by_id["SIBSwissTree"]
    assert "orthofinder2026correction" in by_id


def test_v2_adds_only_reviewed_igraph_and_relocates(tmp_path, monkeypatch):
    root = Path(__file__).resolve().parents[2] / "benchmark_tools/results"
    name = "publication_bibliography_selection_20260919_v2.json"
    selection = json.loads((root / name).read_text())
    original = json.loads((root / "publication_bibliography_selection_20260919.json").read_text())
    assert selection[:-1] == original
    shutil.copyfile(root / name, tmp_path / name)
    for source in selection:
        shutil.copyfile(root / source["file"], tmp_path / source["file"])
    monkeypatch.chdir(tmp_path)
    output = tmp_path / "bibliography.json"
    result = assemble(tmp_path / name, output, tmp_path / "provenance.json")
    actual = json.loads(output.read_text())
    assert result["count"] == 37
    assert actual[:-1] == json.loads((root / "publication_bibliography_20260919.csl.json").read_text())
    item = actual[-1]
    assert item["id"] == "Csardi2006igraph" and "DOI" not in item
    assert item["issued"] == {"date-parts": [[2006]]}
    assert item["volume"] == "Complex Systems" and item["page"] == "1695"
    assert item["author"] == [{"given": "G\u00e1bor", "family": "Cs\u00e1rdi"},
                              {"given": "Tam\u00e1s", "family": "Nepusz"}]
    provenance = json.loads((root / "publication_igraph_reference_provenance_20260919.json").read_text())
    assert provenance["output"]["sha256"] == selection[-1]["sha256"]
    assert provenance["status"] == "manually_transcribed_official_citation_guidance"


def test_v3_relocates_and_changes_only_reviewed_suffix(tmp_path):
    root = Path(__file__).resolve().parents[2] / "benchmark_tools/results"
    name = "publication_bibliography_selection_20260919_v3.json"
    shutil.copyfile(root / name, tmp_path / name)
    for source in json.loads((root / name).read_text()):
        shutil.copyfile(root / source["file"], tmp_path / source["file"])
    output = tmp_path / "out.json"
    result = assemble(tmp_path / name, output, tmp_path / "provenance.json")
    actual = json.loads(output.read_text())
    assert result["count"] == 37
    row = next(r for r in actual if r["id"] == "orthohmm2024preprint")
    assert row["author"][1].pop("suffix") == "III"
    assert actual == json.loads((root / "publication_bibliography_20260919_v2.csl.json").read_text())
