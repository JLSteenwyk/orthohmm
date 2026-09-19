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
