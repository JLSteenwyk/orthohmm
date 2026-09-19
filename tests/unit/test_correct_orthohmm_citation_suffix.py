from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools import correct_orthohmm_citation_suffix as module


def records():
    root = Path(__file__).resolve().parents[2] / "benchmark_tools/results"
    return json.loads((root / "publication_references_20260918.csl.json").read_text())


def test_only_one_suffix_changes():
    original = records()
    before = deepcopy(original)
    result = module.corrected(original)
    assert original == before
    row = next(r for r in result if r["id"] == "orthohmm2024preprint")
    assert row["author"][1].pop("suffix") == "III"
    assert result == original


@pytest.mark.parametrize("fault", ["absent", "duplicate", "doi", "order", "type", "already_corrected"])
def test_wrong_citation_rejected(fault):
    data = records()
    row = next(r for r in data if r["id"] == "orthohmm2024preprint")
    if fault == "absent": data.remove(row)
    elif fault == "duplicate": data.append(deepcopy(row))
    elif fault == "doi": row["DOI"] = "other"
    elif fault == "order": row["author"].reverse()
    elif fault == "type": row["genre"] = "journal article"
    else: row["author"][1]["suffix"] = "III"
    with pytest.raises(ValueError): module.corrected(data)


def test_export_provenance_and_altered_source_rejected(tmp_path, monkeypatch):
    raw, author, api = [tmp_path / n for n in ("raw.json", "author.html", "api.json")]
    raw.write_text(json.dumps(records()))
    author.write_text("reviewed author source fixture")
    api.write_text("{}")
    for constant, path in (("RAW_SHA", raw), ("AUTHOR_SHA", author), ("API_SHA", api)):
        monkeypatch.setattr(module, constant, module.record(path)["sha256"])
    output, provenance = tmp_path / "output.json", tmp_path / "provenance.json"
    result = module.export(raw, author, api, output, provenance)
    assert result["change"]["after"]["suffix"] == "III"
    assert json.loads(output.read_text()) == module.corrected(records())
    with pytest.raises(FileExistsError): module.export(raw, author, api, output, provenance)
    author.write_text("altered")
    with pytest.raises(ValueError, match="bytes changed"):
        module.export(raw, author, api, tmp_path / "new.json", tmp_path / "new-provenance.json")
    assert not (tmp_path / "new.json").exists()
