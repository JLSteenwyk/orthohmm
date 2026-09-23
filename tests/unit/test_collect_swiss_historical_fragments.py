import json

import pytest

from benchmark_tools import collect_swiss_historical_fragments as panel
from benchmark_tools.prepare_ob_candidate_neighborhood import record


@pytest.mark.parametrize("problem", [None, "missing", "duplicate", "description", "version", "source"])
def test_extraction(tmp_path, problem):
    path = tmp_path / "input.fa"
    description = "sp|P12345|TEST Protein OX=9606 SV=1"
    if problem == "version":
        description = description.replace(" SV=1", "")
    text = f">{description}\nAAA\n"
    path.write_text(text*2 if problem == "duplicate" else "" if problem == "missing" else text)
    descriptors = {"P12345": dict(description=description, length=3, input_id="sp|P12345|TEST", source_file="input.fa")}
    if problem == "description":
        descriptors["P12345"]["description"] += " changed"
    if problem == "source":
        descriptors["P12345"]["source_file"] = "other.fa"
    if problem:
        with pytest.raises(ValueError):
            panel.extract({"family": ["P12345"]}, descriptors, [record(path)])
    else:
        result = panel.extract({"family": ["P12345"]}, descriptors, [record(path)])
        assert result["P12345"]["taxid"] == "9606"
        assert result["P12345"]["sequence_version"] == 1


def test_collection_preserves_failure_and_no_overwrite(tmp_path, monkeypatch):
    genes = {g: dict(sequence_version=1, sequence_sha256="a"*64, taxid="9606") for g in ["A", "B"]}
    calls, delays = [], []
    def acquire(accession, *args):
        calls.append(accession)
        path = args[-1]
        path.mkdir()
        if accession == "A":
            raise ValueError("Sequence mismatch")
        (path / "audit.json").write_text("{}")
        return dict(selection_class="baseline_release")
    monkeypatch.setattr(panel, "acquire", acquire)
    monkeypatch.setattr(panel.time, "sleep", delays.append)
    output = tmp_path / "collection"
    result = panel.collect(dict(genes=genes, records=[]), output)
    assert calls == ["A", "B"] and delays == [1.0, 1.0]
    assert result["matched"] == result["missing"] == 1
    assert result["entries"]["A"]["error"] == "Sequence mismatch"
    assert result["annotation_panel_admitted"] is False
    assert json.loads((output / "status.json").read_text()) == result
    with pytest.raises(FileExistsError):
        panel.collect(dict(genes=genes, records=[]), output)
    with pytest.raises(ValueError):
        panel.collect(dict(genes=genes, records=[]), tmp_path / "bad", delay=0)
