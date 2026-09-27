import gzip
import json

import pytest

from benchmark_tools import run_qfo_scored_pair_panel as panel
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def fixture(tmp_path, monkeypatch, defect=None):
    methods = []
    for i in range(8):
        records, outputs = [], []
        for metric in ("GO", "EC"):
            directory = tmp_path/str(i)/"results"/metric
            directory.mkdir(parents=True)
            endpoint = directory/(metric+".json")
            endpoint.write_text("{}")
            records.append(record(endpoint))
            with gzip.open(directory/"pairs_raw.txt.gz", "wt") as stream:
                stream.write(f"# {metric} Similarities between orthologs from test\n"
                    f"# Computing timestamp: test\n"
                    f"# Protein ID 1<tab>Protein ID 2<tab>{metric} Similarity\n"
                    "A\tB\t0.500000\n")
            outputs.append(record(directory/"pairs_raw.txt.gz"))
        execution = tmp_path/str(i)/"execution.json"
        execution.write_text(json.dumps(dict(outputs=outputs)))
        admission = tmp_path/str(i)/"admission.json"
        admission.write_text(json.dumps(dict(metric_files=records, execution_report=record(execution))))
        methods.append(dict(key=str(i), label=str(i), prediction_semantics="fixture",
            status="admitted", admission=record(admission), scores=dict(GO=.5, EC=.5),
            details={m:dict(assessed_relations=1) for m in ("GO", "EC")}))
    if defect == "mean":
        methods[0]["scores"]["GO"] = .6
    if defect == "count":
        methods[0]["details"]["GO"]["assessed_relations"] = 2
    if defect == "status":
        methods[0]["status"] = "pending"
    if defect == "duplicate":
        methods[0]["key"] = methods[1]["key"]
    manifest = tmp_path/"manifest.json"
    manifest.write_text(json.dumps(dict(methods=methods)))
    monkeypatch.setattr(panel, "MANIFEST", "manifest.json")
    monkeypatch.setattr(panel, "MANIFEST_SHA", record(manifest)["sha256"])
    return manifest


def test_complete_panel(tmp_path, monkeypatch):
    fixture(tmp_path, monkeypatch)
    output = tmp_path/"result.json"
    result = panel.run(tmp_path, output)
    assert len(result["comparisons"]) == 56
    assert len(result["endpoints"]) == 16
    assert all(r["result"]["original_mean_difference"] == 0 for r in result["comparisons"])
    assert result["uncertainty_admitted"] is False
    with pytest.raises(FileExistsError):
        panel.run(tmp_path, output)


@pytest.mark.parametrize("defect", ["mean", "count", "status", "duplicate"])
def test_bad_panel_rejected_without_report(tmp_path, monkeypatch, defect):
    fixture(tmp_path, monkeypatch, defect)
    output = tmp_path/"result.json"
    with pytest.raises(ValueError):
        panel.run(tmp_path, output)
    assert not output.exists()


def test_changed_manifest(tmp_path, monkeypatch):
    manifest = fixture(tmp_path, monkeypatch)
    manifest.write_text(manifest.read_text()+"\n")
    with pytest.raises(ValueError, match="manifest"):
        panel.run(tmp_path, tmp_path/"result.json")


def test_same_count_and_mean_changed_pairs_rejected(tmp_path, monkeypatch):
    fixture(tmp_path, monkeypatch)
    raw = tmp_path/"0/results/GO/pairs_raw.txt.gz"
    with gzip.open(raw, "rt") as stream:
        text = stream.read()
    with gzip.open(raw, "wt") as stream:
        stream.write(text.replace("A\tB\t", "C\tD\t"))
    with pytest.raises(ValueError, match="historical execution"):
        panel.run(tmp_path, tmp_path/"result.json")


def test_checked_record_execution_route(tmp_path):
    directory = tmp_path/"qfo_corrected_assessment_v1/sonic"
    directory.mkdir(parents=True)
    path = directory/"results.json"
    path.write_text('{"outputs": []}')
    identity, data = panel.execution_binding(dict(checked_records=[record(path)]))
    assert identity == record(path)
    assert data == dict(outputs=[])
    with pytest.raises(ValueError, match="Ambiguous"):
        panel.execution_binding(dict(checked_records=[]))
