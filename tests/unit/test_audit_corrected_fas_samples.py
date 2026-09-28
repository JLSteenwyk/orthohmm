import gzip
import json

import pytest

from benchmark_tools import audit_corrected_fas_samples as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))
    return record(path)


@pytest.fixture
def fixture(tmp_path, monkeypatch):
    methods = []
    for i in range(8):
        base = tmp_path / str(i)
        endpoint = base / "results/FAS/FAS.json"
        participant = dict(participant_id=str(i), metric_x=10, metric_y=0.5, stderr_y=0.25)
        endpoint_ref = write(endpoint, dict(datalink=dict(inline_data=dict(
            visualization=dict(x_axis="NR_ORTHOLOGS", y_axis="FAS"),
            challenge_participants=[participant]))))
        raw = endpoint.parent / "sample_raw.txt.gz"
        with gzip.open(raw, "wt") as stream:
            stream.write("Acc1\tAcc2\tFAS\nA\tB\t0.25\nA\tC\t0.75\n")
        execution = write(base / "execution.json", dict(outputs=[record(raw)]))
        admission = write(base / "admission.json", dict(execution_report=execution, metric_files=[endpoint_ref]))
        methods.append(dict(key=str(i), participant=str(i), status="admitted", admission=admission,
                            scores=dict(FAS=0.5), details=dict(FAS=dict(assessed_relations=10))))
    manifest = tmp_path / "manifest.json"
    ref = write(manifest, dict(methods=methods))
    monkeypatch.setattr(module, "MANIFEST", "manifest.json")
    monkeypatch.setattr(module, "MANIFEST_SHA", ref["sha256"])
    return tmp_path, methods


def test_complete_corrected_panel(fixture):
    root, _ = fixture
    result = module.audit(root)
    assert len(result["methods"]) == 8
    assert not result["uncertainty_admitted"]
    for row in result["methods"]:
        assert row["mean"] == 0.5
        assert row["sample_fraction"] == 0.2
        assert row["proteins_in_multiple_sample_pairs"] == 1
        assert row["historical_execution_output_bound"]


def test_comparator_participant_from_admission(fixture, monkeypatch):
    root, methods = fixture
    methods[0].pop("participant")
    path = root / "0/admission.json"
    admission = json.loads(path.read_text())
    admission["assessment"] = dict(participant="0")
    methods[0]["admission"] = write(path, admission)
    ref = write(root / "manifest.json", dict(methods=methods))
    monkeypatch.setattr(module, "MANIFEST_SHA", ref["sha256"])
    assert len(module.audit(root)["methods"]) == 8


@pytest.mark.parametrize("problem", ["manifest", "raw", "missing_raw_pin", "endpoint", "mean", "count", "sem", "axes", "participant", "duplicate_method"])
def test_rejects_inconsistent_provenance(fixture, monkeypatch, problem):
    root, methods = fixture
    admission_path = root / "0/admission.json"
    admission = json.loads(admission_path.read_text())
    if problem == "manifest":
        monkeypatch.setattr(module, "MANIFEST_SHA", "0" * 64)
    elif problem == "raw":
        with gzip.open(root / "0/results/FAS/sample_raw.txt.gz", "wt") as stream:
            stream.write("Acc1\tAcc2\tFAS\nA\tB\t0.5\nA\tC\t0.5\n")
    elif problem == "missing_raw_pin":
        admission["execution_report"] = write(root / "0/execution.json", dict(outputs=[]))
    elif problem == "endpoint":
        admission["metric_files"] *= 2
    elif problem == "duplicate_method":
        methods[1]["key"] = methods[0]["key"]
    else:
        path = root / "0/results/FAS/FAS.json"
        data = json.loads(path.read_text())
        inline = data["datalink"]["inline_data"]
        if problem == "axes":
            inline["visualization"]["y_axis"] = "GO"
        elif problem == "participant":
            inline["challenge_participants"][0]["participant_id"] = "another_method"
        else:
            inline["challenge_participants"][0][{"mean": "metric_y", "count": "metric_x", "sem": "stderr_y"}[problem]] = 0.1
        admission["metric_files"] = [write(path, data)]
    if problem != "manifest":
        methods[0]["admission"] = write(admission_path, admission)
        manifest = write(root / "manifest.json", dict(methods=methods))
        monkeypatch.setattr(module, "MANIFEST_SHA", manifest["sha256"])
    with pytest.raises(ValueError):
        module.audit(root)
