import copy
import json
from pathlib import Path

import pytest

from benchmark_tools import audit_qfo_corrected_factorial_swiss as module
from benchmark_tools import audit_qfo_factorial_swiss, bootstrap_qfo_factorial
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_audit_qfo_factorial_swiss import synthetic
from tests.unit.test_export_qfo_corrected_factorial import fixture


def execution_fixture():
    report = fixture(1)
    report["scheduler"] = {"State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon",
                           "AllocCPUS": "8", "JobIDRaw": "456"}
    report["pairs_manifest"] = {"path": "/pairs", "sha256": "fixture", "bytes": 1}
    execution = {"status": "process_succeeded_pending_independent_admission", "exit_code": 0,
        "job_id": "456", "index": 1, "cell": report["cell"], "pairs_manifest": report["pairs_manifest"],
        "stage": report["conversion"], "outputs": [{"path": "/fixture/SwissTrees/test.raw.txt.gz"}]}
    return report, execution


@pytest.mark.parametrize("problem", [None, "running", "failed", "job", "index", "cell", "pairs", "stage", "missing", "duplicate"])
def test_execution_binding(problem):
    report, execution = execution_fixture()
    if problem == "running":
        report["scheduler"]["State"] = "RUNNING"
    elif problem == "failed":
        execution["exit_code"] = 1
    elif problem == "job":
        execution["job_id"] = "999"
    elif problem in ("index", "cell", "pairs", "stage"):
        key = {"pairs": "pairs_manifest"}.get(problem, problem)
        execution[key] = "wrong"
    elif problem == "missing":
        execution["outputs"] = []
    elif problem == "duplicate":
        execution["outputs"] *= 2
    if problem:
        with pytest.raises(ValueError):
            module.execution_raw(report, execution)
    else:
        assert module.execution_raw(report, execution) == execution["outputs"][0]


def test_complete_corrected_count_audit(tmp_path, monkeypatch):
    entries, baseline = synthetic(tmp_path)
    reference = tmp_path / "reference"
    reference.write_text("synthetic reference\n")
    baseline["reference"] = record(reference)
    for target in (audit_qfo_factorial_swiss, bootstrap_qfo_factorial):
        monkeypatch.setattr(target, "REFERENCE_SHA", baseline["reference"]["sha256"])

    def save(name, value):
        path = tmp_path / name
        path.write_text(json.dumps(value))
        return record(path)

    inventory = {"cells": []}
    for index, entry in enumerate(entries):
        report = fixture(index)
        participant = report["assessment"]["participant"]
        native = copy.deepcopy(entry["assessment"])
        native["participant"] = participant
        for metric in native["native_assessments"]:
            metric["participant_id"] = participant
        report["assessment"].update(native)
        metrics = {m["metrics"]["metric_id"]: m["metrics"]["value"] for m in native["native_assessments"]
                   if m["challenge_id"] == "SwissTrees"}
        endpoint = report["assessment"]["endpoints"]["SwissTrees"]
        endpoint["native_participant"].update(metric_x=metrics["TPR"], metric_y=metrics["PPV"])
        endpoint["score"] = 2 * metrics["TPR"] * metrics["PPV"] / (metrics["TPR"] + metrics["PPV"])
        report["assessment"]["secondary_six_metric_mean"] = (2.5 + endpoint["score"]) / 6
        report["pairs_manifest"] = save(f"pairs_{index}.json", report["conversion"])
        report["scheduler"] = {"State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon",
                               "AllocCPUS": "8", "JobIDRaw": str(500 + index)}
        raw = tmp_path / str(index) / "SwissTrees/test.raw.txt.gz"
        raw.parent.mkdir(parents=True)
        raw.write_bytes(Path(entry["raw_file"]["path"]).read_bytes())
        execution = {"status": "process_succeeded_pending_independent_admission", "exit_code": 0,
            "job_id": str(500 + index), "index": index, "cell": report["cell"],
            "pairs_manifest": report["pairs_manifest"], "stage": report["conversion"], "outputs": [record(raw)]}
        report["execution_report"] = save(f"execution_{index}.json", execution)
        report["checked_records"] = [report["execution_report"]]
        inventory["cells"].append({"cell": report["cell"], "admission": save(f"admission_{index}.json", report)})
    base_record = save("baseline.json", baseline)
    monkeypatch.setattr(module, "BASE_COUNTS_SHA", base_record["sha256"])
    inv = save("inventory.json", inventory)
    result = module.audit(Path(inv["path"]), inv["sha256"], Path(base_record["path"]))
    assert result["status"] == "corrected_qfo_factorial_swiss_counts_verified"
    assert result["reference_relation_count"] == 270
    assert result["uncertainty_admitted"] is False
    assert result["cells"][0]["families"][0]["counts_without_prior"] != result["cells"][2]["families"][0]["counts_without_prior"]
    with pytest.raises(ValueError):
        bootstrap_qfo_factorial.validated_values(result)
    raw.write_bytes(b"changed")
    with pytest.raises(ValueError):
        module.audit(Path(inv["path"]), inv["sha256"], Path(base_record["path"]))
