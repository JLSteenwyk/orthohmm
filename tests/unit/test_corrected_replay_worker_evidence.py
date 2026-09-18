import json

import pytest

from benchmark_tools.checked_replay_payload_worker import corrected_evidence
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def evidence(tmp_path):
    checkpoint = tmp_path / "checkpoint"
    checkpoint.mkdir()
    names = checkpoint / "gene_names.txt"
    names.write_text("a\nb\n")
    manifest = checkpoint / "manifest.json"
    manifest.write_text("{}")
    admission = tmp_path / "admission.json"
    admission.write_text(json.dumps({"status": "corrected_high_sensitivity_native_evidence_admitted",
        "accuracy_evaluated": False, "checkpoint_manifest": record(manifest),
        "content": {"genes": 984137, "species_ownership": {str(i): i for i in range(78)},
                    "checked_records": [record(names)]}}))
    plan = tmp_path / "plan.json"
    plan.write_text(json.dumps({"status": "corrected_replay_command_frozen_unrun",
        "accuracy_evaluated": False, "execution_authorized": False,
        "admission": record(admission), "checkpoint_manifest": record(manifest),
        "checked_records": [record(names), record(manifest), record(admission)]}))
    return plan, admission, names


def test_corrected_evidence_binds_names(tmp_path):
    plan, admission, names = evidence(tmp_path)
    _, plan_record, admitted, name_record = corrected_evidence(plan, record(plan)["sha256"])
    assert plan_record == record(plan) and admitted == record(admission) and name_record == record(names)


@pytest.mark.parametrize("problem", ["plan_hash", "plan_status", "old_count", "wrong_species",
    "missing_names", "duplicate_names_record", "changed_names", "changed_admission", "checkpoint"])
def test_reject_changed_evidence(tmp_path, problem):
    path, admission, names = evidence(tmp_path)
    plan = json.loads(path.read_text())
    data = json.loads(admission.read_text())
    if problem == "plan_status":
        plan["status"] = "preparing"
    elif problem == "old_count":
        data["content"]["genes"] = 976504
    elif problem == "wrong_species":
        data["content"]["species_ownership"] = {}
    elif problem == "missing_names":
        data["content"]["checked_records"] = []
    elif problem == "duplicate_names_record":
        data["content"]["checked_records"] *= 2
    elif problem == "checkpoint":
        plan["checkpoint_manifest"] = {**plan["checkpoint_manifest"], "sha256": "changed"}
    admission.write_text(json.dumps(data))
    plan["admission"] = record(admission)
    plan["checked_records"][-1] = record(admission)
    path.write_text(json.dumps(plan))
    expected = record(path)["sha256"]
    if problem == "changed_names":
        names.write_text("a\nforeign\n")
    elif problem == "changed_admission":
        admission.write_text("{}")
    elif problem == "plan_hash":
        expected = "wrong"
    with pytest.raises(ValueError):
        corrected_evidence(path, expected)
