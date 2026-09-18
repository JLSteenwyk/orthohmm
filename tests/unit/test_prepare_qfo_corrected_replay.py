import json
from pathlib import Path
import sys

import pytest

from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.prepare_qfo_corrected_replay import command_for, prepare, validate_admission, METHOD


def fixture():
    primary = {"input_directory": "/corrected", "methods": {METHOD: {"output": "/new/hmm"}}}
    inputs = [{"path": f"/corrected/s{i}.fasta", "bytes": 10, "sha256": str(i)} for i in range(78)]
    manifest = {"path": "/new/hmm/orthohmm_working_res/high_sensitivity_checkpoint/manifest.json",
                "bytes": 100, "sha256": "new-checkpoint"}
    admission = {"status": "corrected_high_sensitivity_native_evidence_admitted", "accuracy_evaluated": False,
        "scheduler": {"State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon", "AllocCPUS": "32"},
        "content": {"status": "high_sensitivity_output_content_verified", "accuracy_evaluated": False,
                    "genes": 984137, "species_ownership": {r["path"]: i for i, r in enumerate(inputs)},
                    "numeric_checkpoint": {"manifest": manifest}},
        "checkpoint_manifest": manifest, "checked_records": [*inputs, manifest]}
    return admission, primary, inputs


def test_admitted_corrected_checkpoint():
    path, sha = validate_admission(*fixture())
    assert str(path).startswith("/new/hmm/") and sha == "new-checkpoint"


@pytest.mark.parametrize("change", ["old_count", "unadmitted", "wrong_checkpoint", "old_inputs",
    "missing_input_hash", "duplicate_species", "duplicate_input", "running", "scored"])
def test_reject_incompatible_evidence(change):
    admission, primary, inputs = fixture()
    if change == "old_count":
        admission["content"]["genes"] = 976504
    elif change == "unadmitted":
        admission["status"] = "process_succeeded_pending_native_admission"
    elif change == "wrong_checkpoint":
        admission["checkpoint_manifest"] = {**admission["checkpoint_manifest"], "path": "/old/manifest.json"}
    elif change == "old_inputs":
        primary["input_directory"] = "/old"
    elif change == "missing_input_hash":
        admission["checked_records"] = admission["checked_records"][1:]
    elif change == "duplicate_species":
        admission["content"]["species_ownership"][inputs[0]["path"]] = 1
    elif change == "duplicate_input":
        inputs[0] = inputs[1]
    elif change == "running":
        admission["scheduler"]["State"] = "RUNNING"
    else:
        admission["accuracy_evaluated"] = True
    with pytest.raises(ValueError):
        validate_admission(admission, primary, inputs)


def test_scientific_flags_unchanged():
    from benchmark_tools.run_qfo_publication_replay import command_for as original
    root, launcher, output = Path("/root"), Path("/launcher"), Path("/output")
    old = original(root, launcher, output)
    new = command_for(Path(sys.executable), launcher, output, Path("/corrected"), Path("/new/checkpoint"), "new-sha")
    for flag in ("--accuracy-checkpoint", "--checkpoint-sha256", "--fasta-directory"):
        old[old.index(flag) + 1] = new[new.index(flag) + 1]
    assert old == new
    assert "--official-benchmark" not in new


def test_prepare_manifest_with_real_file_hashes(tmp_path, monkeypatch):
    results = tmp_path / "benchmark_tools/results"
    results.mkdir(parents=True)
    fasta = tmp_path / "corrected"
    fasta.mkdir()
    inputs = []
    for i in range(78):
        path = fasta / f"s{i}.fasta"
        path.write_text(f">g{i}\nACDE\n")
        inputs.append(record(path))
    stage = fasta / "staging_manifest.json"
    stage.write_text(json.dumps({"input_fastas": inputs}))
    native = tmp_path / "native"
    checkpoint = native / "orthohmm_working_res/high_sensitivity_checkpoint"
    checkpoint.mkdir(parents=True)
    manifest = checkpoint / "manifest.json"
    manifest.write_text("{}")
    admission, primary, _ = fixture()
    primary.update(input_directory=str(fasta), inputs=[record(stage), *inputs])
    primary["methods"][METHOD] = {"output": str(native), "native_argv": [sys.executable]}
    admission["checkpoint_manifest"] = record(manifest)
    admission["content"].update(species_ownership={r["path"]: i for i, r in enumerate(inputs)},
        numeric_checkpoint={"manifest": record(manifest)}, checked_records=[record(manifest), *inputs])
    admission["checked_records"] = [record(manifest), *inputs]
    plan = results / "qfo_corrected_primary_commands_20260918.json"
    plan.write_text(json.dumps(primary))
    source = tmp_path / "admission.json"
    source.write_text(json.dumps(admission))
    monkeypatch.setattr("benchmark_tools.prepare_qfo_corrected_replay.PLAN_SHA", record(plan)["sha256"])
    monkeypatch.setattr("benchmark_tools.prepare_qfo_corrected_replay.verify", lambda *args: {"status": "test-runtime"})
    destination = tmp_path / "prepared.json"
    result = prepare(tmp_path, source, record(source)["sha256"], tmp_path / "replay", destination)
    assert result["execution_authorized"] is False
    assert result["accuracy_evaluated"] is False
    assert result["checkpoint_manifest"] == record(manifest)
    assert len(result["input_fastas"]) == 78
    assert json.loads(destination.read_text()) == result
    assert not (tmp_path / "replay").exists()
    with pytest.raises(FileExistsError):
        prepare(tmp_path, source, record(source)["sha256"], tmp_path / "replay", destination)
