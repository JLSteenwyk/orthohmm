"""Recheck retained corrected replay stages, without admitting the parent job."""

import json
from pathlib import Path
import sys

from benchmark_tools.checked_replay_payload_worker import FILES, STAGES
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.validate_checked_replay_payload import validate


def audit(root, executor, directory, worker, replay, plan_record, sequence_variant=None, *, cpm_arm=None):
    if sequence_variant is not None and sequence_variant not in ("all_hits", "top100"):
        raise ValueError("Unknown sequence variant")
    if cpm_arm is not None:
        if cpm_arm not in ("control", "cpm_low", "cpm_high") or sequence_variant is not None:
            raise ValueError("CPM stage audit requires a recognized corrected HMM arm")
    expected_stages = STAGES if sequence_variant is None else STAGES[:2]
    calls = worker["calls"]
    if (len(calls) != len(expected_stages) or [(c["index"], c["stage"]) for c in calls] != list(enumerate(expected_stages))
            or any(type(c["index"]) is not int or c["status"] != "checked"
                   or type(c["exit_code"]) is not int or c["exit_code"] != 0
                   or c["accuracy_evaluated"] is not False for c in calls)):
        raise ValueError("Incomplete or reordered checked clustering stages")
    labels = ["multipass", "multipass_refined", "strict_profiles", "strict_profiles_refined"]
    if sequence_variant is not None:
        labels = labels[:2]
    if [s["label"] for s in replay["stages"]] != labels:
        raise ValueError("Unexpected replay stage inventory")
    records, summaries = [], []
    for row in calls:
        stage = directory / "clustering" / f"cluster_{row['index']}_{row['stage']}"
        payload = stage / "payload"
        execution_record = record(stage / "execution.json")
        if json.loads((stage / "execution.json").read_text()) != row:
            raise ValueError("Preserved execution differs from parent worker")
        manifest_record = record(stage / "payload_manifest.json")
        if row["manifest"] != manifest_record:
            raise ValueError("Stage manifest changed")
        manifest = json.loads((stage / "payload_manifest.json").read_text())
        if (manifest["index"] != row["index"] or manifest["stage"] != row["stage"]
                or manifest["accuracy_evaluated"] is not False
                or manifest["output_directory"] != str(directory / "replay")):
            raise ValueError("Stage manifest identity differs")
        inputs = [record(payload / name) for name in FILES]
        if manifest["inputs"] != inputs or len(manifest["original_payload"]) != len(FILES):
            raise ValueError("Preserved payload inventory differs")
        original_dirs = {Path(r["path"]).parent for r in manifest["original_payload"]}
        if len(original_dirs) != 1:
            raise ValueError("Original payload has multiple directories")
        original = next(iter(original_dirs))
        if manifest["original_command"] != [sys.executable, "-m", "orthohmm.leiden_worker", str(original)]:
            raise ValueError("Original worker command differs")
        for name, prior, copied in zip(FILES, manifest["original_payload"], inputs):
            if (Path(prior["path"]).name != name
                    or any(prior[k] != copied[k] for k in ("bytes", "sha256"))):
                raise ValueError("Original/copied payload content differs")
        expected = [sys.executable, str(executor / "benchmark_tools/checked_replay_payload_worker.py"),
                    "--root", str(root), "--payload", str(payload), "--manifest", manifest_record["path"],
                    "--manifest-sha256", manifest_record["sha256"]]
        if sequence_variant is None:
            expected += ["--corrected-plan", plan_record["path"], "--corrected-plan-sha256", plan_record["sha256"]]
            if cpm_arm is not None:
                expected += ["--cpm-arm", cpm_arm]
        else:
            expected += ["--sequence-plan", plan_record["path"], "--sequence-plan-sha256", plan_record["sha256"],
                         "--sequence-variant", sequence_variant]
        if row["command"] != expected:
            raise ValueError("Checked worker command differs")
        overrides = {k: "1" for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}
        inherited = {**overrides, "OMP_NUM_THREADS": "32" if row["index"] >= 2 else "1"}
        if row["thread_environment"] != {"inherited": inherited, "child_overrides": overrides}:
            raise ValueError("Clustering thread isolation differs")
        partition = record(stage / "partition.txt")
        if row["partition"] != partition:
            raise ValueError("Retained partition changed")
        provenance = (dict(corrected_plan=plan_record) if sequence_variant is None else
                      dict(sequence_plan=plan_record, sequence_variant=sequence_variant))
        if cpm_arm is not None:
            provenance["cpm_arm"] = cpm_arm
        fresh = validate(payload, manifest, root, executor, retained_partition=partition, **provenance)
        # The in-run callback preceded the retained copy. All other evidence
        # must match the fresh retrospective check exactly.
        original_validation = {**fresh, "provenance_checked": fresh["provenance_checked"][:-1]}
        if fresh["provenance_checked"][-1] != partition or row["validation"] != original_validation:
            raise ValueError("Recorded validation differs from fresh native audit")
        if fresh["genes"] != 984137 or fresh["saved_graph"]["vertices"] != 984137:
            raise ValueError("Wrong corrected gene universe")
        summaries.append({"index": row["index"], "stage": row["stage"], "partition": partition,
                          "genes": fresh["genes"], "groups": fresh["groups"], "saved_graph": fresh["saved_graph"]})
        records.extend([execution_record, manifest_record, partition, record(stage / "worker.log"),
                        *fresh["provenance_checked"]])
    if (replay["counts"]["rbnh_edges"] != summaries[0]["saved_graph"]["edges"]
            or replay["counts"]["multipass_edges"] != summaries[1]["saved_graph"]["edges"]):
        raise ValueError("Replay edge counts differ from native payloads")
    stages = {s["label"]: s["output"] for s in replay["stages"]}
    copies = ((1, "multipass"), (3, "strict_profiles")) if sequence_variant is None else ((1, "multipass"),)
    for index, label in copies:
        output = stages[label]
        check(output)
        if any(output[k] != summaries[index]["partition"][k] for k in ("bytes", "sha256")):
            raise ValueError("Replay stage not copied from corresponding checked partition")
        records.append(output)
    unique = {}
    for item in records:
        if item["path"] in unique and unique[item["path"]] != item:
            raise ValueError("Conflicting stage provenance records")
        unique[item["path"]] = item
    for item in unique.values():
        check(item)
    status = "corrected_replay_stages_verified" if sequence_variant is None else "sequence_replay_stages_verified"
    if cpm_arm is not None:
        status = "cpm_replay_stages_verified"
    return {"status": status,
            "accuracy_evaluated": False,
            "clustering": summaries, "checked_records": list(unique.values()),
            "limitations": ["Stage evidence only: parent scheduler, plan, runtime and refined outputs require separate admission.",
                            "Recorded native observations are rechecked; this is not retrospective live-memory inspection."]}
