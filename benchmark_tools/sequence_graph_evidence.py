"""Validate sequence-control provenance without substituting HMM admission."""

import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from prepare_ob_candidate_neighborhood import check, record


def sequence_evidence(path, sha256, label):
    plan_record = record(path)
    if plan_record["sha256"] != sha256:
        raise ValueError("Sequence graph plan hash differs")
    plan = json.loads(path.read_text())
    if (plan["status"] != "corrected_sequence_graph_commands_frozen_unrun"
            or set(plan["variants"]) != {"all_hits", "top100"}
            or label not in plan["variants"]
            or any(plan[key] is not False for key in (
                "execution_authorized", "accuracy_evaluated", "graph_feasibility_admitted"))):
        raise ValueError("Wrong sequence graph plan")
    variant = plan["variants"][label]
    checkpoint = Path(variant["checkpoint_manifest"]["path"]).parent
    output = Path(variant["output_root"])
    # Keep worker validation independent of the already-imported frozen benchmark_tools package.
    command = [sys.executable, str(Path(plan["cwd"]) / "benchmark_tools/replay_high_sensitivity.py"),
               "--accuracy-checkpoint", str(checkpoint), "--checkpoint-sha256", variant["checkpoint_manifest"]["sha256"],
               "--output-directory", str(output / "replay"), "--json", str(output / "replay.json"),
               "--cpu", "32", "--matrix", "BLOSUM62", "--cpm-resolution", "0.1", "--leiden-seed", "4"]
    if (not output.is_absolute() or not checkpoint.is_absolute()
            or variant["expected_genes"] != 984137 or variant["expected_species"] != 78
            or variant["cpu"] != 32 or type(variant["requested_memory_gib"]) is not int
            or variant["requested_memory_gib"] <= 0
            or variant["cap"] != (None if label == "all_hits" else 100)
            or variant["expected_clustering_calls"] != ["initial", "multipass"]
            or variant["expected_stages"] != ["multipass", "multipass_refined"]
            or variant["gene_names"]["path"] != str(checkpoint / "gene_names.txt")
            or variant["native_command"] != command):
        raise ValueError("Sequence variant settings or identities differ")
    admission_record = plan["payload_report"]
    check(admission_record)
    admission = json.loads(Path(admission_record["path"]).read_text())
    if (admission["status"] != "admitted_qfo_graph_payload_estimated"
            or admission["accuracy_evaluated"] is not False
            or admission["graph_launched"] is not False):
        raise ValueError("Sequence payload estimate not admitted")
    estimated = admission["variants"][label]
    if (estimated["checkpoint_audit"]["manifest"] != variant["checkpoint_manifest"]
            or variant["gene_names"] not in estimated["checkpoint_files"]
            or estimated["estimate"]["genes"] != 984137
            or estimated["estimate"]["species_slot_extent"] != 78):
        raise ValueError("Sequence checkpoint differs from admitted estimate")
    for item in (admission_record, variant["checkpoint_manifest"], variant["gene_names"]):
        if item not in plan["checked_records"]:
            raise ValueError("Sequence evidence absent from checked inventory")
    for item in [plan_record, *plan["checked_records"], plan["source"], *plan["helpers"], plan["python"]]:
        check(item)
    return plan, plan_record, admission_record, variant["gene_names"]


def sequence_settings(manifest, metadata, variant):
    if (manifest["stage"] not in ("initial", "multipass")
            or manifest["index"] != ("initial", "multipass").index(manifest["stage"])
            or manifest["output_directory"] != str(Path(variant["output_root"]) / "replay")
            or metadata != {"cpm_resolution": .1, "seed": 4, "include_isolates": True,
                            "output_directory": manifest["output_directory"]}):
        raise ValueError("Wrong sequence-control stage or output/settings")
