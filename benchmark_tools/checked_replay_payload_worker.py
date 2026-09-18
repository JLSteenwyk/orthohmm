"""Checked native worker for a manifest-pinned payload from a full QfO replay."""

import argparse
import json
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from prepare_ob_candidate_neighborhood import check, record

ADMISSION_SHA = "78f51f5ce703caf39307c5e518ad737acdf655dbe0926a34f2c20bbdec6d03f1"
STAGES = ("initial", "multipass", "profile_base", "profile_expanded")
FILES = ("gene_names.txt", "sources.npy", "targets.npy", "weights.npy", "metadata.json")


def corrected_evidence(plan_path, expected_sha):
    plan_record = record(plan_path)
    if plan_record["sha256"] != expected_sha:
        raise ValueError("Corrected replay plan hash differs")
    plan = json.loads(plan_path.read_text())
    if (plan.get("status") != "corrected_replay_command_frozen_unrun"
            or plan.get("accuracy_evaluated") is not False or plan.get("execution_authorized") is not False):
        raise ValueError("Unexpected corrected replay plan state")
    check(plan["admission"])
    admission = json.loads(Path(plan["admission"]["path"]).read_text())
    if (admission.get("status") != "corrected_high_sensitivity_native_evidence_admitted"
            or admission.get("accuracy_evaluated") is not False
            or admission["content"]["genes"] != 984137
            or len(admission["content"]["species_ownership"]) != 78
            or admission["checkpoint_manifest"] != plan["checkpoint_manifest"]):
        raise ValueError("Corrected HMM admission differs")
    names_path = Path(plan["checkpoint_manifest"]["path"]).parent / "gene_names.txt"
    matches = [r for r in admission["content"]["checked_records"] if r["path"] == str(names_path)]
    if len(matches) != 1:
        raise ValueError("Require uniquely admitted corrected gene-name record")
    for item in [plan_record, plan["admission"], plan["checkpoint_manifest"], *plan["checked_records"], matches[0]]:
        check(item)
    return plan, plan_record, plan["admission"], matches[0]


def validate_payload(manifest, payload, admitted_names):
    if (manifest["stage"] not in STAGES or manifest["index"] != STAGES.index(manifest["stage"])
            or manifest["accuracy_evaluated"] is not False
            or manifest["inputs"] != [record(payload / name) for name in FILES]):
        raise ValueError("Payload inventory or stage differs")
    names = manifest["inputs"][0]
    if any(names[key] != admitted_names[key] for key in ("bytes", "sha256")):
        raise ValueError("Replay gene order differs from admitted QfO graph")
    metadata = json.loads((payload / "metadata.json").read_text())
    expected = {"cpm_resolution": .1, "seed": 4, "include_isolates": True,
                "output_directory": manifest["output_directory"]}
    if metadata != expected or not Path(metadata["output_directory"]).is_absolute():
        raise ValueError("Payload settings or output directory differ")
    for name in ("native_boundary.json", "constructor_adapter.json", "worker_before.json", "checked_payload_provenance.json"):
        if (payload / name).exists():
            raise FileExistsError("Refuse already observed worker payload")
    return metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--payload", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--corrected-plan", type=Path)
    parser.add_argument("--corrected-plan-sha256")
    parser.add_argument("--sequence-plan", type=Path)
    parser.add_argument("--sequence-plan-sha256")
    parser.add_argument("--sequence-variant", choices=("all_hits", "top100"))
    args = parser.parse_args()
    if (args.corrected_plan is None) != (args.corrected_plan_sha256 is None):
        parser.error("Corrected plan and hash must be supplied together")
    sequence_flags = (args.sequence_plan, args.sequence_plan_sha256, args.sequence_variant)
    if any(value is not None for value in sequence_flags) and not all(value is not None for value in sequence_flags):
        parser.error("Sequence plan, hash and variant must be supplied together")
    if args.corrected_plan is not None and args.sequence_plan is not None:
        parser.error("HMM and sequence provenance are mutually exclusive")
    root, payload = args.root.resolve(), args.payload.resolve()
    manifest_record = record(args.manifest)
    if manifest_record["sha256"] != args.manifest_sha256:
        raise ValueError("Changed parent payload manifest")
    manifest = json.loads(args.manifest.read_text())
    corrected_plan_record = None
    sequence_plan_record = None
    if args.sequence_plan is not None:
        from sequence_graph_evidence import sequence_evidence, sequence_settings
        plan, sequence_plan_record, admission_record, names_record = sequence_evidence(
            args.sequence_plan.resolve(), args.sequence_plan_sha256, args.sequence_variant)
        metadata = validate_payload(manifest, payload, names_record)
        sequence_settings(manifest, metadata, plan["variants"][args.sequence_variant])
    elif args.corrected_plan is not None:
        plan, corrected_plan_record, admission_record, names_record = corrected_evidence(
            args.corrected_plan.resolve(), args.corrected_plan_sha256)
        if manifest["output_directory"] != str(Path(plan["output_root"]) / "replay"):
            raise ValueError("Corrected clustering output outside frozen replay")
        validate_payload(manifest, payload, names_record)
    else:
        admission_path = root / "benchmark_tools/results/qfo_checked_repeats_verified_20260917.json"
        admission_record = record(admission_path)
        if admission_record["sha256"] != ADMISSION_SHA:
            raise ValueError("Changed checked-repeat admission")
        admission = json.loads(admission_path.read_text())
        if admission["status"] != "checked_repeats_verified" or admission["all_three_partitions_equal"] is not True:
            raise ValueError("Checked initial-graph repeats were not admitted")
        for item in admission["provenance_checked"]:
            check(item)
        validate_payload(manifest, payload, admission["native_report"]["graph_inputs"][0])
        if manifest["index"] == 0:
            for actual, expected in zip(manifest["inputs"][:4], admission["native_report"]["graph_inputs"]):
                if any(actual[key] != expected[key] for key in ("bytes", "sha256")):
                    raise ValueError("Regenerated initial graph differs from admitted input")
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    if Path.cwd() != launcher:
        raise ValueError("Wrong frozen worker working directory")
    from repeat_qfo_saved_graph import worker, set_worker_affinity
    from checked_python_pair_worker import python_pair_constructor
    set_worker_affinity([min(os.sched_getaffinity(0))])
    import numpy as np
    arrays = [np.load(payload / name, mmap_mode="r", allow_pickle=False) for name in FILES[1:4]]
    if (any(a.ndim != 1 or a.shape != arrays[0].shape for a in arrays)
            or [a.dtype for a in arrays] != [np.dtype("int32"), np.dtype("int32"), np.dtype("float64")]
            or not len(arrays[0])):
        raise ValueError("Invalid replay graph arrays")
    vertices = len((payload / "gene_names.txt").read_text().splitlines())
    if any(np.any(a < 0) or np.any(a >= vertices) for a in arrays[:2]) or not np.isfinite(arrays[2]).all():
        raise ValueError("Invalid replay endpoints or weights")
    del arrays
    import igraph
    provenance = {"source": record(__file__), "manifest": manifest_record, "admission": admission_record,
                  "inputs": manifest["inputs"], "stage": manifest["stage"], "accuracy_evaluated": False,
                  "helpers": [record(Path(__file__).with_name(name)) for name in
                              ("repeat_qfo_saved_graph.py", "checked_python_pair_worker.py", "probe_leiden_boundary.py")]}
    if corrected_plan_record is not None:
        provenance["corrected_plan"] = corrected_plan_record
    if sequence_plan_record is not None:
        provenance.update(sequence_plan=sequence_plan_record, sequence_variant=args.sequence_variant)
        provenance["helpers"].append(record(Path(__file__).with_name("sequence_graph_evidence.py")))
    (payload / "checked_payload_provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    with python_pair_constructor(igraph, payload / "constructor_adapter.json"):
        worker(launcher, payload, native_boundary=True)
    raise RuntimeError("Frozen native worker unexpectedly returned")


if __name__ == "__main__":
    main()
