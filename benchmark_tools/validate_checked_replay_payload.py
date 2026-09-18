"""Validate a completed native payload before the frozen replay consumes it."""

import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from prepare_ob_candidate_neighborhood import check, record
from probe_leiden_boundary import saved_fingerprint


def validate(payload, manifest, root, executor, corrected_plan=None, retained_partition=None,
             sequence_plan=None, sequence_variant=None):
    if ((sequence_plan is None) != (sequence_variant is None)
            or (sequence_plan is not None and corrected_plan is not None)):
        raise ValueError("Require exactly one complete provenance mode")
    import numpy as np
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    observed = json.loads((payload / "worker_before.json").read_text())
    metadata = json.loads((payload / "metadata.json").read_text())
    if (observed["status"] != "before_native_clustering" or observed["accuracy_evaluated"] is not False
            or observed["metadata"] != metadata or observed["cwd"] != str(launcher)
            or observed["inputs"] != manifest["inputs"][:4] or len(observed["cpu_affinity"]) != 1
            or not observed["native_libraries"]):
        raise ValueError("Native worker identity differs")
    overrides = {"PYTHONPATH": str(launcher), "PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1",
                 "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    differences = {key: {"expected": value, "observed": observed["environment"].get(key)}
                   for key, value in overrides.items() if observed["environment"].get(key) != value}
    if differences:
        raise ValueError("Worker environment differs: " + json.dumps(differences, sort_keys=True))
    for name in ("orthohmm.leiden_worker", "orthohmm.externals", "orthohmm.helpers"):
        if observed["modules"][name] != record(launcher / (name.replace(".", "/") + ".py")):
            raise ValueError("Wrong scientific worker module")
    provenance = json.loads((payload / "checked_payload_provenance.json").read_text())
    if sequence_plan is not None:
        from sequence_graph_evidence import sequence_evidence, sequence_settings
        plan, plan_record, admission_record, names_record = sequence_evidence(
            Path(sequence_plan["path"]), sequence_plan["sha256"], sequence_variant)
        if (plan_record != sequence_plan or provenance.get("sequence_plan") != plan_record
                or provenance.get("sequence_variant") != sequence_variant or "corrected_plan" in provenance):
            raise ValueError("Sequence plan provenance differs")
        sequence_settings(manifest, metadata, plan["variants"][sequence_variant])
        if any(manifest["inputs"][0][key] != names_record[key] for key in ("bytes", "sha256")):
            raise ValueError("Sequence gene order differs after clustering")
    elif corrected_plan is None:
        admission_record = record(root / "benchmark_tools/results/qfo_checked_repeats_verified_20260917.json")
        if "corrected_plan" in provenance or "sequence_plan" in provenance or "sequence_variant" in provenance:
            raise ValueError("Unexpected corrected-release provenance")
    else:
        if "sequence_plan" in provenance or "sequence_variant" in provenance:
            raise ValueError("Unexpected sequence provenance on HMM payload")
        from checked_replay_payload_worker import corrected_evidence
        plan, plan_record, admission_record, names_record = corrected_evidence(
            Path(corrected_plan["path"]), corrected_plan["sha256"])
        if plan_record != corrected_plan or provenance.get("corrected_plan") != plan_record:
            raise ValueError("Corrected plan provenance differs")
        # The preflight payload validator also rejects observation files; here
        # only check the preserved names/settings after native execution.
        if any(manifest["inputs"][0][key] != names_record[key] for key in ("bytes", "sha256")):
            raise ValueError("Corrected gene order differs after clustering")
        if metadata != {"cpm_resolution": .1, "seed": 4, "include_isolates": True,
                        "output_directory": str(Path(plan["output_root"]) / "replay")}:
            raise ValueError("Corrected output/settings differ")
    helpers = [record(executor / "benchmark_tools" / name) for name in
               ("repeat_qfo_saved_graph.py", "checked_python_pair_worker.py", "probe_leiden_boundary.py")]
    if sequence_plan is not None:
        helpers.append(record(executor / "benchmark_tools/sequence_graph_evidence.py"))
    if (provenance["source"] != record(executor / "benchmark_tools/checked_replay_payload_worker.py")
            or provenance["helpers"] != helpers or observed["observer"] != helpers[0]
            or provenance["manifest"] != record(payload.parent / "payload_manifest.json")
            or provenance["inputs"] != manifest["inputs"] or provenance["stage"] != manifest["stage"]
            or provenance["accuracy_evaluated"] is not False
            or provenance["admission"] != admission_record):
        raise ValueError("Checked payload provenance differs")
    saved = saved_fingerprint(payload)
    boundary = json.loads((payload / "native_boundary.json").read_text())
    arguments = {"initial_membership": None, "weights": "weight", "n_iterations": 2,
                 "max_comm_size": 0, "seed": 4, "kwargs": {"resolution_parameter": .1},
                 "partition_type": "leidenalg.VertexPartition.CPMVertexPartition"}
    if boundary != {"accuracy_evaluated": False, "calls": [{"arguments": arguments, "before": saved,
            "saved": saved, "after": saved, "status": "optimizer_returned"}]}:
        raise ValueError("Native optimizer graph or settings differ")
    arrays = [np.load(payload / name, mmap_mode="r", allow_pickle=False) for name in ("sources.npy", "targets.npy")]
    digest = hashlib.sha256()
    for start in range(0, len(arrays[0]), 100000):
        digest.update(np.column_stack([a[start:start + 100000] for a in arrays]).tobytes(order="C"))
    adapter = json.loads((payload / "constructor_adapter.json").read_text())
    if adapter != {"format": "python_pairs", "accuracy_evaluated": False, "calls": [{"status": "constructor_returned",
            "shape": [saved["edges"], 2], "dtype": "int32", "ordered_input_bytes_sha256": digest.hexdigest(),
            "n": saved["vertices"], "directed": False}]}:
        raise ValueError("Constructor input observation differs")
    records = [*manifest["inputs"], *helpers, provenance["source"], provenance["manifest"], provenance["admission"],
               *observed["modules"].values(), *observed["native_libraries"], observed["python"],
               *[record(payload / name) for name in ("worker_before.json", "native_boundary.json",
                    "constructor_adapter.json", "checked_payload_provenance.json")]]
    for item in records:
        check(item)
    universe = set((payload / "gene_names.txt").read_text().splitlines())
    partition = Path(metadata["output_directory"]) / "orthohmm_working_res/orthohmm_edges_clustered.txt"
    if retained_partition is not None:
        partition = payload.parent / "partition.txt"
        if retained_partition["path"] != str(partition.resolve()):
            raise ValueError("Retained partition is outside this checked stage")
        check(retained_partition)
        records.append(retained_partition)
    seen, groups = set(), 0
    for line in partition.read_text().splitlines():
        members = line.split()
        if not members:
            continue
        genes = set(members)
        if len(genes) != len(members) or seen & genes or not genes <= universe:
            raise ValueError("Invalid partition membership")
        seen.update(genes)
        groups += 1
    if seen != universe or len(universe) != saved["vertices"]:
        raise ValueError("Incomplete partition or duplicate gene universe")
    if retained_partition is not None:
        check(retained_partition)
    return {"status": "payload_checked", "accuracy_evaluated": False, "groups": groups, "genes": len(seen),
            "saved_graph": saved, "worker": observed, "provenance_checked": records}
