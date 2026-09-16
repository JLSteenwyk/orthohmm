"""Prepare a fixed candidate-threshold neighborhood without evaluating accuracy."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import pickle
import shutil
import sys
import time

PREPARED_SHA = "5c325f4d77865e0c7571fe4bb4df0be0959977a49e1d22169f192518700f9382"
ARMS = (("control", {}), ("norm_low", {"min_norm": .024}), ("norm_high", {"min_norm": .036}),
        ("margin_low", {"min_margin": 1.2}), ("margin_high", {"min_margin": 1.8}))


def record(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def check(item):
    if record(item["path"]) != item:
        raise ValueError("Frozen input/source identity changed: " + item["path"])


def controlled_expansion(module, parameters, label, arguments):
    choices = dict(ARMS)
    if label not in choices or parameters["min_norm"] != .03 or parameters["min_margin"] != 1.5:
        raise ValueError("Unexpected candidate arm or baseline parameters")
    original = module.merge_supported_satellite_candidate_clusters
    calls = []
    applied = {**parameters, **choices[label]}
    def intercept(*args, **kwargs):
        observed = {k: v for k, v in kwargs.items() if k != "merge_trace"}
        if observed != parameters or kwargs.get("merge_trace") is None:
            raise ValueError("Frozen wrapper differs from admitted profile parameters")
        calls.append(dict(applied))
        return original(*args, **{**kwargs, **choices[label]})
    module.merge_supported_satellite_candidate_clusters = intercept
    try:
        result = module._expand_phylogeny_candidates(*arguments, profile="satellite_v2")
    finally:
        module.merge_supported_satellite_candidate_clusters = original
    if calls != [applied] or result["parameters"] != parameters:
        raise ValueError("Unexpected candidate-engine invocation or wrapper report")
    result.pop("_membership_constraints", None)
    return {"label": label, "delta": choices[label], "applied_parameters": applied,
            "engine_fixed_profile_report": result, "engine_calls": len(calls)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    root, output = args.root.resolve(), args.output.resolve()
    if output.exists():
        raise FileExistsError(output)
    prepared_path = root / "benchmark_tools/results/orthobench_factorial_prepared_20260916.json"
    if record(prepared_path)["sha256"] != PREPARED_SHA:
        raise ValueError("Changed frozen factorial preparation")
    prepared = json.loads(prepared_path.read_text())
    arm = prepared["candidate_arms"]["p1_c1"]
    inputs = [*prepared["core_sources"], *prepared["fasta_inputs"], prepared["cache"],
              arm["seed_partition"], arm["candidate_partition"], arm["membership_constraints"]]
    for item in inputs:
        check(item)
    frozen = root / "benchmarks/work/publication_method_native_v2"
    sys.path.insert(0, str(frozen))
    import orthohmm.orthohmm as engine
    if Path(engine.__file__).resolve() != frozen / "orthohmm/orthohmm.py":
        raise ValueError("Wrong frozen scientific engine")
    # Resolve analysis helpers from this executor after pinning the scientific package.
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from Bio import SeqIO
    from benchmark_tools.prepare_orthobench_factorial import indexed_species
    from benchmark_tools.replay_high_sensitivity import hit_arrays
    from benchmark_tools.audit_historical_profile_ablation import read_partition
    from benchmark_tools.replay_phylogeny import load_membership_constraints
    owners = {}
    for item in prepared["fasta_inputs"]:
        for protein in SeqIO.parse(item["path"], "fasta"):
            if protein.id in owners:
                raise ValueError("Duplicate FASTA identifier")
            owners[protein.id] = Path(item["path"]).name
    with Path(prepared["cache"]["path"]).open("rb") as handle:
        payload = pickle.load(handle)
    names, species = indexed_species(payload, owners)
    hits = hit_arrays(payload["all_hits"], {g: i for i, g in enumerate(names)})
    import numpy as np
    if not np.isfinite(hits[2]).all():
        raise ValueError("Nonfinite cached hits")
    del payload
    output.mkdir(parents=True)
    report = {"status": "preparing_unscored", "accuracy_evaluated": False, "publication_ready": False,
              "job_id": os.environ.get("SLURM_JOB_ID"), "source": record(__file__),
              "prepared_manifest": record(prepared_path), "inputs": inputs, "arms": [],
              "helpers": [record(Path(__file__).with_name(name)) for name in
                  ("prepare_orthobench_factorial.py", "replay_high_sensitivity.py", "replay_phylogeny.py")],
              "limitations": ["Only candidate min_norm and min_margin vary; HMM seed groups and cached hits stay fixed.",
                  "The frozen wrapper's nominal parameter report is retained separately from explicitly applied overrides.",
                  "Four threshold variants are part of a six-variant protocol; CPM variants and all downstream reconciliation/scoring remain pending.",
                  "No accuracy-based selection or default promotion; shared-node preparation times are not efficiency benchmarks."]}
    try:
        for label, _ in ARMS:
            directory = output / label
            working = directory / "orthohmm_working_res"
            working.mkdir(parents=True)
            partition_path = working / "orthohmm_edges_clustered.txt"
            shutil.copyfile(arm["seed_partition"]["path"], partition_path)
            started = time.monotonic()
            row = controlled_expansion(engine, arm["expansion"]["parameters"], label,
                                       (str(directory), names, species, hits))
            row["preparation_wall_s"] = time.monotonic() - started
            read_partition(partition_path, set(names))
            constraint_path = working / "phylogeny_candidate_merges.json"
            load_membership_constraints(constraint_path, partition_path)
            row["partition"] = record(partition_path)
            row["constraints"] = record(constraint_path)
            if label == "control":
                if (row["partition"]["sha256"] != arm["candidate_partition"]["sha256"]
                        or row["constraints"]["sha256"] != arm["membership_constraints"]["sha256"]):
                    raise ValueError("Unchanged control does not reproduce frozen candidate/merge bytes")
                row["baseline_byte_equivalent"] = True
            report["arms"].append(row)
            (output / "progress.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        for item in inputs:
            check(item)
        report["status"] = "five_candidate_arms_prepared_unscored"
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        (output / "manifest.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
