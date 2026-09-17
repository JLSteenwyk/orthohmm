"""Expand each admitted CPM-specific HMM seed partition with fixed satellite rules."""

import argparse
import json
import os
from pathlib import Path
import pickle
import sys

REPLAY_SHA = "ceb8f9fe7c12b35317cde99c3d027d4d1080d398403afc04e57fb029ce1b460a"
LABELS = ("control", "cpm_low", "cpm_high")


def seed_arms(replay):
    if (replay["status"] != "cpm_replay_panel_verified_unscored" or replay["accuracy_evaluated"] is not False
            or [row["label"] for row in replay["arms"]] != list(LABELS)):
        raise ValueError("Require complete admitted CPM replay panel")
    return [(row["label"], row["stages"]["strict_profiles_refined"]["output"]) for row in replay["arms"]]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    root, output = args.root.resolve(), args.output.resolve()
    if output.exists():
        raise FileExistsError(output)
    frozen = root / "benchmarks/work/publication_method_native_v2"
    sys.path.insert(0, str(frozen))
    import orthohmm.orthohmm as engine
    if Path(engine.__file__).resolve() != frozen / "orthohmm/orthohmm.py":
        raise ValueError("Wrong scientific engine")
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    import numpy as np
    from Bio import SeqIO
    from benchmark_tools.prepare_ob_candidate_neighborhood import PREPARED_SHA, check, record
    from benchmark_tools.prepare_orthobench_factorial import indexed_species, prepare_partition
    from benchmark_tools.replay_high_sensitivity import hit_arrays
    from benchmark_tools.replay_phylogeny import load_membership_constraints
    from benchmark_tools.run_simulation_methods import read_frozen
    from benchmark_tools.trace_ob_families import partition, validate_merge_reconstruction
    results = root / "benchmark_tools/results"
    replay_path = results / "ob_cpm_replay_verified_20260916.json"
    replay = read_frozen(replay_path, REPLAY_SHA)
    seeds = seed_arms(replay)
    prepared_path = results / "orthobench_factorial_prepared_20260916.json"
    prepared = read_frozen(prepared_path, PREPARED_SHA)
    baseline = prepared["candidate_arms"]["p1_c1"]
    records = [record(replay_path), record(prepared_path), *replay["provenance_checked"],
               baseline["candidate_partition"], baseline["membership_constraints"]]
    for item in records:
        check(item)
    owners = {}
    for item in prepared["fasta_inputs"]:
        for sequence in SeqIO.parse(item["path"], "fasta"):
            if sequence.id in owners:
                raise ValueError("Duplicate FASTA identifier")
            owners[sequence.id] = Path(item["path"]).name
    with Path(prepared["cache"]["path"]).open("rb") as handle:
        payload = pickle.load(handle)
    names, species = indexed_species(payload, owners)
    hits = hit_arrays(payload["all_hits"], {name: index for index, name in enumerate(names)})
    if not np.isfinite(hits[2]).all():
        raise ValueError("Nonfinite normalized hits")
    del payload
    output.mkdir(parents=True)
    report = {"status": "preparing_unscored", "accuracy_evaluated": False,
        "job_id": os.environ.get("SLURM_JOB_ID"), "source": record(__file__), "inputs": records, "arms": [],
        "helpers": [record(Path(__file__).with_name(name)) for name in
            ("prepare_orthobench_factorial.py", "replay_high_sensitivity.py", "replay_phylogeny.py", "trace_ob_families.py")],
        "limitations": ["Only CPM-specific HMM seeds vary; candidate parameters use the unchanged satellite_v2 profile.",
            "Complete inferred phylogeny and independent admission remain required before accuracy evaluation."]}
    try:
        for label, seed in seeds:
            row = prepare_partition(Path(seed["path"]), output / label, names, species, hits, True,
                                    engine._expand_phylogeny_candidates, load_membership_constraints)
            if row["expansion"]["parameters"] != baseline["expansion"]["parameters"]:
                raise ValueError("Candidate expansion parameters changed")
            seed_groups, _ = partition(Path(seed["path"]), "plain", set(names))
            candidate_path = Path(row["candidate_partition"]["path"])
            candidates, _ = partition(candidate_path, "plain", set(names))
            events = load_membership_constraints(Path(row["membership_constraints"]["path"]), candidate_path)
            validate_merge_reconstruction(events, seed_groups, candidates, set())
            row.update(label=label, genes=len(names), seed_groups=len(seed_groups), candidate_groups=len(candidates),
                       reconstructed_merges=len(events))
            if label == "control":
                if any(row[key]["sha256"] != baseline[key]["sha256"] for key in ("candidate_partition", "membership_constraints")):
                    raise ValueError("CPM candidate control differs from frozen baseline")
                row["baseline_byte_equivalent"] = True
            report["arms"].append(row)
            (output / "progress.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        for item in records:
            check(item)
        report["status"] = "cpm_candidates_prepared_unscored"
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        (output / "results.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
