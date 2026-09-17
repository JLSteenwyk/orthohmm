"""Verify CPM-specific candidate seeds, fixed rules and complete merge reconstruction."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
from benchmark_tools.prepare_ob_cpm_candidates import REPLAY_SHA, LABELS, seed_arms
from benchmark_tools.prepare_ob_candidate_neighborhood import PREPARED_SHA, check, record
from benchmark_tools.replay_phylogeny import load_membership_constraints
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.trace_ob_families import partition, validate_merge_reconstruction
from benchmark_tools.verify_ygob_validation import require_completed_job

PREPARATION_SHA = "05c1b299c4992d124fab0c760125af45081937283652c67058dd0d6b053cd756"


def check_arm(row, label, seed, parameters):
    if (row["label"] != label or row["seed_partition"] != seed or row["candidate_expansion"] is not True
            or row["expansion"]["parameters"] != parameters or row["expansion"]["profile"] != "satellite_v2"):
        raise ValueError("Changed CPM seed or candidate-expansion rules")


def admit(root, output):
    if output.exists():
        raise FileExistsError(output)
    accounting = subprocess.check_output(["sacct", "-j", "21322", "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = require_completed_job(accounting, 21322)
    path = root / "benchmarks/results/ob_cpm_candidates_v1/results.json"
    report = read_frozen(path, PREPARATION_SHA)
    if (report["status"] != "cpm_candidates_prepared_unscored" or report["accuracy_evaluated"] is not False
            or report["job_id"] != "21322" or [row["label"] for row in report["arms"]] != list(LABELS)):
        raise ValueError("Incomplete or changed CPM candidate preparation")
    results = root / "benchmark_tools/results"
    replay_path = results / "ob_cpm_replay_verified_20260916.json"
    replay = read_frozen(replay_path, REPLAY_SHA)
    seeds = seed_arms(replay)
    prepared_path = results / "orthobench_factorial_prepared_20260916.json"
    prepared = read_frozen(prepared_path, PREPARED_SHA)
    baseline = prepared["candidate_arms"]["p1_c1"]
    expected_inputs = [record(replay_path), record(prepared_path), *replay["provenance_checked"],
                       baseline["candidate_partition"], baseline["membership_constraints"]]
    if report["inputs"] != expected_inputs:
        raise ValueError("Candidate preparation input inventory changed")
    records = [record(path), report["source"], *report["inputs"], *report["helpers"]]
    for row in report["arms"]:
        records.extend([row["candidate_partition"], row["membership_constraints"]])
    for item in records:
        check(item)
    universe = set()
    for item in prepared["fasta_inputs"]:
        for sequence in SeqIO.parse(item["path"], "fasta"):
            if sequence.id in universe:
                raise ValueError("Duplicate FASTA ID")
            universe.add(sequence.id)
    for row, (label, seed) in zip(report["arms"], seeds):
        check_arm(row, label, seed, baseline["expansion"]["parameters"])
        seed_groups, _ = partition(Path(seed["path"]), "plain", universe)
        candidate_path = Path(row["candidate_partition"]["path"])
        candidates, _ = partition(candidate_path, "plain", universe)
        events = load_membership_constraints(Path(row["membership_constraints"]["path"]), candidate_path)
        validate_merge_reconstruction(events, seed_groups, candidates, set())
        if (row["genes"] != len(universe) or row["seed_groups"] != len(seed_groups)
                or row["candidate_groups"] != len(candidates) or row["reconstructed_merges"] != len(events)
                or row["expansion"]["candidate_families"] != len(candidates)
                or row["expansion"]["merges"] != len(events) or len(seed_groups) - len(events) != len(candidates)):
            raise ValueError("Native counts differ from reconstructed candidate partition")
        if label == "control":
            if row.get("baseline_byte_equivalent") is not True or any(row[key]["sha256"] != baseline[key]["sha256"]
                    for key in ("candidate_partition", "membership_constraints")):
                raise ValueError("Unchanged candidate control differs from baseline")
    for item in records:
        check(item)
    result = {"status": "cpm_candidates_verified_unscored", "accuracy_evaluated": False,
        "scheduler": scheduler, "arms": report["arms"], "source": record(__file__),
        "provenance_checked": records, "limitations": ["Native preparation and logged merges verified; search decisions not independently recalculated.",
            "Both CPM variants still require inferred phylogeny and native-output admission before six-variant accuracy scoring."]}
    output.mkdir(parents=True)
    (output / "results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    admit(args.root.resolve(), args.output.resolve())
