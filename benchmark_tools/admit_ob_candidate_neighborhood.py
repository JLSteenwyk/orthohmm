"""Independently admit the fixed, unscored candidate-threshold preparation."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
from benchmark_tools.prepare_ob_candidate_neighborhood import ARMS, PREPARED_SHA, check, record
from benchmark_tools.replay_phylogeny import load_membership_constraints
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.trace_ob_families import partition, validate_merge_reconstruction
from benchmark_tools.verify_ygob_validation import require_completed_job

MANIFEST_SHA = "29db21af15d34f83b3fb9a0aa4b96f1871bf4a10b021fff8cc98827382249ec7"


def validate_plan(report, baseline):
    if (report["status"] != "five_candidate_arms_prepared_unscored"
            or report["job_id"] != "21314" or report["accuracy_evaluated"] is not False
            or report["publication_ready"] is not False):
        raise ValueError("Unexpected preparation identity or status")
    if [row["label"] for row in report["arms"]] != [label for label, _ in ARMS]:
        raise ValueError("Changed candidate arm plan")
    for row, (_, delta) in zip(report["arms"], ARMS):
        if (row["delta"] != delta or row["applied_parameters"] != {**baseline, **delta}
                or row["engine_calls"] != 1
                or row["engine_fixed_profile_report"]["parameters"] != baseline):
            raise ValueError("Changed applied or nominal parameters")


def admit(root, output):
    if output.exists():
        raise FileExistsError(output)
    path = root / "benchmarks/results/ob_candidate_neighborhood_v1/manifest.json"
    report = read_frozen(path, MANIFEST_SHA)
    accounting = subprocess.check_output(["sacct", "-j", "21314", "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = require_completed_job(accounting, 21314)
    prepared = read_frozen(Path(report["prepared_manifest"]["path"]), PREPARED_SHA)
    baseline = prepared["candidate_arms"]["p1_c1"]
    validate_plan(report, baseline["expansion"]["parameters"])
    expected_inputs = [*prepared["core_sources"], *prepared["fasta_inputs"], prepared["cache"],
        baseline["seed_partition"], baseline["candidate_partition"], baseline["membership_constraints"]]
    if report["inputs"] != expected_inputs:
        raise ValueError("Changed scientific inputs")
    records = [record(path), report["source"], report["prepared_manifest"],
               *report["inputs"], *report["helpers"]]
    for row in report["arms"]:
        records.extend([row["partition"], row["constraints"]])
    for item in records:
        check(item)
    universe = set()
    for item in prepared["fasta_inputs"]:
        for sequence in SeqIO.parse(item["path"], "fasta"):
            if sequence.id in universe:
                raise ValueError("Duplicate FASTA identifier")
            universe.add(sequence.id)
    seeds, _ = partition(Path(baseline["seed_partition"]["path"]), "plain", universe)
    verified = []
    for row in report["arms"]:
        candidate_path = Path(row["partition"]["path"])
        candidates, _ = partition(candidate_path, "plain", universe)
        events = load_membership_constraints(Path(row["constraints"]["path"]), candidate_path)
        validate_merge_reconstruction(events, seeds, candidates, set())
        engine = row["engine_fixed_profile_report"]
        if (engine["candidate_families"] != len(candidates) or engine["merges"] != len(events)
                or engine["seed_families"] != len(seeds) or len(seeds) - len(events) != len(candidates)):
            raise ValueError("Reported counts differ from reconstructed partition")
        if row["label"] == "control":
            if (row.get("baseline_byte_equivalent") is not True
                    or row["partition"]["sha256"] != baseline["candidate_partition"]["sha256"]
                    or row["constraints"]["sha256"] != baseline["membership_constraints"]["sha256"]):
                raise ValueError("Control differs from frozen baseline")
        verified.append({"label": row["label"], "genes": len(universe), "seed_groups": len(seeds),
            "candidate_groups": len(candidates), "reconstructed_merges": len(events)})
    for item in records:
        check(item)
    result = {"status": "candidate_neighborhood_preparation_verified_unscored",
        "accuracy_evaluated": False, "publication_ready": False,
        "source": record(__file__), "preparation": record(path), "scheduler": scheduler,
        "arms": report["arms"], "verification": verified, "provenance_checked": records,
        "helpers": [record(Path(__file__).with_name(name)) for name in
            ("prepare_ob_candidate_neighborhood.py", "trace_ob_families.py", "replay_phylogeny.py")],
        "limitations": ["Merge reconstruction validates retained memberships, not independent recomputation of search evidence or threshold decisions.",
            "Four changed candidate arms require downstream inferred-tree evaluation; two CPM variants remain pending.",
            "No scores inspected, default promotion or claim of robustness from preparation alone."]}
    output.mkdir(parents=True)
    (output / "results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    admit(args.root.resolve(), args.output.resolve())
