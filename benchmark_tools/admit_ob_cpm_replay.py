"""Admit complete unscored CPM replay panel before candidate expansion."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
from benchmark_tools.prepare_ob_candidate_neighborhood import PREPARED_SHA, check, record
from benchmark_tools.run_ob_cpm_neighborhood import ARMS, RUNTIME_SHA, replay_command
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.build_publication_runtime import verify_runtime
from benchmark_tools.run_publication_replay_check import OUTPUTS, compare_stages
from benchmark_tools.audit_historical_profile_ablation import read_partition
from benchmark_tools.verify_ygob_validation import require_completed_job


def check_metrics(metrics, row, cache, source, frozen):
    expected = {"accuracy_profile": "high_sensitivity", "cpm_resolution": row["resolution"],
        "jackknife_profile_thresholds": False, "jackknife_single_copy_profiles": False,
        "leiden_seed": 4, "matrix": "BLOSUM62", "profile_expansion": True,
        "profile_iterations": 1, "profile_min_species": 1}
    if (metrics["parameters"] != expected or metrics["command"] != row["command"]
            or metrics["input"] != cache or metrics["source"] != source or metrics["cwd"] != str(frozen)
            or metrics["git"] != {"commit": "7f3a9e40dd7e79f842cc2c11fb8b548f9a802806", "dirty": False}):
        raise ValueError("Replay metadata differs from fixed design")
    counts = metrics["counts"]
    if any(counts[key] != value for key, value in (("genes", 251378), ("species", 12), ("significant_hits", 18235373))):
        raise ValueError("Replay input counts differ")
    iterations = metrics["profile_iterations"]
    if len(iterations) != 1 or iterations[0]["iteration"] != 1 or iterations[0]["profiles_built"] <= 0:
        raise ValueError("Missing positive HMM profile-build evidence")
    for key in ("profiles_built", "profile_candidates", "significant_profile_hits", "strict_profile_edges", "calibrated_profiles"):
        if iterations[0][key] != counts[key]:
            raise ValueError("Profile iteration accounting differs")
    if [stage["label"] for stage in metrics["stages"]] != list(OUTPUTS):
        raise ValueError("Missing or reordered native stages")
    if any(set(stage) != {"label", "clusters", "output"} for stage in metrics["stages"]):
        raise ValueError("Unexpected stage fields or benchmark scoring")


def admit(root, output):
    if output.exists():
        raise FileExistsError(output)
    accounting = subprocess.check_output(["sacct", "-j", "21319", "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = require_completed_job(accounting, 21319)
    directory = root / "benchmarks/results/ob_cpm_neighborhood_v1"
    path = directory / "results.json"
    report = json.loads(path.read_text())
    if (report["status"] != "three_cpm_arms_replayed_unscored" or report["accuracy_evaluated"] is not False
            or report["job_id"] != "21319" or [row["label"] for row in report["arms"]] != [name for name, _ in ARMS]):
        raise ValueError("CPM panel incomplete or changed")
    executor = root / "benchmarks/work/publication_ob_cpm_neighborhood_v1"
    revision = subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip()
    expected_revision = subprocess.check_output(["git", "-C", str(root), "rev-parse", "1bae2d2^{commit}"], text=True).strip()
    if revision != expected_revision or report["source"] != record(executor / "benchmark_tools/run_ob_cpm_neighborhood.py"):
        raise ValueError("Changed replay executor")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    prepared_path = root / "benchmark_tools/results/orthobench_factorial_prepared_20260916.json"
    prepared = read_frozen(prepared_path, PREPARED_SHA)
    check(prepared["replay_verification"])
    baseline = json.loads(Path(prepared["replay_verification"]["path"]).read_text())
    check(baseline["preflight"])
    preflight = json.loads(Path(baseline["preflight"]["path"]).read_text())
    frozen = root / "benchmarks/work/publication_method_native_v2"
    runtime = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    read_frozen(runtime, RUNTIME_SHA)
    verify_runtime(runtime, frozen)
    expected = {name: baseline["stages"][name]["output"] for name in OUTPUTS}
    inputs = [record(prepared_path), record(runtime), prepared["replay_verification"], baseline["preflight"],
        preflight["replay_source"], *prepared["core_sources"], *prepared["fasta_inputs"], prepared["cache"], *expected.values()]
    if report["inputs"] != inputs or report["environment_overrides"] != {**prepared["environment_overrides"], "PYTHONPATH": str(frozen)}:
        raise ValueError("Changed replay inputs or environment")
    records = [record(path), report["source"], *inputs]
    for item in records:
        check(item)
    universe = set()
    for item in prepared["fasta_inputs"]:
        for sequence in SeqIO.parse(item["path"], "fasta"):
            if sequence.id in universe:
                raise ValueError("Duplicate FASTA identifier")
            universe.add(sequence.id)
    for row, (label, resolution) in zip(report["arms"], ARMS):
        if (row["status"] != "replayed_unscored" or row["exit_code"] != 0 or row["resolution"] != float(resolution)
                or row["command"] != replay_command(preflight["command"], label, directory / label)):
            raise ValueError("Changed or unsuccessful replay arm")
        check(row["metrics"])
        records.append(row["metrics"])
        metrics = json.loads(Path(row["metrics"]["path"]).read_text())
        check_metrics(metrics, row, prepared["cache"], preflight["replay_source"], frozen)
        for stage in metrics["stages"]:
            name = stage["label"]
            expected_path = directory / label / "replay" / OUTPUTS[name]
            check(stage["output"])
            groups = read_partition(expected_path, universe)
            if (stage["output"] != record(expected_path) or stage["clusters"] != len(groups)
                    or row["stages"][name] != {"groups": len(groups), "output": stage["output"]}):
                raise ValueError("Native stage inventory/count/hash mismatch")
            records.append(stage["output"])
        if label == "control":
            comparison = compare_stages(expected, directory / label / "replay", universe)
            if comparison != row["control_equivalence"] or not all(r["byte_equal"] and r["partition_equal"] for r in comparison.values()):
                raise ValueError("Unchanged control equivalence failed")
    for item in records:
        check(item)
    verify_runtime(runtime, frozen)
    result = {"status": "cpm_replay_panel_verified_unscored", "accuracy_evaluated": False,
        "scheduler": scheduler, "arms": report["arms"], "source": record(__file__),
        "provenance_checked": records, "limitations": ["Native grouping/profile replay verified; candidate expansion and inferred phylogeny still pending.",
            "HMM build evidence is from hash-verified native metrics/runtime, not independent reproduction of every profile score.",
            "Shared-node incremental costs; no accuracy scoring or default selection."]}
    output.mkdir(parents=True)
    (output / "results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    admit(args.root.resolve(), args.output.resolve())
