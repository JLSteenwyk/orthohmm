"""Complete the missing Three Kingdoms candidate/membership factorial cells."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time

import numpy as np
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_accuracy_checkpoint import audit
from benchmark_tools.audit_historical_profile_ablation import read_partition
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.prepare_orthobench_factorial import prepare_partition
from benchmark_tools.run_membership_hillclimb import (
    PRODUCTS, compare_golden, score_three_kingdoms, validate_partition, write_json,
)
from benchmark_tools.run_search_recovery import host_activity, runtime_records
from benchmark_tools.summarize_matched_resources import verbose_time
from benchmark_tools.validate_profile_runtime import require_profile_runtime

PLAN_NAME = "hillclimb_three_kingdoms_supplement_20261012.json"


def validate_ownership(names, species, owners):
    if not names or len(set(names)) != len(names) or set(names) != set(owners):
        raise ValueError("Checkpoint/FASTA gene universes differ")
    if species.shape != (len(names),) or species.dtype != np.dtype("int32"):
        raise ValueError("Invalid species vector")
    classes = {}
    for gene, label in zip(names, species):
        if label < 0:
            raise ValueError("Negative species label")
        classes.setdefault(int(label), set()).add(owners[gene])
    if any(len(files) != 1 for files in classes.values()):
        raise ValueError("A checkpoint species mixes FASTA owners")
    mapping = {label: next(iter(files)) for label, files in classes.items()}
    if len(set(mapping.values())) != len(mapping) or set(mapping.values()) != set(owners.values()):
        raise ValueError("Species classes are not bijective with FASTA files")
    return {"genes": len(names), "species": len(mapping), "classes": mapping}


def candidate_control_matches(candidate, constraints, golden_candidate, golden_constraints):
    matches = {"candidate_partition": record(candidate)["sha256"] == golden_candidate["sha256"],
               "membership_constraints": record(constraints)["sha256"] == golden_constraints["sha256"]}
    if not all(matches.values()):
        raise ValueError("Unchanged candidate control failed exact retained-output check: " + str(matches))
    return matches


def replay_command(source, execution, folder, resolution):
    return [sys.executable, str(source / "benchmark_tools/replay_high_sensitivity.py"),
        "--accuracy-checkpoint", execution["checkpoint"],
        "--checkpoint-sha256", execution["plan"]["checkpoint_manifest_sha256"],
        "--output-directory", str(folder / "graph"), "--json", str(folder / "graph_metrics.json"),
        "--fasta-directory", execution["fasta_directory"], "--cpu", str(execution["cpu"]),
        "--matrix", "BLOSUM62", "--cpm-resolution", str(resolution),
        "--leiden-seed", "4", "--profile-iterations", "1", "--profile-min-species", "1"]


def phylogeny_command(source, execution, folder, candidates, constraints, rule, checkpoint):
    return [sys.executable, str(source / "benchmark_tools/replay_phylogeny.py"),
        "--fasta-directory", execution["fasta_directory"],
        "--candidate-clusters", str(candidates), "--membership-constraints", str(constraints),
        "--output-directory", str(folder / "native"), "--json", str(folder / "native_metrics.json"),
        "--cpu", str(execution["cpu"]), "--species-tree-mode", "infer",
        "--aligner", execution["aligner"], "--tree-builder", execution["tree_builder"],
        "--root-rule", "species_overlap", "--pair-rule", "positive_paralogy",
        "--species-tree-rooting", "min_variance", "--satellite-membership-rule", rule,
        "--checkpoint-source", str(checkpoint)]


def sources(source):
    paths = [*sorted((source / "benchmark_tools").glob("*.py")),
             source / "benchmark_tools/run_three_kingdoms_membership_supplement.slurm",
             *sorted((source / "orthohmm").rglob("*.py")),
             *sorted((source / "orthohmm/search/csrc").glob("*.c")),
             *sorted((source / "orthohmm/search/csrc").glob("*.so"))]
    for name in ("pair_align.so", "hmm_viterbi.so", "kmer_prefilter.so"):
        if not (source / "orthohmm/search/csrc" / name).is_file():
            raise ValueError("Missing CPU extension: " + name)
    return [record(path) for path in paths]


def prepare(repo, source, output, cpu):
    if output.exists() or not 1 <= cpu <= 16:
        raise ValueError("Require fresh output and 1-16 CPUs")
    subprocess.run(["git", "-C", str(source), "diff", "--exit-code", "HEAD", "--",
                    "orthohmm", "benchmark_tools"], check=True)
    plan_path = source / "benchmark_tools" / PLAN_NAME
    plan = json.loads(plan_path.read_text())
    metadata_path = repo / plan["baseline_metadata"]
    metadata = json.loads(metadata_path.read_text())["metadata"]
    if (metadata["accuracy_profile"] != "high_sensitivity" or metadata["cpm_resolution"] != .1
            or metadata["phylogeny_candidates"] != "satellite_v2"):
        raise ValueError("Wrong Three Kingdoms baseline settings")
    checkpoint = Path(metadata["high_sensitivity_checkpoint"])
    if record(checkpoint / "manifest.json")["sha256"] != plan["checkpoint_manifest_sha256"]:
        raise ValueError("Wrong native HMM hit checkpoint")
    baseline = Path(metadata["output_directory"])
    working = baseline / "orthohmm_working_res"
    golden = {name: record(baseline / "orthohmm_phylogeny" / name)
              for name in (*PRODUCTS, "species_tree.rooted.nwk")}
    candidate = record(working / "phylogeny_candidate_superfamilies.txt")
    constraints = record(working / "phylogeny_candidate_merges.json")
    fasta = Path(metadata["fasta_directory"])
    fastas = [record(path) for path in sorted(fasta.iterdir()) if path.is_file()]
    if len(fastas) != 12 or any(Path(item["path"]).suffix != ".fasta" for item in fastas):
        raise ValueError("Wrong Three Kingdoms FASTA inventory")
    tools = []
    for name in ("aligner", "tree_builder"):
        resolved = shutil.which(metadata[name])
        if resolved is None:
            raise FileNotFoundError(metadata[name])
        tools.append(record(resolved))
    reference = repo / "three_kingdoms/busco/reference_orthogroups.txt"
    inputs = [record(plan_path), record(metadata_path), candidate, constraints, *golden.values(),
        *tools, *fastas, record(reference), record(repo / "three_kingdoms/score_against_busco.py"),
        *[record(path) for path in sorted(checkpoint.iterdir())]]
    for item in inputs:
        check(item)
    report = {"status": "prepared_unrun", "repo": str(repo), "source_root": str(source),
        "source_revision": subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip(),
        "plan": plan, "output": str(output), "cpu": cpu, "checkpoint": str(checkpoint),
        "baseline_output": str(baseline), "baseline_candidate": candidate,
        "baseline_constraints": constraints, "baseline_candidate_parameters": metadata["phylogeny_candidate_profile"]["parameters"],
        "golden": golden, "fastas": fastas, "fasta_directory": str(fasta),
        "aligner": metadata["aligner"], "tree_builder": metadata["tree_builder"],
        "inputs": inputs, "sources": sources(source), "runtime": runtime_records(),
        "profile_runtime": require_profile_runtime(source), "tool_path": os.environ["PATH"],
        "defaults_promoted": False, "accuracy_evaluated": False}
    output.mkdir(parents=True)
    write_json(output / "execution_plan.json", report)
    return report


def timed_stage(command, folder, name, source, env):
    timing, log = folder / (name + ".time.txt"), folder / (name + ".log")
    with log.open("x") as stream:
        result = subprocess.run(["/usr/bin/time", "-v", "-o", str(timing), *command],
                                cwd=source, env=env, stdout=stream, stderr=subprocess.STDOUT)
    receipt = {"command": command, "returncode": result.returncode,
               "log": record(log), "time": record(timing)}
    if not result.returncode:
        receipt["resources"] = verbose_time(timing.read_text())
    return receipt


def run(execution_path):
    from orthohmm.accuracy import load_accuracy_checkpoint
    from orthohmm.orthohmm import _expand_phylogeny_candidates
    from benchmark_tools.replay_phylogeny import load_membership_constraints

    execution_identity = record(execution_path)
    execution = json.loads(execution_path.read_text())
    source, output, repo = (Path(execution[key]) for key in ("source_root", "output", "repo"))
    if Path(__file__).resolve() != source / "benchmark_tools/run_three_kingdoms_membership_supplement.py":
        raise ValueError("Execute frozen supplement source")
    if Path(sys.modules["orthohmm.orthohmm"].__file__).resolve() != source / "orthohmm/orthohmm.py":
        raise ValueError("Wrong candidate engine")
    if (output / "results.json").exists():
        raise FileExistsError("An existing attempt must be preserved, not restarted in place")
    checked = [execution_identity, *execution["inputs"], *execution["sources"], *execution["runtime"]["files"]]
    for item in checked:
        check(item)
    if runtime_records() != execution["runtime"] or os.environ.get("CUDA_VISIBLE_DEVICES") != "":
        raise ValueError("Changed runtime or CUDA visibility")
    env = {**os.environ, "PATH": execution["tool_path"], "PYTHONPATH": str(source),
           "PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    report = {"status": "running", "execution_plan": execution_identity, "job_id": os.environ.get("SLURM_JOB_ID"),
              "source_revision": execution["source_revision"], "arms": [], "defaults_promoted": False,
              "scope": execution["plan"]["limitations"]}
    write_json(output / "results.json", report)
    control_valid, broader_valid = False, False
    try:
        report["host_activity_before"] = host_activity()
        checkpoint = Path(execution["checkpoint"])
        report["checkpoint_audit"] = audit(checkpoint, execution["plan"]["checkpoint_manifest_sha256"])
        names, species, queries, targets, scores = load_accuracy_checkpoint(checkpoint, verify=False)
        owners = {}
        for item in execution["fastas"]:
            for sequence in SeqIO.parse(item["path"], "fasta"):
                if sequence.id in owners:
                    raise ValueError("Duplicate FASTA gene")
                owners[sequence.id] = Path(item["path"]).name
        ownership = validate_ownership(names, species, owners)
        if ownership["genes"] != 443217 or ownership["species"] != 12:
            raise ValueError("Wrong Three Kingdoms sequence universe")
        report["checkpoint_fasta_ownership"] = ownership
        universe = set(names)
        for arm in execution["plan"]["arms"]:
            row = {**arm, "status": "running", "scores": {}}
            report["arms"].append(row)
            write_json(output / "results.json", report)
            label = arm["label"]
            if label != "control_direct" and not control_valid:
                row.update(status="not_run_failed_control", reason="Exact control was not verified")
                write_json(output / "results.json", report)
                continue
            if label == "broader_path" and not broader_valid:
                row.update(status="not_run_failed_broader_direct", reason="Matching broader direct arm failed")
                write_json(output / "results.json", report)
                continue
            folder = output / label
            folder.mkdir(exist_ok=False)
            started = time.monotonic()
            try:
                if label != "broader_path":
                    row["graph"] = timed_stage(replay_command(source, execution, folder, arm["cpm_resolution"]),
                                                folder, "graph", source, env)
                    if row["graph"]["returncode"]:
                        raise RuntimeError("Graph/profile replay failed")
                    graph = json.loads((folder / "graph_metrics.json").read_text())
                    for stage in graph["stages"]:
                        check(stage["output"])
                        read_partition(Path(stage["output"]["path"]), universe)
                    seed = folder / "graph/orthogroups_profiles_refined.txt"
                    usage, expansion_started = resource.getrusage(resource.RUSAGE_SELF), time.monotonic()
                    candidate_root = folder / "candidate"
                    row["candidate"] = prepare_partition(seed, candidate_root, names, species,
                        (queries, targets, scores), True, _expand_phylogeny_candidates, load_membership_constraints)
                    after = resource.getrusage(resource.RUSAGE_SELF)
                    row["candidate_resources"] = {"wall_s": time.monotonic() - expansion_started,
                        "cpu_s": after.ru_utime + after.ru_stime - usage.ru_utime - usage.ru_stime,
                        "parent_peak_rss_GiB": after.ru_maxrss / 1024 ** 2}
                    if row["candidate"]["expansion"]["parameters"] != execution["baseline_candidate_parameters"]:
                        raise ValueError("Satellite expansion parameters changed")
                    candidates = Path(row["candidate"]["candidate_partition"]["path"])
                    constraints = Path(row["candidate"]["membership_constraints"]["path"])
                    if label == "control_direct":
                        row["candidate_golden_matches"] = candidate_control_matches(candidates, constraints,
                            execution["baseline_candidate"], execution["baseline_constraints"])
                    checkpoint_source = Path(execution["baseline_output"])
                else:
                    direct = next(item for item in report["arms"] if item["label"] == "broader_direct")
                    row["candidate"] = direct["candidate"]
                    candidates = Path(row["candidate"]["candidate_partition"]["path"])
                    constraints = Path(row["candidate"]["membership_constraints"]["path"])
                    checkpoint_source = output / "broader_direct/native"
                row["native"] = timed_stage(phylogeny_command(source, execution, folder, candidates,
                    constraints, arm["membership_rule"], checkpoint_source), folder, "native", source, env)
                if row["native"]["returncode"]:
                    raise RuntimeError("Native phylogenetic inference failed")
                phylogeny = folder / "native/orthohmm_phylogeny"
                prediction = folder / "native/orthohmm_working_res/orthohmm_edges_clustered.txt"
                row["coverage"] = validate_partition(candidates, prediction)
                if label == "control_direct":
                    row["golden_matches"] = compare_golden(phylogeny, execution["golden"])
                elif label == "broader_path":
                    previous_tree = checkpoint_source / "orthohmm_phylogeny/species_tree.rooted.nwk"
                    if record(previous_tree)["sha256"] != record(phylogeny / "species_tree.rooted.nwk")["sha256"]:
                        raise ValueError("Membership-only change altered species tree")
                row["outputs"] = [record(phylogeny / name) for name in (*PRODUCTS, "species_tree.rooted.nwk")]
                row["native_metrics"] = record(folder / "native_metrics.json")
                row["native_summary"] = json.loads((phylogeny / "reconciliation_summary.json").read_text())
                row["membership_audit"] = json.loads((phylogeny / "provenance_manifest.json").read_text())["membership_reconciliation"]
                row["status"] = "native_outputs_verified_scoring_pending"
                write_json(output / "results.json", report)
                row.update(status="scored_development_result", **score_three_kingdoms(repo, prediction, folder))
                control_valid |= label == "control_direct"
                broader_valid |= label == "broader_direct"
            except Exception as error:
                row.update(status="failed", error_type=type(error).__name__, error=str(error))
            finally:
                row["total_incremental_wall_s"] = time.monotonic() - started
                for item in checked:
                    check(item)
                row["inputs_reverified_unchanged"] = True
                write_json(output / "results.json", report)
        report["host_activity_after"] = host_activity()
        if runtime_records() != execution["runtime"]:
            raise ValueError("Runtime changed during supplement")
        report["status"] = ("finished_all_arms_scored" if all(row["status"] == "scored_development_result"
                                                             for row in report["arms"]) else "finished_with_explicit_gaps")
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        for row in report["arms"]:
            if row["status"] == "scored_development_result":
                row["status"] = "invalidated_by_execution_failure"
        raise
    finally:
        write_json(output / "results.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    preparation = commands.add_parser("prepare")
    for name in ("repo", "source", "output"):
        preparation.add_argument("--" + name, type=Path, required=True)
    preparation.add_argument("--cpu", type=int, default=16)
    execution = commands.add_parser("run")
    execution.add_argument("--execution-plan", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "prepare":
        prepare(args.repo.resolve(), args.source.resolve(), args.output.resolve(), args.cpu)
    else:
        run(args.execution_plan.resolve())


if __name__ == "__main__":
    main()
