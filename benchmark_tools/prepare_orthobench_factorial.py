"""Build controlled OrthoBench candidates after verified cached-stage replay."""

import argparse
import importlib.metadata
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_historical_profile_ablation import read_partition, verify_file
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.verify_ygob_validation import require_completed_job


def plan_cells(output, fasta, launcher, cpu):
    cells = []
    for profile in (False, True):
        for expansion in (False, True):
            candidate = output / "candidates" / f"p{int(profile)}_c{int(expansion)}" / "orthohmm_working_res"
            for reconciliation in (False, True):
                label = f"p{int(profile)}_c{int(expansion)}_r{int(reconciliation)}"
                row = {"label": label, "profile_expansion": profile, "candidate_expansion": expansion,
                       "reconciliation": reconciliation, "candidate_partition": str(candidate / "orthohmm_edges_clustered.txt"),
                       "prediction_format": "root_hogs" if reconciliation else "space_separated_groups",
                       "runtime_kind": "incremental_cached_replay"}
                if reconciliation:
                    target = output / "cells" / label
                    argv = [sys.executable, str(launcher), "--fasta-directory", str(fasta),
                            "--candidate-clusters", row["candidate_partition"], "--output-directory", str(target),
                            "--json", str(output / "cells" / f"{label}.json"), "--cpu", str(cpu),
                            "--species-tree-mode", "infer", "--aligner", "mafft", "--tree-builder", "FastTree",
                            "--root-rule", "species_overlap", "--pair-rule", "positive_paralogy",
                            "--species-tree-rooting", "min_variance"]
                    if expansion:
                        argv.extend(["--membership-constraints", str(candidate / "phylogeny_candidate_merges.json")])
                    row.update(argv=argv, prediction=str(target / "orthohmm_phylogeny/orthohmm_root_hogs.tsv"))
                else:
                    row["prediction"] = row["candidate_partition"]
                cells.append(row)
    return cells


def indexed_species(payload, owners):
    names = sorted(str(g) for g in payload["all_gene_ids"])
    if len(names) != len(set(names)) or set(names) != set(owners):
        raise ValueError("Cache gene universe differs from validated FASTAs")
    labels = {str(payload["gene_to_species"][g]) for g in names}
    label_owners = {label: {owners[g] for g in names if str(payload["gene_to_species"][g]) == label} for label in labels}
    if any(len(v) != 1 for v in label_owners.values()) or len({next(iter(v)) for v in label_owners.values()}) != len(labels):
        raise ValueError("Cache species classes differ from FASTA ownership")
    indices = {label: index for index, label in enumerate(sorted(labels))}
    species = np.array([indices[str(payload["gene_to_species"][g])] for g in names], dtype=np.int32)
    return names, species


def prepare_partition(seed, target, names, species, hits, expand, engine, validate_constraints):
    universe = set(names)
    read_partition(seed, universe)
    working = target / "orthohmm_working_res"
    working.mkdir(parents=True, exist_ok=False)
    partition = working / "orthohmm_edges_clustered.txt"
    shutil.copy2(seed, partition)
    result = {"seed_partition": file_provenance(seed), "candidate_expansion": expand}
    if expand:
        details = engine(str(target), names, species, hits, profile="satellite_v2")
        details.pop("_membership_constraints", None)
        trace = working / "phylogeny_candidate_merges.json"
        validate_constraints(trace, partition)
        result.update(expansion=details, membership_constraints=file_provenance(trace))
    read_partition(partition, universe)
    result["candidate_partition"] = file_provenance(partition)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replay-root", type=Path, required=True)
    parser.add_argument("--replay-job", type=int, required=True)
    parser.add_argument("--frozen-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cpu", type=int, default=32)
    args = parser.parse_args()
    if args.cpu != 32 or args.output.exists():
        raise ValueError("Require frozen CPU budget and absent output directory")
    accounting = subprocess.check_output(["sacct", "-j", str(args.replay_job), "--parsable2",
                                         "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = require_completed_job(accounting, args.replay_job)
    replay = args.replay_root.resolve()
    verification = json.loads((replay / "verification.json").read_text())
    if verification["status"] != "equivalent" or verification["exit_code"] != 0:
        raise ValueError("Replay has not established stage equivalence")
    verify_file(replay / "preflight.json", verification["preflight"])
    preflight = json.loads((replay / "preflight.json").read_text())
    if preflight["job_id"] != str(args.replay_job) or preflight["accuracy_scoring_requested"]:
        raise ValueError("Unexpected replay job or scoring configuration")
    root = Path(__file__).resolve().parent.parent
    frozen = args.frozen_root.resolve()
    commit = subprocess.check_output(["git", "-C", str(frozen), "rev-parse", "HEAD"], text=True).strip()
    if not commit.startswith("7f3a9e4") or commit != preflight["source_commit"]:
        raise ValueError("Wrong frozen core revision")
    # Launcher worktree may have newer benchmark scripts, but the core must be byte-identical.
    expected_sources = {Path(item["path"]).relative_to(frozen) for item in preflight["source_manifest"]}
    actual_sources = {p.relative_to(root) for p in (root / "orthohmm").rglob("*")
                      if p.is_file() and p.suffix in {".py", ".c", ".cu", ".h"}}
    if expected_sources != actual_sources:
        raise ValueError("Launcher worktree core source file set differs")
    for item in preflight["source_manifest"]:
        source = Path(item["path"])
        verify_file(source, item)
        verify_file(root / source.relative_to(frozen), item)
    audit_path = Path(preflight["historical_audit"]["path"])
    verify_file(audit_path, preflight["historical_audit"])
    audit = json.loads(audit_path.read_text())
    owners = {}
    parents = set()
    for item in audit["inputs"]["fastas"]:
        path = Path(item["path"])
        verify_file(path, item)
        parents.add(path.parent)
        for record in SeqIO.parse(path, "fasta"):
            if record.id in owners:
                raise ValueError("Duplicate FASTA ID")
            owners[record.id] = path.name
    if len(parents) != 1:
        raise ValueError("Expected one validated FASTA directory")
    stages = verification["stages"]
    for item in stages.values():
        if not item["partition_equal"]:
            raise ValueError("Non-equivalent stage in replay")
        verify_file(Path(item["output"]["path"]), item["output"])
    cache = audit["inputs"]["cache"]
    verify_file(Path(cache["path"]), cache)
    output = args.output.resolve()
    cells = plan_cells(output, next(iter(parents)), root / "benchmark_tools/replay_phylogeny.py", args.cpu)
    output.mkdir(parents=True)
    report = {"schema_version": 1, "status": "preparing", "accuracy_computed": False,
              "scheduler": scheduler, "replay_verification": file_provenance(replay / "verification.json"),
              "core_commit": commit, "core_source_equivalence_verified": True, "cells": cells,
              "initial_search": "frozen HMM normalized-hit cache; profile-off is not HMM-free",
              "cache": cache, "candidate_arms": {}, "source": file_provenance(Path(__file__)),
              "fasta_inputs": audit["inputs"]["fastas"], "core_sources": preflight["source_manifest"],
              "launcher": file_provenance(root / "benchmark_tools/replay_phylogeny.py"),
              "environment_overrides": {"PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"},
              "packages": sorted((d.metadata["Name"], d.version) for d in importlib.metadata.distributions()),
              "remaining_controls": ["unconstrained membership replay", "matched sequence-search control", "QfO factorial"]}
    path = output / "manifest.json"
    try:
        path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        from benchmark_tools.replay_high_sensitivity import load_hits, hit_arrays
        from benchmark_tools.replay_phylogeny import load_membership_constraints
        from orthohmm.orthohmm import _expand_phylogeny_candidates
        if Path(sys.modules["orthohmm.orthohmm"].__file__).resolve() != root / "orthohmm/orthohmm.py":
            raise ValueError("Candidate engine loaded from an unverified core location")
        report["helper_sources"] = [file_provenance(root / "benchmark_tools" / name) for name in
                                     ("replay_high_sensitivity.py", "replay_phylogeny.py", "audit_historical_profile_ablation.py")]
        payload = load_hits(Path(cache["path"]))
        names, species = indexed_species(payload, owners)
        hits = hit_arrays(payload["all_hits"], {g: i for i, g in enumerate(names)})
        if not np.isfinite(hits[2]).all():
            raise ValueError("Nonfinite normalized search scores")
        del payload
        for profile, stage in ((False, "multipass_refined"), (True, "strict_profiles_refined")):
            seed = Path(stages[stage]["output"]["path"])
            for expansion in (False, True):
                label = f"p{int(profile)}_c{int(expansion)}"
                started = time.perf_counter()
                record = prepare_partition(seed, output / "candidates" / label, names, species, hits,
                                           expansion, _expand_phylogeny_candidates, load_membership_constraints)
                record["incremental_preparation_seconds"] = time.perf_counter() - started
                report["candidate_arms"][label] = record
                path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        report["status"] = "prepared_not_reconciled"
    except Exception as error:
        report.update(status="failed", error=str(error))
        raise
    finally:
        path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
