"""Execute one frozen reconciliation cell without reading reference labels."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_historical_profile_ablation import verify_file
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.prepare_orthobench_factorial import plan_cells
from benchmark_tools.run_simulation_methods import read_frozen, verify_environment, execution_environment, execute
from benchmark_tools.verify_ygob_validation import require_completed_job


def select_cell(manifest, index):
    if manifest.get("status") != "prepared_not_reconciled" or manifest.get("accuracy_computed") is not False:
        raise ValueError("Expected label-blind prepared factorial")
    launcher = Path(manifest["launcher"]["path"])
    partition = Path(manifest["candidate_arms"]["p0_c0"]["candidate_partition"]["path"])
    root = partition.parents[3]
    parents = {Path(item["path"]).parent for item in manifest["fasta_inputs"]}
    if len(parents) != 1:
        raise ValueError("Expected one FASTA input directory")
    expected = plan_cells(root, next(iter(parents)), launcher, 32)
    if manifest["cells"] != expected:
        raise ValueError("Factorial commands differ from the prespecified design")
    cells = [cell for cell in expected if cell["reconciliation"]]
    if not 0 <= index < len(cells):
        raise ValueError("Unknown reconciliation cell")
    return cells[index], root, launcher.parent.parent


def verify_prepared(manifest, cell, launcher_root, preparation_job):
    accounting = subprocess.check_output(["sacct", "-j", str(preparation_job), "--parsable2",
                                          "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = require_completed_job(accounting, preparation_job)
    verify_file(Path(manifest["source"]["path"]), manifest["source"])
    for item in [manifest["launcher"], *manifest["helper_sources"], *manifest["fasta_inputs"]]:
        verify_file(Path(item["path"]), item)
    replay = manifest["replay_verification"]
    verify_file(Path(replay["path"]), replay)
    if json.loads(Path(replay["path"]).read_text())["status"] != "equivalent":
        raise ValueError("Replay equivalence is not established")
    roots = {Path(item["path"]).parent.parent for item in manifest["core_sources"]
             if Path(item["path"]).name == "orthohmm.py" and Path(item["path"]).parent.name == "orthohmm"}
    if len(roots) != 1:
        raise ValueError("Ambiguous prepared core source root")
    frozen_root = next(iter(roots))
    expected = set()
    for item in manifest["core_sources"]:
        source = Path(item["path"])
        verify_file(source, item)
        relative = source.relative_to(frozen_root)
        target = launcher_root / relative
        verify_file(target, item)
        expected.add(target)
    actual = {p for p in (launcher_root / "orthohmm").rglob("*") if p.is_file() and p.suffix in {".py", ".c", ".cu", ".h"}}
    if actual != expected:
        raise ValueError("Replay core file set differs from prepared source inventory")
    fasta = Path(manifest["fasta_inputs"][0]["path"]).parent
    actual_fastas = {p.resolve() for p in fasta.iterdir() if p.suffix.lower() in {".fa", ".faa", ".fasta", ".fsa"}}
    if actual_fastas != {Path(item["path"]).resolve() for item in manifest["fasta_inputs"]}:
        raise ValueError("FASTA input file set changed")
    label = f"p{int(cell['profile_expansion'])}_c{int(cell['candidate_expansion'])}"
    arm = manifest["candidate_arms"][label]
    for key in ("seed_partition", "candidate_partition", "membership_constraints"):
        if key in arm:
            verify_file(Path(arm[key]["path"]), arm[key])
    if cell["candidate_partition"] != arm["candidate_partition"]["path"]:
        raise ValueError("Cell does not use its own candidate arm")
    return scheduler


def unconstrained_cell(cell, output):
    """Prespecified satellite diagnostic: change only constraints and destinations."""
    if cell["label"] != "p1_c1_r1" or not all(cell[k] for k in
                                               ("profile_expansion", "candidate_expansion", "reconciliation")):
        raise ValueError("Unconstrained diagnostic requires the full satellite cell")
    argv = list(cell["argv"])
    if argv.count("--membership-constraints") != 1:
        raise ValueError("Expected exactly one membership constraint argument")
    index = argv.index("--membership-constraints")
    omitted = argv[index + 1]
    del argv[index:index + 2]
    argv.append("--unconstrained-membership")
    label = "p1_c1_r1_unconstrained_v2"
    target = output / "cells" / label
    argv[argv.index("--output-directory") + 1] = str(target)
    argv[argv.index("--json") + 1] = str(output / "cells" / (label + ".json"))
    return {**cell, "label": label, "argv": argv, "parent_cell": cell["label"],
            "omitted_membership_constraints": omitted,
            "prediction": str(target / "orthohmm_phylogeny/orthohmm_root_hogs.tsv")}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--environment-manifest", type=Path, required=True)
    parser.add_argument("--environment-sha256", required=True)
    parser.add_argument("--preparation-job", type=int, required=True)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--check-only", action="store_true")
    parser.add_argument("--unconstrained", action="store_true",
                        help="Separate prespecified p1_c1_r1 membership-filter diagnostic")
    args = parser.parse_args()
    manifest = read_frozen(args.manifest, args.manifest_sha256)
    cell, output, launcher_root = select_cell(manifest, args.index)
    scheduler = verify_prepared(manifest, cell, launcher_root, args.preparation_job)
    if args.unconstrained:
        cell = unconstrained_cell(cell, output)
    environment = read_frozen(args.environment_manifest, args.environment_sha256)
    verify_environment(environment)
    env, resolved = execution_environment(environment)
    env["PYTHONPATH"] = str(launcher_root)
    env.update(manifest["environment_overrides"])
    if args.check_only:
        print(f"Verified {cell['label']}: no inference or scoring")
        return
    argv = cell["argv"]
    target = Path(argv[argv.index("--output-directory") + 1])
    metrics = Path(argv[argv.index("--json") + 1])
    evidence = output / "execution" / cell["label"]
    method = {"argv": argv, "output": str(target), "metrics": str(metrics)}
    provenance = {"prepared_manifest": file_provenance(args.manifest), "preparation_scheduler": scheduler,
                  "environment_manifest": file_provenance(args.environment_manifest), "resolved_tools": resolved,
                  "launcher": manifest["launcher"], "executor": file_provenance(Path(__file__)),
                  "execution_helper": file_provenance(Path(__file__).with_name("run_simulation_methods.py")),
                  "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
                  "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
                  "runtime_kind": "incremental_cached_reconciliation; shared-machine timings",
                  "cwd": str(launcher_root),
                  "scope": "Reconciliation only; profile expansion is a verified upstream checkpoint"}
    if args.unconstrained:
        provenance.update(parent_cell=cell["parent_cell"],
                          omitted_membership_constraints=cell["omitted_membership_constraints"],
                          diagnostic="No membership filter; not an additional factorial cell")
    verified = {"status": "ready", "inputs": [{**r, "absolute_path": r["path"]} for r in manifest["fasta_inputs"]]}
    verification_cwd = Path.cwd()
    try:
        os.chdir(launcher_root)
        result = execute({"label": cell["label"], "methods": {cell["label"]: method}},
                         [cell["label"]], env, evidence, verified, provenance)
    finally:
        # Distribution discovery must use the same cwd before and after inference.
        os.chdir(verification_cwd)
    verify_prepared(manifest, cell, launcher_root, args.preparation_job)
    verify_environment(environment)
    if result.get("failed_methods"):
        raise SystemExit("Reconciliation failed; evidence preserved, no score assigned")
    print(f"{cell['label']}: process finished; native validation and scoring remain separate")


if __name__ == "__main__":
    main()
