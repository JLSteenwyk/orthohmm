"""Run one QfO reconciliation cell after frozen candidate preparation succeeds."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_orthobench_factorial import plan_cells
from benchmark_tools.prepare_qfo_factorial import select_seeds, verify_inputs
from benchmark_tools.prepare_qfo_recovered_pairs import ADMISSION_SHA, stage_outputs
from benchmark_tools.run_simulation_methods import read_frozen, verify_environment, execution_environment, execute
from benchmark_tools.verify_qfo_replay_launcher import verify
from benchmark_tools.verify_ygob_validation import require_completed_job

ENVIRONMENT_SHA = "bf677728cf9312edd9755d9ac1ceefbce3abb8c3848d43413631dc00ca00eb2f"


def select_cell(manifest, index):
    if manifest.get("status") != "qfo_four_candidate_arms_prepared_unscored" or manifest.get("accuracy_computed") is not False:
        raise ValueError("Preparation is incomplete or scored")
    if set(manifest["candidate_arms"]) != {"p0_c0", "p0_c1", "p1_c0", "p1_c1"}:
        raise ValueError("Incomplete candidate arms")
    partition = Path(manifest["candidate_arms"]["p0_c0"]["candidate_partition"]["path"])
    output = partition.parents[3]
    parents = {Path(r["path"]).parent for r in manifest["input_fastas"]}
    if len(parents) != 1:
        raise ValueError("Ambiguous FASTA directory")
    prepared_executor = Path(manifest["source"]["path"]).parent.parent
    expected = plan_cells(output, next(iter(parents)), prepared_executor / "benchmark_tools/replay_phylogeny.py", 32)
    if manifest["cells"] != expected:
        raise ValueError("Cell commands differ from frozen design")
    cells = [r for r in expected if r["reconciliation"]]
    if type(index) is not int or not 0 <= index < 4:
        raise ValueError("Unknown reconciliation cell")
    return cells[index], output, prepared_executor


def native_command(cell, frozen_launcher, prepared_executor):
    original = cell["argv"]
    expected = prepared_executor / "benchmark_tools/replay_phylogeny.py"
    if original[:2] != [sys.executable, str(expected)]:
        raise ValueError("Unexpected planned interpreter or launcher")
    sources = []
    for name in ("replay_phylogeny.py", "orthobench_stage_diagnostics.py"):
        a, b = record(prepared_executor / "benchmark_tools" / name), record(frozen_launcher / "benchmark_tools" / name)
        if (a["sha256"], a["bytes"]) != (b["sha256"], b["bytes"]):
            raise ValueError("Frozen launcher or helper differs from planned source")
        sources.append({"prepared": a, "executed": b})
    return [original[0], str(frozen_launcher / "benchmark_tools/replay_phylogeny.py"), *original[2:]], sources


def verify_prepared(manifest, cell, job):
    if manifest["job_id"] != str(job):
        raise ValueError("Wrong preparation job")
    accounting = subprocess.check_output(["sacct", "-j", str(job), "--parsable2",
                                          "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = require_completed_job(accounting, job)
    for item in [manifest["source"], manifest["admission"], *manifest["helpers"], *manifest["input_fastas"]]:
        check(item)
    admission = read_frozen(Path(manifest["admission"]["path"]), ADMISSION_SHA)
    seeds = dict(select_seeds(stage_outputs(admission)))
    fasta = Path(manifest["input_fastas"][0]["path"]).parent
    verify_inputs([record(p) for p in sorted(fasta.glob("*.fasta"))], manifest["input_fastas"])
    for label, arm in manifest["candidate_arms"].items():
        profile, expanded = label.startswith("p1"), label.endswith("c1")
        if arm["seed_partition"] != seeds[profile] or arm["candidate_expansion"] is not expanded:
            raise ValueError("Candidate arm uses wrong seed or expansion setting")
        for key in ("seed_partition", "candidate_partition", "membership_constraints"):
            if key in arm:
                check(arm[key])
        if expanded != ("membership_constraints" in arm):
            raise ValueError("Wrong candidate membership policy")
    label = f"p{int(cell['profile_expansion'])}_c{int(cell['candidate_expansion'])}"
    arm = manifest["candidate_arms"][label]
    if cell["candidate_partition"] != arm["candidate_partition"]["path"]:
        raise ValueError("Cell candidate differs from prepared arm")
    if cell["candidate_expansion"]:
        argv = cell["argv"]
        if argv[argv.index("--membership-constraints") + 1] != arm["membership_constraints"]["path"]:
            raise ValueError("Cell constraints differ from prepared arm")
    runtime_path = Path(manifest["runtime_before"]["runtime_manifest"]["path"])
    current = verify(Path(manifest["core_root"]), Path(manifest["launcher_root"]), runtime_path)
    if current != manifest["runtime_before"] or current != manifest["runtime_after"]:
        raise ValueError("Runtime differs from successful candidate preparation")
    return scheduler


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--environment-manifest", type=Path, required=True)
    parser.add_argument("--preparation-job", type=int, required=True)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    manifest = read_frozen(args.manifest, args.manifest_sha256)
    cell, output, prepared_executor = select_cell(manifest, args.index)
    scheduler = verify_prepared(manifest, cell, args.preparation_job)
    launcher = Path(manifest["launcher_root"])
    argv, launcher_sources = native_command(cell, launcher, prepared_executor)
    environment = read_frozen(args.environment_manifest, ENVIRONMENT_SHA)
    verify_environment(environment)
    env, resolved = execution_environment(environment)
    env.update(manifest["environment_overrides"])
    env["PYTHONPATH"] = str(launcher)
    if args.check_only:
        print(f"Verified {cell['label']}; no inference or scoring")
        return
    evidence = output / "execution" / cell["label"]
    method = {"argv": argv, "output": argv[argv.index("--output-directory") + 1],
              "metrics": argv[argv.index("--json") + 1]}
    provenance = {"prepared_manifest": record(args.manifest), "preparation_scheduler": scheduler,
                  "environment": record(args.environment_manifest), "resolved_tools": resolved,
                  "planned_argv": cell["argv"], "executed_argv": argv, "launcher_source_equivalence": launcher_sources,
                  "executor": record(__file__), "execution_helper": record(Path(__file__).with_name("run_simulation_methods.py")),
                  "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
                  "cwd": str(launcher), "scope": "incremental cached reconciliation on shared workstation; no accuracy scoring"}
    inputs = {"status": "ready", "inputs": [{**r, "absolute_path": r["path"]} for r in manifest["input_fastas"]]}
    cwd = Path.cwd()
    try:
        os.chdir(launcher)
        result = execute({"label": cell["label"], "methods": {cell["label"]: method}}, [cell["label"]], env, evidence, inputs, provenance)
    finally:
        os.chdir(cwd)
    verify_prepared(manifest, cell, args.preparation_job)
    verify_environment(environment)
    for item in [provenance["executor"], provenance["execution_helper"], provenance["prepared_manifest"]]:
        check(item)
    for pair in launcher_sources:
        check(pair["prepared"])
        check(pair["executed"])
    with (evidence / "postflight.json").open("x") as stream:
        json.dump({"status": "inputs_runtime_sources_reverified", "failed_methods": result["failed_methods"],
                   "native_outputs_validated": False, "accuracy_evaluated": False}, stream, indent=2)
        stream.write("\n")
    if result["failed_methods"]:
        raise SystemExit("Reconciliation failed; outputs retained without scoring")
    print(f"{cell['label']}: process succeeded, native admission and scoring pending")


if __name__ == "__main__":
    main()
