"""Prepare and execute the fresh retained-order QfO phylogeny arm."""

import argparse
from dataclasses import asdict
import json
import os
from pathlib import Path
import shutil
import sys

if __name__ == "__main__":
    # Scientific imports precede repository helpers in isolated worker execution.
    import orthohmm
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.run_qfo_order_replay import record, save

PROTOCOL_SHA = "4067cd3ff24f2380e489364ed3f03303dbcdde7c518b8a5c73b56e26532b43fb"


def validate(path, sha):
    if record(path)["sha256"] != sha:
        raise ValueError("Changed phylogeny plan")
    plan = json.loads(path.read_text())
    if (plan["arm"] != "retained_order" or plan["attempts"] != 1
            or plan["checkpoint_reuse"] is not False or plan["accuracy_evaluated"] is not False):
        raise ValueError("Wrong fresh-arm scope")
    for item in plan["checked_records"]:
        if record(item["path"]) != item:
            raise ValueError("Changed pinned artifact: " + item["path"])
    return plan


def config_arguments(plan):
    return dict(mode="reconcile", species_tree_mode="infer", species_tree=None,
                aligner=plan["aligner"], tree_builder=plan["tree_builder"],
                root_duplication_rule="species_overlap", pair_orthology_rule="positive_paralogy",
                species_tree_rooting="min_variance")


def prepare(repo, directory):
    from benchmark_tools.readback_qfo_order_replay import admit
    if directory.exists():
        raise FileExistsError(directory)
    candidate_dir = repo / "benchmarks/work/qfo_order_replay_20260927"
    admitted = admit(repo, candidate_dir)
    candidate = json.loads((candidate_dir / "plan.json").read_text())
    tool_plan = repo / "benchmarks/work/publication_full_recovery_orthobench_20260926/plan.json"
    if record(tool_plan)["sha256"] != "57878c0712b3216376ea1d431b1d8e9e8e3ae8eb58e3bd3069e15cf7dddb7005":
        raise ValueError("Changed validated tool plan")
    tools = json.loads(tool_plan.read_text())
    command = tools["command"]
    aligner, tree_builder = (command[command.index(flag) + 1] for flag in ("--aligner", "--tree-builder"))
    prefix = Path(aligner).parent.parent
    tool_records = [r for r in tools["checked_records"] if Path(r["path"]).is_relative_to(prefix)
                    or r["path"] == tree_builder]
    if record(aligner) not in tool_records or record(tree_builder) not in tool_records:
        raise ValueError("Missing validated external tools")
    prepared_path = repo / "benchmarks/results/qfo_corrected_factorial_v1/manifest.json"
    prepared = json.loads(prepared_path.read_text())
    fastas = prepared["input_fastas"]
    protocol = repo / "benchmark_tools/results/QFO_ORDER_DEPENDENCY_TRACE_20260927.md"
    if record(protocol)["sha256"] != PROTOCOL_SHA:
        raise ValueError("Changed downstream protocol")
    historical = repo / "benchmarks/results/qfo_corrected_factorial_v1/cells/p1_c1_r1"
    inputs = candidate_dir / "retained_order/orthohmm_working_res"
    checked = [*candidate["checked_records"], *admitted["checked_records"], *fastas,
               *tool_records, record(tool_plan), record(protocol), record(__file__),
               record(historical.with_suffix(".json"))]
    # Record historical comparison products now, without reading benchmark scores.
    checked.extend(record(historical / "orthohmm_phylogeny" / name) for name in
        ("orthohmm_root_hogs.tsv", "orthohmm_pairwise_orthologs.tsv",
         "orthohmm_pairwise_orthologs_confidence.tsv", "provenance_manifest.json",
         "reconciliation_summary.json", "species_tree.rooted.nwk"))
    checked = list({r["path"]: r for r in checked}.values())
    plan = dict(repo=str(repo), directory=str(directory), python=candidate["python"],
        arm="retained_order", attempts=1, checkpoint_reuse=False, accuracy_evaluated=False,
        environment=tools["environment"], aligner=aligner, tree_builder=tree_builder,
        checked_records=checked, fastas=fastas, input_directory=str(Path(fastas[0]["path"]).parent),
        partition=record(inputs / "orthohmm_edges_clustered.txt"),
        constraints=record(inputs / "phylogeny_candidate_merges.json"),
        expected_genes=984137, expected_species=78, cpu=32,
        historical=str(historical), resources=dict(cpus=32, memory_gib=128, hours=6),
        publication_ready=False, timing_scope="shared-host descriptive only")
    for item in checked:
        if record(item["path"]) != item:
            raise ValueError("Changed preparation artifact")
    directory.mkdir(parents=True)
    save(directory / "candidate_admission.json", admitted)
    save(directory / "plan.json", plan)
    return record(directory / "plan.json")


def worker(path, sha):
    import orthohmm
    from orthohmm.phylogeny import PhylogenyConfig
    from orthohmm.phylogeny_pipeline import run_phylogeny_stage
    plan = validate(path, sha)
    if (sys.executable != plan["python"] or not Path(orthohmm.__file__).resolve().is_relative_to(
            Path(sys.prefix).resolve())):
        raise ValueError("Wrong installed interpreter")
    if any(os.environ.get(k) != v for k, v in plan["environment"].items()):
        raise ValueError("Changed native environment")
    if plan["expected_genes"] == 984137 and (os.environ.get("SLURM_CPUS_PER_TASK") != "32"
                                           or not os.environ.get("SLURM_JOB_ID")):
        raise ValueError("Require full-data scheduled allocation")
    target = Path(plan["directory"]) / "native"
    target.mkdir(exist_ok=False)
    started = dict(plan=record(path), source=record(__file__), executable=sys.executable,
        environment={k: os.environ[k] for k in plan["environment"]}, config=config_arguments(plan),
        job_id=os.environ.get("SLURM_JOB_ID"), attempts=1, checkpoint_reuse=False)
    save(target / "started.json", started)
    try:
        working = target / "inference/orthohmm_working_res"
        working.mkdir(parents=True)
        shutil.copy2(plan["partition"]["path"], working / "orthohmm_edges_clustered.txt")
        shutil.copy2(plan["constraints"]["path"], working / "phylogeny_candidate_merges.json")
        constraints = json.loads(Path(plan["constraints"]["path"]).read_text())
        files = [Path(r["path"]).name for r in plan["fastas"]]
        result = run_phylogeny_stage(plan["input_directory"], str(target / "inference"), files,
            PhylogenyConfig(**config_arguments(plan)), plan["cpu"], membership_constraints=constraints or None)
        summary = asdict(result)
        if (summary["checkpoint_hits"] != 0 or summary["remapped_checkpoint_hits"] != 0
                or summary["species_tree_checkpoint_hit"] is not False):
            raise ValueError("Unexpected checkpoint reuse")
        validate(path, sha)
        save(target / "complete.json", dict(status="fresh_phylogeny_complete_pending_readback",
            started=record(target / "started.json"), summary=summary, plan=record(path),
            outputs=[record(p) for p in sorted((target / "inference").rglob("*")) if p.is_file()],
            accuracy_evaluated=False))
    except BaseException as error:
        save(target / "failure.json", dict(type=type(error).__name__, error=str(error), retry=False))
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("prepare", "worker"))
    parser.add_argument("--repo", type=Path)
    parser.add_argument("--directory", type=Path)
    parser.add_argument("--plan", type=Path)
    parser.add_argument("--sha256")
    args = parser.parse_args()
    if args.mode == "prepare":
        print(json.dumps(prepare(args.repo.resolve(), args.directory.resolve()), indent=2))
    else:
        worker(args.plan.resolve(), args.sha256)
