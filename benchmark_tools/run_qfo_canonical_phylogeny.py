"""Run canonical QfO phylogeny only after fresh-arm reuse admission."""

import argparse
from dataclasses import asdict
import json
import os
from pathlib import Path
import shutil
import sys

if __name__ == "__main__":
    import orthohmm
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.run_qfo_order_replay import record, save
from benchmark_tools.run_qfo_fresh_phylogeny import config_arguments
from benchmark_tools.qfo_phylogeny_cache import select_cache, copy_cache


def validate(path, sha):
    if record(path)["sha256"] != sha:
        raise ValueError("Changed canonical plan")
    plan = json.loads(path.read_text())
    if (plan["arm"] != "canonical_order" or plan["attempts"] != 1
            or plan["reuse_policy"] != "admitted_input_identical_raw_trees"
            or plan["species_tree_cache_reuse"] is not False or plan["accuracy_evaluated"] is not False):
        raise ValueError("Wrong canonical scope")
    for item in plan["checked_records"]:
        if record(item["path"]) != item:
            raise ValueError("Changed canonical input: " + item["path"])
    return plan


def prepare(repo, directory):
    from benchmark_tools.admit_qfo_raw_tree_reuse import admit
    from benchmark_tools.readback_qfo_order_replay import admit as admit_candidates
    if directory.exists():
        raise FileExistsError(directory)
    fresh = repo / "benchmarks/work/qfo_fresh_phylogeny_20260927"
    gate = admit(repo, fresh)
    candidate_dir = repo / "benchmarks/work/qfo_order_replay_20260927"
    candidate_admission = admit_candidates(repo, candidate_dir)
    native = json.loads(Path(gate["native_plan"]["path"]).read_text())
    prefix = Path(native["aligner"]).parent.parent
    tools = [r for r in native["checked_records"] if Path(r["path"]).is_relative_to(prefix)
             or r["path"] == native["tree_builder"]]
    working = candidate_dir / "canonical_order/orthohmm_working_res"
    checked = [*gate["checked_records"], *candidate_admission["checked_records"], record(__file__),
               record(repo / "benchmark_tools/qfo_phylogeny_cache.py"),
               record(repo / "benchmark_tools/admit_qfo_raw_tree_reuse.py")]
    checked = list({r["path"]: r for r in checked}.values())
    plan = {key: native[key] for key in ("python", "environment", "fastas", "input_directory", "cpu",
            "expected_genes", "expected_species", "aligner", "tree_builder", "resources", "historical")}
    plan.update(repo=str(repo), directory=str(directory), arm="canonical_order", attempts=1,
        reuse_policy="admitted_input_identical_raw_trees", species_tree_cache_reuse=False,
        accuracy_evaluated=False, publication_ready=False, checked_records=checked, tool_records=tools,
        partition=record(working / "orthohmm_edges_clustered.txt"),
        constraints=record(working / "phylogeny_candidate_merges.json"))
    directory.mkdir(parents=True)
    save(directory / "reuse_admission.json", gate)
    plan["reuse_admission"] = record(directory / "reuse_admission.json")
    plan["checked_records"].append(plan["reuse_admission"])
    save(directory / "plan.json", plan)
    validate(directory / "plan.json", record(directory / "plan.json")["sha256"])
    return record(directory / "plan.json")


def worker(path, sha):
    import orthohmm
    from orthohmm import phylogeny_pipeline as pipeline
    from orthohmm.phylogeny import PhylogenyConfig
    plan = validate(path, sha)
    if (sys.executable != plan["python"] or not Path(orthohmm.__file__).resolve().is_relative_to(
            Path(sys.prefix).resolve())):
        raise ValueError("Wrong installed canonical interpreter")
    if any(os.environ.get(k) != v for k, v in plan["environment"].items()):
        raise ValueError("Changed canonical environment")
    if (os.environ.get("SLURM_CPUS_PER_TASK") != "32" or not os.environ.get("SLURM_JOB_ID")):
        raise ValueError("Canonical full run requires scheduled allocation")
    gate = json.loads(Path(plan["reuse_admission"]["path"]).read_text())
    if (gate["status"] != "fresh_qfo_raw_tree_reuse_admitted" or gate["native_job"] != 22329
            or gate["readback_job"] != 22332):
        raise ValueError("Wrong reuse admission")
    target = Path(plan["directory"]) / "native"
    target.mkdir(exist_ok=False)
    save(target / "started.json", dict(plan=record(path), source=record(__file__), executable=sys.executable,
        environment={k: os.environ[k] for k in plan["environment"]}, config=config_arguments(plan),
        reuse_admission=plan["reuse_admission"], job_id=os.environ["SLURM_JOB_ID"], attempts=1))
    try:
        working = target / "inference/orthohmm_working_res"
        working.mkdir(parents=True)
        shutil.copy2(plan["partition"]["path"], working / "orthohmm_edges_clustered.txt")
        shutil.copy2(plan["constraints"]["path"], working / "phylogeny_candidate_merges.json")
        files = [Path(r["path"]).name for r in plan["fastas"]]
        seqs, species = pipeline._load_sequence_data(plan["input_directory"], files)
        if len(seqs) != plan["expected_genes"] or len(set(species.values())) != plan["expected_species"]:
            raise ValueError("Wrong canonical universe")
        clusters = pipeline._read_clusters(Path(plan["partition"]["path"]))
        families = {f"Family{i:07d}": genes for i, genes in enumerate(clusters)}
        config = PhylogenyConfig(**config_arguments(plan))
        cache = select_cache(Path(gate["source_directory"]), families, gate["admitted_outputs"],
            lambda name, genes: pipeline._family_tree_input_hash(name, genes, seqs, config), plan["tool_records"])
        save(target / "cache_manifest.json", cache)
        copied = copy_cache(cache, target / "inference/orthohmm_phylogeny")
        save(target / "cache_copies.json", copied)
        # Release the selection copy before the native stage loads all sequences again.
        del seqs, species, families, clusters
        constraints = json.loads(Path(plan["constraints"]["path"]).read_text())
        summary = asdict(pipeline.run_phylogeny_stage(plan["input_directory"], str(target / "inference"),
            files, config, plan["cpu"], membership_constraints=constraints or None))
        if (summary["checkpoint_hits"] != len(cache["included_families"])
                or summary["remapped_checkpoint_hits"] != 0 or summary["species_tree_checkpoint_hit"] is not False):
            raise ValueError("Unexpected cache reuse or species-tree reuse")
        validate(path, sha)
        save(target / "complete.json", dict(status="canonical_phylogeny_complete_pending_readback",
            plan=record(path), started=record(target / "started.json"), summary=summary,
            cache=record(target / "cache_manifest.json"), copies=record(target / "cache_copies.json"),
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
