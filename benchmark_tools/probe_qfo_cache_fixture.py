"""Exercise filtered raw-tree reuse on the retained 16-gene fixture only."""

import argparse
from dataclasses import asdict
import json
import os
from pathlib import Path
import shutil
import sys


def run(repo, directory):
    from orthohmm import phylogeny_pipeline as pipeline
    from orthohmm.phylogeny import PhylogenyConfig
    if not Path(pipeline.__file__).resolve().is_relative_to(Path(sys.prefix).resolve()):
        raise ValueError("Require installed isolated implementation")
    sys.path.insert(0, str(repo))
    from benchmark_tools.run_qfo_fresh_phylogeny import config_arguments, record, save
    from benchmark_tools.qfo_phylogeny_cache import check, copy_cache, select_cache
    source = repo / "benchmarks/work/qfo_fresh_phylogeny_fixture_20260927"
    plan_path = source / "plan.json"
    plan = json.loads(plan_path.read_text())
    if plan["expected_genes"] != 16 or plan["expected_species"] != 4 or sys.executable != plan["python"]:
        raise ValueError("Wrong fixture")
    if any(os.environ.get(k) != v for k, v in plan["environment"].items()):
        raise ValueError("Wrong environment")
    checked = [record(plan_path), record(source / "native/complete.json"), record(__file__),
               record(repo / "benchmark_tools/qfo_phylogeny_cache.py"), *plan["checked_records"]]
    complete = json.loads((source / "native/complete.json").read_text())
    checked.extend(complete["outputs"])
    for item in checked:
        check(item)
    files = [Path(r["path"]).name for r in plan["fastas"]]
    sequences, gene_species = pipeline._load_sequence_data(plan["input_directory"], files)
    if len(sequences) != 16 or len(set(gene_species.values())) != 4:
        raise ValueError("Wrong fixture universe")
    genes = pipeline._read_clusters(Path(plan["partition"]["path"]))
    families = {f"Family{i:07d}": members for i, members in enumerate(genes)}
    config = PhylogenyConfig(**config_arguments(plan))
    tool_prefix = Path(plan["aligner"]).parent.parent
    tools = [r for r in plan["checked_records"] if Path(r["path"]).is_relative_to(tool_prefix)
             or r["path"] == plan["tree_builder"]]
    cache = select_cache(source / "native/inference/orthohmm_phylogeny", families, complete["outputs"],
        lambda name, members: pipeline._family_tree_input_hash(name, members, sequences, config), tools)
    directory.mkdir(parents=True, exist_ok=False)
    save(directory / "started.json", dict(checked_records=checked, cache=cache, config=config_arguments(plan)))
    try:
        working = directory / "inference/orthohmm_working_res"
        working.mkdir(parents=True)
        shutil.copy2(plan["partition"]["path"], working / "orthohmm_edges_clustered.txt")
        shutil.copy2(plan["constraints"]["path"], working / "phylogeny_candidate_merges.json")
        copies = copy_cache(cache, directory / "inference/orthohmm_phylogeny")
        result = asdict(pipeline.run_phylogeny_stage(plan["input_directory"], str(directory / "inference"),
            files, config, 2, membership_constraints=None))
        if (result["checkpoint_hits"] != 1 or result["remapped_checkpoint_hits"] != 0
                or result["species_tree_checkpoint_hit"] is not False):
            raise ValueError("Expected one raw-tree reuse and fresh species tree")
        for item in checked:
            check(item)
        save(directory / "complete.json", dict(status="cache_fixture_complete_pending_readback",
            summary=result, cache=cache, copies=copies, started=record(directory / "started.json"),
            outputs=[record(p) for p in sorted((directory / "inference").rglob("*")) if p.is_file()]))
    except BaseException as error:
        save(directory / "failure.json", dict(type=type(error).__name__, error=str(error), retry=False))
        raise
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--directory", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.repo.resolve(), args.directory.resolve()), indent=2))
