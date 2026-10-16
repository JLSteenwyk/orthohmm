"""Execute only original broad-panel refinement and candidate construction."""

import argparse
import json
from pathlib import Path
import shutil
import sys


class UnusedSearchHits:
    def values(self):
        raise RuntimeError("Broad-panel recovery must not access a search-hit mapping")


def verify_metadata(rows, names, species, minimum_species):
    import numpy as np
    if [str(row["name"]) for row in rows] != list(names):
        raise ValueError("FASTA metadata and checkpoint gene order differ")
    classes, observed = {}, []
    for row in rows:
        name = str(row["spp"])
        classes.setdefault(name, len(classes))
        observed.append(classes[name])
    if len(classes) < minimum_species or not np.array_equal(species, observed):
        raise ValueError("Require original broad-panel species ownership and order")
    return {"genes": len(names), "species": len(classes), "classes": classes}


def run(path):
    plan = json.loads(path.read_text())
    science, output = Path(plan["scientific_source"]), Path(plan["output"])
    if Path(__file__).resolve() != Path(plan["source_root"]) / "benchmark_tools/continue_cap400_downstream_worker.py":
        raise ValueError("Wrong frozen continuation worker")
    # Scientific imports must precede any executor-package imports.
    sys.path.insert(0, str(science))
    import numpy as np
    import orthohmm.orthohmm as engine
    from orthohmm import accuracy, helpers
    from orthohmm.files import fetch_fasta_files
    from benchmark_tools.audit_historical_profile_ablation import read_partition
    from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
    from benchmark_tools.replay_phylogeny import load_membership_constraints
    from benchmark_tools.run_membership_hillclimb import write_json

    modules = (engine, accuracy, helpers)
    if any(not Path(module.__file__).resolve().is_relative_to(science) for module in modules):
        raise ValueError("Nonoriginal scientific import")
    for pin in [*plan["checked"], *plan["sources"]]:
        check(pin)
    checkpoint = Path(plan["checkpoint"])
    names, species, queries, targets, scores = accuracy.load_accuracy_checkpoint(checkpoint, verify=True)
    files = fetch_fasta_files(plan["fasta"])
    rows = helpers.get_sequence_lengths(plan["fasta"], files)
    rows = rows[np.argsort(np.asarray([str(row["name"]) for row in rows]), kind="stable")]
    ownership = verify_metadata(rows, names, species, engine.DEFAULT_COPY_SPLIT_MIN_DATASET_SPECIES)
    if (ownership["genes"], ownership["species"]) != (984137, 78):
        raise ValueError("Wrong complete corrected QfO universe")
    payload = Path(plan["payload"])
    if (payload / "gene_names.txt").read_bytes() != (checkpoint / "gene_names.txt").read_bytes():
        raise ValueError("Recovered graph and checkpoint gene order differ")
    edges = helpers.IndexedEdges(names, *[np.load(payload / name, mmap_mode="r", allow_pickle=False)
        for name in ("sources.npy", "targets.npy", "weights.npy")])
    graph = output / "graph"
    working = graph / "orthohmm_working_res"
    working.mkdir(parents=True, exist_ok=False)
    partition = working / "orthohmm_edges_clustered.txt"
    shutil.copyfile(plan["partition"]["path"], partition)
    universe = set(names)
    read_partition(partition, universe)
    refinement = engine._refine_cluster_file(str(graph), rows, UnusedSearchHits(), edges,
                                            1e-4, refinement_profile="default")
    refined = graph / "orthogroups_refined.txt"
    shutil.copyfile(partition, refined)
    refined_groups = len(read_partition(refined, universe))
    candidate = engine._expand_phylogeny_candidates(str(graph), names, species,
        (queries, targets, scores), profile="satellite_v2")
    candidate.pop("_membership_constraints", None)
    candidates = working / "phylogeny_candidate_superfamilies.txt"
    constraints = working / "phylogeny_candidate_merges.json"
    candidate_groups = len(read_partition(candidates, universe))
    load_membership_constraints(constraints, candidates)
    products = [record(p) for p in (refined, candidates, constraints,
                working / "phylogeny_candidate_seeds.tsv", partition)]
    for pin in [*plan["checked"], *plan["sources"], *products]:
        check(pin)
    result = {"status": "original_refinement_and_satellite_candidates_verified",
        "execution_plan": record(path), "scientific_revision": plan["scientific_revision"],
        "modules": [record(module.__file__) for module in modules], "ownership": ownership,
        "refinement": refinement, "refined_groups": refined_groups,
        "candidate_groups": candidate_groups, "candidate": candidate, "products": products,
        "initial_hmm_search_rerun": False, "profile_expansion_rerun": False,
        "graph_clustering_rerun": False, "scores": {}}
    write_json(output / "graph_report.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execution-plan", required=True, type=Path)
    run(parser.parse_args().execution_plan.resolve())
