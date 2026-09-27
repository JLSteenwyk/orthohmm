"""Trace candidate-order changes before planning downstream QfO inference."""

import argparse
from collections import Counter
import json
from pathlib import Path
import sys


def family_records(path):
    groups = [tuple(sorted(line.split())) for line in path.read_text().splitlines() if line.split()]
    flat = [g for group in groups for g in group]
    if len(flat) != len(set(flat)):
        raise ValueError("Duplicate partition gene")
    return [(f"Family{i:07d}", genes) for i, genes in enumerate(groups)]


def trace_changes(left, right, left_constraints, right_constraints):
    left_by_genes = {genes: name for name, genes in left}
    right_by_genes = {genes: name for name, genes in right}
    common = left_by_genes.keys() & right_by_genes.keys()
    changed = dict(left=[dict(family=name, genes=list(genes)) for name, genes in left
                         if genes not in right_by_genes],
                   right=[dict(family=name, genes=list(genes)) for name, genes in right
                          if genes not in left_by_genes])
    semantic = lambda row: (tuple(sorted(row["source_genes"])), tuple(sorted(row["target_genes"])))
    lc, rc = Counter(map(semantic, left_constraints)), Counter(map(semantic, right_constraints))
    differences = {}
    for side, delta in (("left_only", lc - rc), ("right_only", rc - lc)):
        differences[side] = [dict(source_genes=list(source), target_genes=list(target), count=n)
                            for (source, target), n in sorted(delta.items())]
    return dict(changed_families=changed, constraint_multiset_difference=differences,
                shared_families=len(common), shared_families_with_changed_ids=sum(
                    left_by_genes[g] != right_by_genes[g] for g in common))


def run(repo, output):
    # Use the installed science package before exposing repository helpers.
    from orthohmm import phylogeny_pipeline as pipeline
    if not Path(pipeline.__file__).resolve().is_relative_to(Path(sys.prefix).resolve()):
        raise ValueError("Require isolated installed scientific package")
    sys.path.insert(0, str(repo))
    from benchmark_tools.run_qfo_order_replay import ARMS, record, save, validate_plan
    if output.exists():
        raise FileExistsError(output)
    directory = repo / "benchmarks/work/qfo_order_replay_20260927"
    plan_path = directory / "plan.json"
    sha = "bbbb8b04e10d72839dc51de5b1b0484e2bc61a6802a189616cb64f09872c5c6b"
    plan = validate_plan(plan_path, sha)
    if sys.executable != plan["python"]:
        raise ValueError("Wrong interpreter")
    source = record(pipeline.__file__)
    if source not in plan["checked_records"]:
        raise ValueError("Unpinned scientific selector")
    readback_path = repo / "benchmark_tools/results/qfo_order_replay_readback_22328.json"
    if record(readback_path)["sha256"] != "10fdcb758574438492f939efa53324ef0d0df98691d37ae38ed90ea9624f9dee":
        raise ValueError("Changed admitted readback")
    readback = json.loads(readback_path.read_text())
    checked = [record(readback_path), source, record(__file__), *readback["admission"]["checked_records"]]
    prepared_path = repo / "benchmarks/results/qfo_corrected_factorial_v1/manifest.json"
    if record(prepared_path) not in plan["checked_records"]:
        raise ValueError("Unpinned input preparation")
    prepared = json.loads(prepared_path.read_text())
    checked.extend(prepared["input_fastas"])
    for item in checked:
        if record(item["path"]) != item:
            raise ValueError("Changed trace input")
    fastas = [Path(r["path"]) for r in prepared["input_fastas"]]
    if len({p.parent for p in fastas}) != 1 or len(fastas) != 78:
        raise ValueError("Wrong FASTA layout")
    sequences, species = pipeline._load_sequence_data(str(fastas[0].parent), [p.name for p in fastas])
    if len(sequences) != plan["expected_genes"]:
        raise ValueError("Wrong sequence universe")
    families, constraints, markers = {}, {}, {}
    for arm in ARMS:
        working = directory / arm / "orthohmm_working_res"
        families[arm] = family_records(working / "orthohmm_edges_clustered.txt")
        if {g for _, genes in families[arm] for g in genes} != sequences.keys():
            raise ValueError("Partition universe differs")
        constraints[arm] = json.loads((working / "phylogeny_candidate_merges.json").read_text())
        markers[arm] = pipeline.select_species_tree_families(
            families[arm], species, sequences, sorted(set(species.values())), max_families=200)
    left, right = ARMS
    result = trace_changes(families[left], families[right], constraints[left], constraints[right])
    result["markers"] = dict(exact_records_equal=markers[left] == markers[right],
        ordered_gene_sets_equal=[g for _, g in markers[left]] == [g for _, g in markers[right]],
        arms={arm: [dict(family=name, genes=list(genes)) for name, genes in rows]
              for arm, rows in markers.items()})
    for side in ("left", "right"):
        for row in result["changed_families"][side]:
            row["requires_reconciliation"] = pipeline.family_requires_reconciliation(row["genes"], species)
            row["species"] = len({species[g] for g in row["genes"]})
    for item in checked:
        if record(item["path"]) != item:
            raise ValueError("Changed trace input during inspection")
    validate_plan(plan_path, sha)
    result.update(status="candidate_dependencies_traced", checked_records=checked,
        plan=record(plan_path), accuracy_evaluated=False, downstream_equivalence_proven=False,
        limitation="Marker input equality is not inferred-tree equality or full historical artifact admission")
    save(output, result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.repo.resolve(), args.output.resolve())
    print(json.dumps({k: v for k, v in result.items() if k not in ("markers", "checked_records")}, indent=2))
