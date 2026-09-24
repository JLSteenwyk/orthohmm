"""Evaluate the frozen copy-split predicate on retained VATB-containing groups."""

import argparse
import ast
from collections import Counter
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_simulation_methods import read_frozen

TRACE_SHA = "797d4e3f019948657379c8a261b35fe582f1ba9d65839afb164f92fa1f44ae4c"
REFINEMENT_SHA = "991f1eb6a5f73d0442529ed19095a34b7c6ba8bff8dfe43a1127e24ec73fb26d"
CHECKPOINT_SHA = "30d1718e0ad1ec080b9783b0f94f3e8855063459c4b2999d86f93eb8d1b4c8e1"


def frozen_predicate(source):
    tree = ast.parse(source)
    constants = {}
    predicate = None
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            name = node.targets[0].id
            if name.startswith("DEFAULT_") and isinstance(node.value, ast.Constant):
                constants[name] = ast.literal_eval(node.value)
        if isinstance(node, ast.FunctionDef) and node.name == "_should_split_large_copy_cluster":
            predicate = node
    if predicate is None:
        raise ValueError("Missing frozen copy guard")
    namespace = {"Counter": Counter}
    exec(compile(ast.Module(body=[predicate], type_ignores=[]), "frozen-copy-guard", "exec"), namespace)
    return constants, namespace[predicate.name]


def singleton_members(path, wanted):
    found, singletons = set(), set()
    with path.open() as stream:
        for line in stream:
            genes = line.split()
            overlap = wanted.intersection(genes)
            if overlap & found or len(genes) != len(set(genes)):
                raise ValueError("Duplicate partition member")
            found.update(overlap)
            if len(genes) == 1:
                singletons.update(overlap)
    if found != wanted:
        raise ValueError("Missing affected group members")
    return len(singletons)


def run(root, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    trace_path = root / "benchmark_tools/results/qfo_vatb_partition_trace_20260923.json"
    trace = read_frozen(trace_path, TRACE_SHA)
    source = root / "benchmarks/work/publication_qfo_replay_native_v1/orthohmm/refinement.py"
    if record(source)["sha256"] != REFINEMENT_SHA:
        raise ValueError("Changed frozen refinement")
    defaults, predicate = frozen_predicate(source.read_text())
    checkpoint = root / "benchmarks/results/qfo_corrected_primary_v1/orthohmm_high_sensitivity/orthohmm_working_res/high_sensitivity_checkpoint"
    manifest = read_frozen(checkpoint / "manifest.json", CHECKPOINT_SHA)
    records = [record(__file__), record(trace_path), record(source), record(checkpoint / "manifest.json")]
    for name in ("gene_names.txt", "gene_to_species.npy"):
        item = record(checkpoint / name)
        if {key: item[key] for key in ("bytes", "sha256")} != manifest["files"][name]:
            raise ValueError("Changed saved species mapping")
        records.append(item)
    names = (checkpoint / "gene_names.txt").read_text().splitlines()
    species = np.load(checkpoint / "gene_to_species.npy", allow_pickle=False)
    if species.shape != (984137,) or len(names) != 984137 or len(set(names)) != len(names):
        raise ValueError("Wrong mapping universe")
    mapping = dict(zip(names, map(int, species)))
    n_species = len(set(mapping.values()))
    if n_species != 78:
        raise ValueError("Wrong dataset species count")
    parameters = dict(min_size=defaults["DEFAULT_COPY_SPLIT_MIN_SIZE"],
        min_species_count=max(defaults["DEFAULT_COPY_SPLIT_MIN_SPECIES_COUNT"], defaults["DEFAULT_BROAD_COPY_SPLIT_MIN_SPECIES_COUNT"]),
        min_dataset_species=defaults["DEFAULT_COPY_SPLIT_MIN_DATASET_SPECIES"])
    rows = []
    for before_label, after_label in (("multipass", "multipass_refined"), ("strict_profiles", "strict_profiles_refined")):
        before = next(row for row in trace["stages"] if row["arm"] == "cpm_low" and row["stage"] == before_label)
        after = next(row for row in trace["stages"] if row["arm"] == "cpm_low" and row["stage"] == after_label)
        for item in (before["partition"], after["partition"]):
            check(item)
            records.append(item)
        selected = {g["line"]: g for g in before["groups"]}
        with Path(before["partition"]["path"]).open() as stream:
            for number, line in enumerate(stream, 1):
                if number not in selected:
                    continue
                genes = line.split()
                counts = Counter(mapping[gene] for gene in genes)
                decision = predicate(len(genes), counts, n_species, **parameters)
                singles = singleton_members(Path(after["partition"]["path"]), set(genes))
                rows.append(dict(stage=before_label, line=number, group_genes=len(genes),
                    reference_genes=len(selected[number]["reference_genes"]), species_counts=dict(sorted(counts.items())),
                    max_species_copies=max(counts.values()), predicate_result=decision,
                    observed_singleton_members=singles, predicted_all_singletons_match=(singles == len(genes)) if decision else None))
    for item in records:
        check(item)
    result = dict(status="vatb_frozen_copy_guard_evaluated", parameters=parameters, dataset_species=78,
        groups=rows, checked_records=records, publication_ready=False,
        limitations=["Post-hoc predicate evaluation plus full affected-group output check; not a dynamic branch execution trace.",
                     "No full refinement, search, clustering or phylogeny rerun; no defaults changed.",
                     "False predicate alone does not exclude other splitting rules or explain candidate expansion."])
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.root.resolve(), args.output.absolute())
