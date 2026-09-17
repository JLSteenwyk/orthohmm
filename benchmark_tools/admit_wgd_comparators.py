"""Validate native WGD comparator outputs before biological endpoint scoring."""

import argparse
import json
from pathlib import Path
import shlex
import subprocess
import sys

from Bio import Phylo

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_wgd_orthohmm import check_receipt
from benchmark_tools.read_wgd_native_groups import read_species_table, read_orthofinder_root_ids, read_orthohmm
from benchmark_tools.run_wgd_application import DIRECTORIES, check_copies, pinned
from benchmark_tools.score_wgd_application import membership
from benchmark_tools.snapshot_orthohmm_input_order import record
from benchmark_tools.validate_scaling_outputs import input_universe, unique_path, validate_orthofinder


def check_root_tree(path, species):
    tree = Phylo.read(path, "newick")
    leaves = [leaf.name for leaf in tree.get_terminals()]
    if tree.root.name != "N0" or len(leaves) != len(set(leaves)) or set(leaves) != set(species):
        raise ValueError("Root-HOG species-tree scope differs from four-species inputs")


def check_sonic_log(text, selected, species_count):
    lines = text.splitlines()
    expected = ["SonicParanoid 2.0.9 will be executed with the following parameters:",
                "Input directory: " + selected["copy_inputs_to"],
                f"Input proteomes:\t{species_count}", "Threads:\t32",
                "Alignment tool:\tdiamond", "Run mode:\tdefault (Diamond [--very-sensitive])",
                "MCL inflation:\t1.50", "Perfom only graph-based orthology:\tFalse"]
    if any(lines.count(value) != 1 for value in expected):
        raise ValueError("SonicParanoid native configuration differs or is ambiguous")
    # A downstream native stage repeats this path; conflicting repetitions fail.
    output_lines = [line for line in lines if line.startswith("Main output directory: ")]
    if not output_lines or any(line != "Main output directory: " + selected["argv"][4] for line in output_lines):
        raise ValueError("SonicParanoid native output directory differs")
    totals = [line.split("\t")[-1] for line in lines if line.startswith("Total elapsed time (seconds):\t")]
    if len(totals) != 1:
        raise ValueError("SonicParanoid completion marker missing or ambiguous")
    import math
    elapsed = float(totals[0])
    if not math.isfinite(elapsed) or elapsed < 0:
        raise ValueError("Invalid SonicParanoid native elapsed time")


def admit(repo, index, accounting):
    if index not in (2, 3):
        raise ValueError("Comparator index must be2or3")
    spec_path = repo / "benchmark_tools/results/biological_wgd_execution_20260917.json"
    spec = pinned({"path": str(spec_path), "sha256": "c43704020c56ead316678461f5be3e8d4efc43a56ff5ced4f0e3cfd3c189025c"})
    plan = pinned(spec["command_plan"])
    selected = plan["runs"][index]
    root = Path(plan["output_root"]) / DIRECTORIES[index]
    receipt_path = root / "execution.json"
    receipt = json.loads(receipt_path.read_text())
    check_receipt(receipt, selected, root, spec, accounting, f"21661_{index}")
    time_path = root / "native.time.log"
    lines = time_path.read_text().splitlines()
    timed = [shlex.split(line.strip().removeprefix('Command being timed: "').removesuffix('"'))
             for line in lines if line.strip().startswith("Command being timed:")]
    if timed != [selected["argv"]] or [line.strip() for line in lines if "Exit status:" in line] != ["Exit status: 0"]:
        raise ValueError("GNU-time command/exit status differs")
    inputs = pinned(plan["inputs"])
    owners, species = input_universe({**inputs, "proteomes": 4})
    copies = check_copies(Path(selected["copy_inputs_to"]), inputs["inputs"])
    if index == 2:
        native = validate_orthofinder({"native_argv": selected["argv"], "dataset": inputs,
                                      "configuration": {"output": str(root), "copy_inputs_to": selected["copy_inputs_to"]}}, owners, species)
        path = unique_path(root, "**/WorkingDirectory/N0.ids.tsv")
        ids_path = unique_path(root, "**/WorkingDirectory/SequenceIDs.txt")
        tree_path = unique_path(root, "**/Species_Tree/SpeciesTree_rooted_node_labels.txt")
        check_root_tree(tree_path, species)
        groups = read_orthofinder_root_ids(path, ids_path, owners, {s: s for s in species})
        final_path = unique_path(root, "**/Orthogroups/Orthogroups.txt")
        final_groups = read_orthohmm(final_path, "named_groups", owners)
        root_sets = {frozenset(genes) for genes in groups.values() if len(genes) > 1}
        final_sets = {frozenset(genes) for genes in final_groups.values() if len(genes) > 1}
        if root_sets != final_sets:
            raise ValueError("Restored root-HOG non-singletons differ from native final groups")
        native["root_restoration"] = "N0.ids.tsv restored exactly through SequenceIDs; no singleton supplementation"
        native["checked_files"] += [record(path), record(ids_path), record(tree_path), record(final_path)]
    else:
        log_path = root / "native.log"
        check_sonic_log(log_path.read_text(), selected, len(species))
        path = unique_path(root / "output", "runs/*/ortholog_groups/ortholog_groups.tsv")
        columns = {Path(r["path"]).name: Path(r["path"]).stem for r in inputs["inputs"]}
        groups = read_species_table(path, "sonicparanoid", owners, columns)
        native = {"input_genes": len(owners), "checked_files": [record(path), record(log_path)]}
    assigned = membership(groups, owners)
    native.update(application_groups=len(groups), assigned_genes=len(assigned),
                  unassigned_genes=sorted(set(owners) - assigned.keys()))
    return {"status": "comparator_application_native_outputs_admitted", "accuracy_evaluated": False,
            "method": selected["method"], "output_semantics": selected["output_semantics"],
            "group_table": record(path), "native": native, "copied_inputs": copies,
            "receipt": record(receipt_path), "time_log": record(time_path), "spec": record(spec_path),
            "validator": record(__file__), "accounting": accounting,
            "limitations": ["Native integrity admission, not biological endpoint evaluation.",
                            "Unassigned input proteins remain missing, not synthetic singleton groups.",
                            "Original shared-host timing is not controlled efficiency evidence."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    accounting = subprocess.check_output(["sacct", "-j", "21661", "--noheader", "--parsable2",
                                         "--format=JobIDRaw,JobID,State,ExitCode", "-X"], text=True)
    result = admit(args.repo.resolve(), args.index, accounting)
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
