"""Construct fresh supplied-tree runs without changing frozen non-tree settings."""

from copy import deepcopy
from pathlib import Path


def _value(argv, flag):
    if argv.count(flag) != 1:
        raise ValueError("Expected exactly one " + flag)
    index = argv.index(flag) + 1
    if index == len(argv) or argv[index].startswith("-"):
        raise ValueError("Missing value for " + flag)
    return index


def fresh_supplied_method(method, original, tree, destination):
    """Return configuration only; callers must verify provenance and native outputs."""
    tree, destination = Path(tree).resolve(), Path(destination).resolve()
    old_output = Path(original["output"]).resolve()
    if (destination == old_output or destination in old_output.parents
            or old_output in destination.parents):
        raise ValueError("New and original output directories overlap")
    if tree == destination or destination in tree.parents:
        raise ValueError("Supplied tree must live outside the inference destination")
    if destination.exists():
        raise FileExistsError(destination)
    if not tree.is_file():
        raise FileNotFoundError(tree)
    result = deepcopy(original)
    argv = result["argv"]
    if method == "orthohmm_satellite_v2":
        if len(argv) < 5 or Path(argv[3]).resolve() != old_output or argv[4] != original["metrics"]:
            raise ValueError("Unexpected positional benchmark command")
        forbidden = {"--species-tree", "--checkpoint-source", "--species_tree", "--species_tree_mode"}
        if any(arg.split("=", 1)[0] in forbidden for arg in argv):
            raise ValueError("Baseline already supplies a tree or checkpoint")
        if argv[_value(argv, "--phylogeny")] != "reconcile":
            raise ValueError("Require phylogenetic reconciliation")
        mode = _value(argv, "--species-tree-mode")
        if argv[mode] != "infer":
            raise ValueError("Require the inferred-tree baseline")
        inputs = Path(argv[2]).resolve()
        metrics = destination.parent / (destination.name + ".json")
        if metrics.exists():
            raise FileExistsError(metrics)
        argv[3:5] = [str(destination), str(metrics)]
        argv[mode] = "supplied"
        argv.extend(["--species-tree", str(tree)])
        result["metrics"] = str(metrics)
    elif method == "orthofinder_full":
        forbidden = {"-s", "--speciestree", "-ft", "-fg", "-fst", "-b", "-o", "--assign", "--core"}
        if any(arg.split("=", 1)[0] in forbidden for arg in argv):
            raise ValueError("Expected fresh FASTA inference, not a tree or restart command")
        fasta = _value(argv, "-f")
        if argv[fasta] != original["copy_inputs_to"]:
            raise ValueError("FASTA argument and copied input directory disagree")
        inputs = Path(original["copy_inputs_from"]).resolve()
        copy = destination / "input"
        argv[fasta] = str(copy)
        argv.extend(["-s", str(tree)])
        result["copy_inputs_to"] = str(copy)
    else:
        raise ValueError("Unsupported phylogenetic method")
    if destination == inputs or destination in inputs.parents or inputs in destination.parents:
        raise ValueError("Input and output directories overlap")
    result.update(output=str(destination), supplied_tree=str(tree),
                  execution_scope="fresh supplied-tree inference; no upstream cache reuse")
    return result
