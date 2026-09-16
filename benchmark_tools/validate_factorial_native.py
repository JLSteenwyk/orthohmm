"""Join execution, native provenance and partition gates for one frozen cell."""

import argparse
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import Phylo, SeqIO
from benchmark_tools.recover_factorial_postflight import audit, PREPARED_HASH, ENVIRONMENT_HASH
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.run_orthobench_factorial_cell import select_cell
from benchmark_tools.replay_phylogeny import load_membership_constraints
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.validate_factorial_partition import validate_partition


def check_native_metadata(metrics, native, cell, constraint_count):
    explicit = bool(cell.get("omitted_membership_constraints"))
    expanded = cell["candidate_expansion"] and not explicit
    supplied = cell.get("species_tree")
    checkpoint = cell.get("checkpoint_source")
    expected = {"aligner": "mafft", "checkpoint_source": str(Path(checkpoint) / "orthohmm_phylogeny") if checkpoint else None, "cpu": 32,
                "explicit_unconstrained_ablation": explicit, "pair_orthology_rule": "positive_paralogy",
                "root_duplication_rule": "species_overlap", "species_tree": supplied,
                "species_tree_mode": "supplied" if supplied else "infer", "species_tree_rooting": "min_variance",
                "tree_builder": "FastTree",
                "satellite_membership_policy": "high_confidence_pair" if expanded else "unconstrained"}
    if metrics.get("status") != "complete" or metrics.get("parameters") != expected:
        raise ValueError("Native completion or parameters mismatch")
    if metrics.get("command") != cell["argv"]:
        raise ValueError("Native command differs from frozen cell")
    if native.get("mode") != "reconcile" or native.get("cpu_budget") != 32:
        raise ValueError("Native reconciliation mode or CPU mismatch")
    for key in ("pair_orthology_rule", "root_duplication_rule", "species_tree_mode", "species_tree_rooting"):
        if native.get(key) != expected[key]:
            raise ValueError(f"Native rule mismatch: {key}")
    expected_source = supplied["path"] if supplied else "internally_inferred_from_orthohmm_single_copy_families"
    if native.get("species_tree_source") != expected_source:
        raise ValueError("Species tree source differs from the declared cell")
    membership = native.get("membership_reconciliation")
    if expanded:
        if (not isinstance(membership, dict) or membership.get("policy") != "high_confidence_pair" or
                membership.get("constraints") != constraint_count or
                membership.get("supported_constraints", -1) + membership.get("detached_constraints", -1) != constraint_count):
            raise ValueError("Native membership accounting mismatch")
    elif membership is not None or constraint_count:
        raise ValueError("Unexpected membership filtering in unexpanded cell")


def validate(root, index):
    integrity = audit(root, index)
    results = root / "benchmark_tools/results"
    prepared = read_frozen(results / "orthobench_factorial_prepared_20260916.json", PREPARED_HASH)
    environment = read_frozen(results / "publication_variable_native_methods_20260916.json", ENVIRONMENT_HASH)
    cell, output, launcher = select_cell(prepared, index)
    return validate_native_cell(prepared, environment, cell, output, launcher, integrity)


def validate_native_cell(prepared, environment, cell, output, launcher, integrity):
    """Native semantics shared by factorial cells and the isolated diagnostic."""
    target = Path(cell["argv"][cell["argv"].index("--output-directory") + 1])
    native_dir = target / "orthohmm_phylogeny"
    metrics_path = Path(cell["argv"][cell["argv"].index("--json") + 1])
    metrics = json.loads(metrics_path.read_text())
    native_path = native_dir / "provenance_manifest.json"
    native = json.loads(native_path.read_text())
    summary_path = native_dir / "reconciliation_summary.json"
    summary = json.loads(summary_path.read_text())
    arm = prepared["candidate_arms"][f"p{int(cell['profile_expansion'])}_c{int(cell['candidate_expansion'])}"]
    candidate = Path(arm["candidate_partition"]["path"])
    constraint = None if cell.get("omitted_membership_constraints") else arm.get("membership_constraints")
    constraints = load_membership_constraints(Path(constraint["path"]), candidate) if constraint else []
    check_native_metadata(metrics, native, cell, len(constraints))
    if metrics["git"] != {"commit": "b66225dc5fc355702575aebe34b644916894f236", "dirty": False}:
        raise ValueError("Native replay source revision mismatch")
    if metrics["source"] != prepared["launcher"] or metrics["cwd"] != str(launcher):
        raise ValueError("Native launcher provenance mismatch")
    if metrics["input"]["candidate_clusters"] != file_provenance(candidate):
        raise ValueError("Native candidate provenance mismatch")
    if metrics["input"]["membership_constraints"] != constraint or metrics["metadata"]["membership_constraints"] != constraint:
        raise ValueError("Native constraint provenance mismatch")
    if native["input_cluster_sha256"] != file_provenance(candidate)["sha256"]:
        raise ValueError("Reconciliation candidate hash mismatch")
    expected_inputs = {Path(r["path"]).name: r["sha256"] for r in prepared["fasta_inputs"]}
    observed = {r["filename"]: r["sha256"] for r in native["input_proteomes"]}
    if observed != expected_inputs or len(native["input_proteomes"]) != len(expected_inputs):
        raise ValueError("Native proteome provenance mismatch")
    if sorted(metrics["input"]["files"]) != sorted(expected_inputs):
        raise ValueError("Replay proteome list mismatch")
    if metrics["counts"] != {k: v for k, v in summary.items() if k != "schema_version"} or native["results"] != summary:
        raise ValueError("Native completion records disagree")
    for key, path in {"manifest": native_path, "summary": summary_path,
                      "root_hogs": native_dir / "orthohmm_root_hogs.tsv"}.items():
        if metrics["outputs"][key] != file_provenance(path):
            raise ValueError("Native output provenance mismatch")
    tree_path = native_dir / "species_tree.rooted.nwk"
    if native["species_tree_sha256"] != file_provenance(tree_path)["sha256"]:
        raise ValueError("Species tree hash mismatch")
    tree = Phylo.read(tree_path, "newick")
    leaves = [leaf.name for leaf in tree.get_terminals()]
    taxa = [r["taxon"] for r in native["input_proteomes"]]
    if sorted(leaves) != sorted(taxa) or sorted(leaves) != native["species_tree_taxa"] or len(set(leaves)) != len(leaves):
        raise ValueError("Species tree taxon coverage mismatch")
    if any(n.branch_length is not None and not math.isfinite(n.branch_length) for n in tree.find_clades()):
        raise ValueError("Nonfinite species-tree branch")
    for role, tool, version in (("aligner", "mafft", "v7.525 (2024/Mar/13)"),
                                ("tree_builder", "FastTree", "FastTree Version 2.2.0 Double precision")):
        if native["tools"][role] != {"path": environment["tool_entrypoints"][tool]["absolute_path"], "version": version}:
            raise ValueError("Native external tool mismatch")
    genes = set()
    for item in prepared["fasta_inputs"]:
        for record in SeqIO.parse(item["path"], "fasta"):
            if record.id in genes:
                raise ValueError("Duplicate FASTA gene identifier")
            genes.add(record.id)
    _, partition = validate_partition(candidate, native_dir / "orthohmm_root_hogs.tsv", genes, summary)
    return {"schema_version": 1, "cell": cell["label"], "status": "native_group_output_verified",
            "integrity": integrity, "partition": partition, "native_manifest": file_provenance(native_path),
            "native_metrics": file_provenance(metrics_path), "verifier": file_provenance(Path(__file__)),
            "membership": native["membership_reconciliation"], "species_tree": file_provenance(tree_path),
            "native_outputs_validated": True, "accuracy_evaluated": False,
            "scope": "Root-HOG group benchmark admission only; no pairwise-truth evaluation or independent tree reconstruction."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--index", type=int, choices=range(4), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = validate(args.root.resolve(), args.index)
    with args.output.open("x") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"{report['cell']}: native group output verified; no accuracy evaluated")


if __name__ == "__main__":
    main()
