"""Separate saved root-lineage and satellite-constraint membership losses."""

import argparse
from collections import Counter, defaultdict
import csv
from itertools import combinations
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
from orthohmm.phylogeny import canonical_species_name
from benchmark_tools.assemble_orthobench_factorial import load_reference_snapshot
from benchmark_tools.audit_historical_profile_ablation import verify_file
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.reconstruct_reconciliation_trace import reconstruct_nodes, reconstruct_bypass, apply_logged_constraints
from benchmark_tools.recover_factorial_postflight import PREPARED_HASH
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.score_ygob_groups import membership
from benchmark_tools.trace_ob_families import FACTORIAL_SHA, partition

TRACE_SHA = "bda00fb593b357bc8f07e43544feae598150ae891ed82c988fb339a92349de0c"


def assemble(root, output):
    if output.exists():
        raise FileExistsError(output)
    results = root / "benchmark_tools/results"
    trace_path = results / "ob_family_trace_verified_20260916.json"
    trace = read_frozen(trace_path, TRACE_SHA)
    prepared = read_frozen(results / "orthobench_factorial_prepared_20260916.json", PREPARED_HASH)
    factorial = read_frozen(results / "orthobench_factorial_results_20260916.json", FACTORIAL_SHA)
    native = factorial["native_validation"]["p1_c1_r1"]
    checked = []
    def check(item, path=None):
        path = Path(path or item.get("absolute_path") or item["path"])
        verify_file(path, item)
        checked.append(file_provenance(path))
        return path
    for name, digest in (("phylogeny.py", "216d97608e6dede8960f3b06fe7036bd23da8fae5677df88583ba9369ec3c1bf"),
                         ("phylogeny_pipeline.py", "44e00316546b5df78354badd2b1a6bb595b685e98b86f686f1a32543d7f15e4f")):
        item = file_provenance(root / "benchmarks/work/publication_method_native_v2/orthohmm" / name)
        if item["sha256"] != digest:
            raise ValueError("Frozen reconciliation source differs from reviewed rules")
        check(item)
    for name in ("ob_family_trace_verified_20260916.json", "orthobench_factorial_prepared_20260916.json", "orthobench_factorial_results_20260916.json"):
        check(file_provenance(results / name))
    status_path = check(native["integrity"]["original_status"])
    status = json.loads(status_path.read_text())["methods"]["p1_c1_r1"]
    if status["status"] != "process_succeeded" or status["exit_code"] != 0:
        raise ValueError("Native inference did not succeed")
    artifacts = {item["absolute_path"]: item for item in status["outputs"]}
    manifest = json.loads(check(native["native_manifest"]).read_text())
    if (manifest["root_duplication_rule"] != "species_overlap" or manifest["pair_orthology_rule"] != "positive_paralogy"
            or manifest["membership_reconciliation"]["policy"] != "high_confidence_pair"):
        raise ValueError("Unsupported reconciliation or membership rules")
    check(native["species_tree"])
    owners = {}
    for item in prepared["fasta_inputs"]:
        path = check(item)
        for protein in SeqIO.parse(path, "fasta"):
            if protein.id in owners:
                raise ValueError("Duplicate gene ID")
            owners[protein.id] = canonical_species_name(path.name)
    references, _, _, ref_records = load_reference_snapshot(results / "orthobench_paired_uncertainty_20260916.json")
    checked.extend(ref_records)
    reference_genes = set().union(*references.values())
    final_path = check(factorial["predictions"]["p1_c1_r1"])
    directory = final_path.parent
    final_groups, candidate_genes, final_family = {}, defaultdict(set), {}
    with final_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["root_hog", "source_family", "genes"]:
            raise ValueError("Unexpected root-HOG schema")
        for row in reader:
            key, family, genes = row["root_hog"], row["source_family"], row["genes"].split(",")
            if key in final_groups or not family:
                raise ValueError("Invalid root-HOG identity")
            final_groups[key] = genes
            final_family[key] = family
            candidate_genes[family].update(genes)
    final_index = membership(final_groups)
    if set(final_index) != set(owners):
        raise ValueError("Incomplete final gene universe")
    candidates, _ = partition(check(factorial["predictions"]["p1_c1_r0"]), "plain", set(owners))
    if {frozenset(g) for g in candidates.values()} != {frozenset(g) for g in candidate_genes.values()}:
        raise ValueError("Source-family labels do not recover candidate partition")
    candidate_index = {gene: family for family, genes in candidate_genes.items() for gene in genes}
    selected = {candidate_index[g] for g in reference_genes}
    rows = defaultdict(list)
    nodes_path = directory / "orthohmm_reconciliation_nodes.tsv"
    check(artifacts[str(nodes_path)])
    with nodes_path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["source_family"] in selected:
                rows[row["source_family"]].append(row)
    constraints_path = check(prepared["candidate_arms"]["p1_c1"]["membership_constraints"])
    constraints = defaultdict(list)
    for i, event in enumerate(json.loads(constraints_path.read_text())):
        families = {candidate_index[g] for key in ("source_genes", "target_genes") for g in event[key]}
        if len(families) != 1:
            raise ValueError("Constraint crosses candidate boundaries")
        family = next(iter(families))
        if family in selected:
            constraints[family].append((i, event))
    candidate_results, pre_index = {}, {}
    for family in sorted(selected):
        genes = candidate_genes[family]
        if rows[family]:
            checkpoint_path = directory / "checkpoints" / (family + ".json")
            checkpoint = json.loads(check(artifacts[str(checkpoint_path)]).read_text())
            if (checkpoint["status"] != "complete" or checkpoint["family_id"] != family
                    or set(checkpoint["genes"]) != genes or checkpoint["species_tree_sha256"] != native["species_tree"]["sha256"]):
                raise ValueError("Tree checkpoint differs from candidate or species tree")
            for suffix, key in (("raw", "raw_tree_sha256"), ("rooted", "rooted_tree_sha256"), ("reconciled", "annotated_tree_sha256")):
                path = directory / "gene_trees" / (family + "." + suffix + ".nwk")
                record = artifacts[str(path)]
                check(record)
                if record["sha256"] != checkpoint[key]:
                    raise ValueError("Tree/checkpoint hash disagreement")
            reconstruction = reconstruct_nodes(rows[family], genes, owners)
        else:
            reconstruction = reconstruct_bypass(genes, owners)
        groups, constraint_evidence = apply_logged_constraints(reconstruction, constraints[family], genes)
        expected = {frozenset(group) for key, group in final_groups.items() if final_family[key] == family}
        if {frozenset(group) for group in groups} != expected:
            raise ValueError("Reconstructed root groups and constraints differ from native final membership")
        for i, group in enumerate(reconstruction["root_groups"]):
            pre_index.update({g: (family, i) for g in group})
        candidate_results[family] = {"genes": len(genes), "reference_genes": sorted(genes & reference_genes),
            "reconciled": bool(rows[family]), "root_duplication_nodes": reconstruction["root_duplication_nodes"],
            "propagated_split_nodes": reconstruction["propagated_split_nodes"],
            "preconstraint_root_groups": [sorted(group) for group in reconstruction["root_groups"]],
            "final_root_groups": [sorted(group) for group in groups], "constraints": constraint_evidence}
    families = {}
    for name, genes in sorted(references.items()):
        counts, lost = Counter(), []
        for a, b in combinations(sorted(genes), 2):
            category = ("different_candidates" if candidate_index[a] != candidate_index[b] else
                        "root_lineage_split" if pre_index[a] != pre_index[b] else
                        "satellite_constraint_split" if final_index[a] != final_index[b] else "retained")
            counts[category] += 1
            if category in {"root_lineage_split", "satellite_constraint_split"}:
                lost.append({"gene_a": a, "gene_b": b, "candidate_family": candidate_index[a], "category": category})
        prior = trace["families"][name]
        if (counts["retained"] != prior["stages"]["root_hogs"]["within_family_pairs"]
                or len(lost) != prior["transitions"]["candidates_to_root_hogs"]["lost"]):
            raise ValueError("Loss decomposition differs from verified pair trace")
        families[name] = {"counts": {key: counts[key] for key in
            ("different_candidates", "root_lineage_split", "satellite_constraint_split", "retained")}, "lost_pairs": lost}
    for item in checked:
        verify_file(Path(item["path"]), item)
    report = {"status": "reference_incident_native_membership_reconstructed", "publication_ready": False,
        "source": file_provenance(Path(__file__)), "reconstruction_source": file_provenance(Path(__file__).with_name("reconstruct_reconciliation_trace.py")),
        "prior_trace": file_provenance(trace_path), "checked_inputs": checked,
        "candidate_families": candidate_results, "reference_families": families,
        "limitations": ["Native node calls are reproduced, not independently validated evolutionary history or rerooted gene trees.",
            "S0000 is the species root under the frozen preorder-label convention; only species_overlap/positive_paralogy rules are supported.",
            "The decomposition follows the executed stage order; pairs split by root-lineage rules are not also counted as constraint losses.",
            "Only reference-incident candidates are reconstructed; this does not audit every genome-wide family or global membership-audit count.",
            "Reference pairs include within-species and low-certainty assignments, and memberships may overlap; counts are not official accuracy statistics.",
            "The recorded node table supplies rooting and mapping decisions; tree error, missing edge evidence and independent biological annotation remain unresolved."]}
    output.mkdir(parents=True)
    (output / "results.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    assemble(args.root.resolve(), args.output.resolve())


if __name__ == "__main__":
    main()
