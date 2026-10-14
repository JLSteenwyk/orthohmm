"""Contrast native membership and pair multisets, not benchmark accuracy."""

import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_membership_hillclimb import write_json
from benchmark_tools.score_ygob_groups import membership, read_predictions

PAIR_COLUMNS = ["gene_a", "species_a", "gene_b", "species_b"]
MASK = (1 << 32) - 1


def pair_multiset(path, index, owners):
    pairs = Counter()
    with path.open() as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames != PAIR_COLUMNS:
            raise ValueError("Unexpected native pair header")
        for row in reader:
            if None in row or any(not value or value.strip() != value for value in row.values()):
                raise ValueError("Malformed native pair row")
            genes = row["gene_a"], row["gene_b"]
            species = row["species_a"], row["species_b"]
            if genes[0] == genes[1] or species[0] == species[1] or any(gene not in index for gene in genes):
                raise ValueError("Self-pair, same-species pair, or gene outside the root universe")
            for gene, owner in zip(genes, species):
                if owners.setdefault(gene, owner) != owner:
                    raise ValueError("Gene species ownership differs within or between native files")
            a, b = sorted(index[gene] for gene in genes)
            if b > MASK:
                raise ValueError("Gene index exceeds lossless pair encoding")
            # Packed IDs keep multi-million-pair contrasts bounded in memory.
            pairs[(a << 32) | b] += 1
    return pairs


def contrast(left, right):
    removed, added = left - right, right - left
    common = sum(min(count, right.get(key, 0)) for key, count in left.items())
    return {"baseline_records": sum(left.values()), "candidate_records": sum(right.values()),
            "common_records": common, "removed_records": sum(removed.values()),
            "added_records": sum(added.values()), "semantic_multisets_equal": not removed and not added}, removed, added


def compare(baseline, candidate, expected_genes=None):
    folders = baseline, candidate
    names = ("orthohmm_root_hogs.tsv", "orthohmm_pairwise_orthologs.tsv", "reconciliation_summary.json")
    pins = [record(folder / name) for folder in folders for name in names]
    pins.extend([record(Path(__file__)), record(Path(__file__).with_name("score_ygob_groups.py"))])
    groups = [read_predictions(folder / names[0], "root_hogs") for folder in folders]
    universes = [set(membership(value)) for value in groups]
    if universes[0] != universes[1] or (expected_genes is not None and len(universes[0]) != expected_genes):
        raise ValueError("Native root gene universes differ or have the wrong size")
    genes = sorted(universes[0])
    index = {gene: i for i, gene in enumerate(genes)}
    roots = [Counter(tuple(sorted(index[gene] for gene in group)) for group in value.values()) for value in groups]
    summaries = [json.loads((folder / names[2]).read_text()) for folder in folders]
    owners = {}
    pairs = [pair_multiset(folder / names[1], index, owners) for folder in folders]
    for summary, root, pair in zip(summaries, roots, pairs):
        if (type(summary.get("root_hogs")) is not int or type(summary.get("ortholog_pairs")) is not int
                or summary["root_hogs"] != sum(root.values()) or summary["ortholog_pairs"] != sum(pair.values())):
            raise ValueError("Native reconciliation counts differ from parsed records")
    root_result, removed_roots, added_roots = contrast(*roots)
    pair_result, removed_pairs, added_pairs = contrast(*pairs)
    changed_pairs = removed_pairs.keys() | added_pairs.keys()
    affected_pair_genes = {i for key in changed_pairs for i in (key >> 32, key & MASK)}
    affected_root_genes = {i for group in removed_roots.keys() | added_roots.keys() for i in group}
    pair_result.update(affected_genes=len(affected_pair_genes),
        baseline_unique_pairs=len(pairs[0]), candidate_unique_pairs=len(pairs[1]),
        baseline_duplicate_rows=sum(pairs[0].values()) - len(pairs[0]),
        candidate_duplicate_rows=sum(pairs[1].values()) - len(pairs[1]))
    root_result["genes_in_changed_groups"] = len(affected_root_genes)
    samples = {}
    for label, values in (("removed", removed_pairs), ("added", added_pairs)):
        samples[label] = [{"pair": [[genes[i], owners[genes[i]]] for i in (key >> 32, key & MASK)],
                           "count": values[key]} for key in sorted(values)[:10]]
    tree_paths = [folder / "species_tree.rooted.nwk" for folder in folders]
    tree_match = None
    if all(path.is_file() for path in tree_paths):
        trees = [record(path) for path in tree_paths]
        pins.extend(trees)
        tree_match = (trees[0]["bytes"], trees[0]["sha256"]) == (trees[1]["bytes"], trees[1]["sha256"])
    for pin in pins:
        check(pin)
    return {"status": "native_membership_and_pair_multiset_contrast_completed", "genes": len(genes),
            "baseline": str(baseline), "candidate": str(candidate), "root_membership": root_result,
            "pair_membership": pair_result, "pair_samples": samples, "species_tree_same_bytes": tree_match,
            "checked": pins, "accuracy_evaluated": False, "defaults_promoted": False, "goal_completed": False,
            "scope": "Ignores row order and arbitrary root/family IDs; preserves gene identity, pair species ownership and multiplicity. Confidence values and tree-event annotations are not compared. No claim of biological correctness or a uniquely identified causal stage."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("baseline", "candidate", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--expected-genes", type=int, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.expected_genes <= 0:
        raise ValueError("Require fresh diagnostic output and a positive expected universe")
    result = compare(args.baseline.resolve(), args.candidate.resolve(), args.expected_genes)
    write_json(args.output, result)
    print(json.dumps({k: result[k] for k in ("status", "genes", "root_membership", "pair_membership")}, indent=2))
