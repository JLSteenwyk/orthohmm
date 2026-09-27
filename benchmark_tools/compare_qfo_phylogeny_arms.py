"""Full native-output contrasts after separate scientific/runtime admission."""

from collections import Counter, defaultdict
from itertools import combinations
import json

import dendropy

from benchmark_tools.audit_installed_orthobench import read_root_hogs, compare_partitions
from benchmark_tools.audit_phylogeny_structure import rows
from benchmark_tools.readback_qfo_fresh_phylogeny import pair_rows, family_partitions
from benchmark_tools.run_qfo_order_replay import record


def species_tree_signature(path):
    trees = dendropy.TreeList.get(path=str(path), schema="newick", preserve_underscores=True)
    if len(trees) != 1 or trees[0].is_rooted is not True:
        raise ValueError("Require one explicitly rooted species tree")
    tree = trees[0]
    leaves = [node.taxon.label for node in tree.leaf_node_iter() if node.taxon is not None]
    if len(leaves) != len(set(leaves)) or len(leaves) != len(list(tree.leaf_node_iter())):
        raise ValueError("Duplicate or missing species labels")
    descendants, clades = {}, set()
    for node in tree.postorder_node_iter():
        members = (frozenset([node.taxon.label]) if node.is_leaf()
                   else frozenset().union(*(descendants[child] for child in node.child_node_iter())))
        descendants[node] = members
        if 1 < len(members) < len(leaves):
            clades.add(members)
    return frozenset(leaves), frozenset(clades)


def pair_family_contrast(left, right, left_family, right_family):
    left, right = iter(left), iter(right)
    a, b = next(left, None), next(right, None)
    counts = Counter(shared=0, left_only=0, right_only=0, annotations_changed=0)
    changed = defaultdict(Counter)

    def family(row, mapping):
        x, y = row[0]
        if x not in mapping or y not in mapping or mapping[x] != mapping[y]:
            raise ValueError("Pair crosses or lacks a candidate family")
        return mapping[x]

    while a is not None or b is not None:
        if b is None or (a is not None and a[0] < b[0]):
            changed[("left", family(a, left_family))]["unique_pairs"] += 1
            counts["left_only"] += 1
            a = next(left, None)
        elif a is None or b[0] < a[0]:
            changed[("right", family(b, right_family))]["unique_pairs"] += 1
            counts["right_only"] += 1
            b = next(right, None)
        else:
            lf, rf = family(a, left_family), family(b, right_family)
            counts["shared"] += 1
            if a[1] != b[1]:
                counts["annotations_changed"] += 1
                changed[("left", lf)]["annotation_changes"] += 1
                changed[("right", rf)]["annotation_changes"] += 1
            a, b = next(left, None), next(right, None)
    return dict(counts, pair_sets_equal=not(counts["left_only"] or counts["right_only"]),
        confidence_tables_equal=not(counts["left_only"] or counts["right_only"] or counts["annotations_changed"]),
        affected_families=[dict(side=side, family=name, **values)
                           for (side, name), values in sorted(changed.items())])


def compare(paths, universe):
    """Inputs must already pass the four readers and completed-run admission.

    This checks complete predictions and records artifact identities, but does
    not replace the runtime or scientific admission gates.
    """
    watched, roots, families, gene_families, trees, summaries, manifests = [], {}, {}, {}, {}, {}, {}
    for label, path in paths.items():
        for name in ("orthohmm_root_hogs.tsv", "orthohmm_pairwise_orthologs_confidence.tsv",
                     "species_tree.rooted.nwk", "reconciliation_summary.json", "provenance_manifest.json"):
            watched.append(record(path / name))
        root_path = path / "orthohmm_root_hogs.tsv"
        roots[label] = read_root_hogs(root_path, universe)
        families[label] = family_partitions(root_path)
        gene_families[label] = {gene: row["source_family"] for row in
            rows(root_path, ["root_hog", "source_family", "genes"]) for gene in row["genes"].split(",")}
        trees[label] = species_tree_signature(path / "species_tree.rooted.nwk")
        summaries[label] = json.loads((path / "reconciliation_summary.json").read_text())
        manifests[label] = json.loads((path / "provenance_manifest.json").read_text())
    contrasts = {}
    for left, right in combinations(paths, 2):
        changed = [name for name in sorted(families[left].keys() | families[right].keys())
                   if families[left].get(name) != families[right].get(name)]
        pairs = pair_family_contrast(
            pair_rows(paths[left] / "orthohmm_pairwise_orthologs_confidence.tsv"),
            pair_rows(paths[right] / "orthohmm_pairwise_orthologs_confidence.tsv"),
            gene_families[left], gene_families[right])
        lt, rt = trees[left], trees[right]
        contrasts[left + "_vs_" + right] = dict(
            partition=compare_partitions(roots[left], roots[right]), changed_root_families=changed,
            pairs=pairs, species_tree=dict(taxa_equal=lt[0] == rt[0], rooted_topology_equal=lt == rt,
                left_only_clades=len(lt[1] - rt[1]), right_only_clades=len(rt[1] - lt[1]),
                bytes_equal=(paths[left] / "species_tree.rooted.nwk").read_bytes() ==
                            (paths[right] / "species_tree.rooted.nwk").read_bytes()))
    for item in watched:
        if record(item["path"]) != item:
            raise ValueError("Prediction changed during comparison")
    return dict(contrasts=contrasts, summaries=summaries,
        membership={label: m["membership_reconciliation"] for label, m in manifests.items()},
        artifacts=watched, source=record(__file__), accuracy_evaluated=False,
        limitation="Topology equality ignores branch lengths/support; caller must admit native runs and scientific reports")
