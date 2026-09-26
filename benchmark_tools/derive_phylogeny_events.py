"""Independent oracle for frozen species-overlap/positive-paralogy rules."""

from collections import Counter
from itertools import combinations, product
import math


NODE_COLUMNS = ["source_family", "node_id", "parent_node_id", "event", "species_tree_node",
                "species", "genes", "pair_event", "event_confidence", "species_overlap_count",
                "mapping_conflict", "branch_support"]


def derive(tree, species_tree, gene_species, family):
    """Return node rows, pre-constraint root groups and native candidate pairs."""
    species_tree.is_rooted = True
    species_labels = {}
    index = 0
    for node in species_tree.preorder_node_iter():
        if node.is_leaf():
            species_labels[node] = node.taxon.label
        else:
            species_labels[node] = f"S{index:04d}"
            index += 1
    taxa = {node.taxon.label for node in species_tree.leaf_node_iter()}
    nodes = list(tree.postorder_node_iter())
    genes, species, mapped, ids = {}, {}, {}, {}
    rows, cuts, pairs = {}, {}, {}
    index = 0
    for node in nodes:
        children = node.child_nodes()
        support = None
        if node.is_leaf():
            label = node.taxon.label
            if label not in gene_species or gene_species[label] not in taxa:
                raise ValueError("Unknown gene/species mapping")
            genes[node] = frozenset([label])
            ids[node] = label
        else:
            genes[node] = frozenset().union(*(genes[child] for child in children))
            if sum(len(genes[child]) for child in children) != len(genes[node]):
                raise ValueError("Repeated gene-tree leaf")
            ids[node] = f"G{index:05d}"
            index += 1
            try:
                support = float(node.label)
                if 1 < support <= 100:
                    support /= 100
                if not math.isfinite(support) or not 0 <= support <= 1:
                    support = None
            except (TypeError, ValueError):
                support = None
        species[node] = {gene_species[gene] for gene in genes[node]}
        mapped[node] = species_tree.mrca(taxon_labels=species[node])
        multiplicities = Counter(taxon for child in children for taxon in species[child])
        overlap = sum(count > 1 for count in multiplicities.values())
        conflict = any(mapped[child] is mapped[node] for child in children)
        if node.is_leaf():
            event, pair_event, confidence = "leaf", "leaf", "not_applicable"
        else:
            event = "duplication" if overlap or conflict else "speciation"
            pair_event = "duplication" if overlap else ("uncertain" if conflict else "speciation")
            confidence = "high" if (overlap >= 2 or (overlap and (support or 0) >= .9)
                                     or (not overlap and not conflict)) else "medium"
        cuts[node] = bool(overlap and mapped[node] is species_tree.seed_node)
        rows[node] = dict(source_family=family, node_id=ids[node], event=event,
            species_tree_node=species_labels[mapped[node]], species=",".join(sorted(species[node])),
            genes=",".join(sorted(genes[node])), pair_event=pair_event, event_confidence=confidence,
            species_overlap_count=str(overlap), mapping_conflict=str(conflict).lower(),
            branch_support="" if support is None else f"{support:.6f}")
        if not node.is_leaf() and pair_event != "duplication":
            for left, right in combinations(children, 2):
                for a, b in product(genes[left], genes[right]):
                    if gene_species[a] != gene_species[b]:
                        pairs[tuple(sorted((a, b)))] = confidence
    for node in nodes:
        rows[node]["parent_node_id"] = ids[node.parent_node] if node.parent_node else ""
    # Split every ancestor of a root-mapped duplication, then keep maximal clades.
    split = set()
    for node in nodes:
        if cuts[node]:
            while node is not None:
                split.add(node)
                node = node.parent_node
    groups = [genes[node] for node in nodes if node not in split
              and (node.parent_node is None or node.parent_node in split)]
    return [rows[node] for node in nodes], groups, pairs


def constrain(groups, pairs, constraints, enabled):
    """Apply source detachment and root-group pair restriction independently."""
    detached, supported = {}, 0
    universe = set().union(*groups)
    for index, constraint in enumerate(constraints):
        source, target = set(constraint["source_genes"]), set(constraint["target_genes"])
        if not source or not target or source & target or not source | target <= universe:
            raise ValueError("Invalid satellite constraint membership")
        if any(pairs.get(tuple(sorted((a, b)))) == "high" for a, b in product(source, target)):
            supported += 1
        else:
            if source & detached.keys():
                raise ValueError("Overlapping detached sources")
            detached.update(dict.fromkeys(source, index))
    refined = []
    for group in groups:
        keys = {detached.get(gene, -1) for gene in group}
        refined.extend(frozenset(gene for gene in group if detached.get(gene, -1) == key) for key in keys)
    by_gene = {gene: index for index, group in enumerate(refined) for gene in group}
    retained = {pair: confidence for pair, confidence in pairs.items()
                if not enabled or by_gene[pair[0]] == by_gene[pair[1]]}
    counts = dict(constraints=len(constraints), supported_constraints=supported,
        detached_constraints=len(constraints) - supported, detached_genes=len(detached),
        root_hogs_added=len(refined) - len(groups), ortholog_pairs_removed=len(pairs) - len(retained))
    return refined, retained, counts
