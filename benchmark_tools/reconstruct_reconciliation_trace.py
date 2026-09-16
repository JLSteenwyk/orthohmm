"""Reconstruct fixed species-overlap root groups from saved reconciliation nodes."""

from collections import Counter, defaultdict
from itertools import combinations, product


def reconstruct_nodes(rows, genes, owners):
    """Return pre-constraint groups and high-confidence pairs, without tree inference."""
    nodes = {}
    children = defaultdict(list)
    for row in rows:
        key = row["node_id"]
        if not key or key in nodes:
            raise ValueError("Empty or duplicate reconciliation node ID")
        members = row["genes"].split(",")
        if not members or len(set(members)) != len(members) or not set(members) <= genes:
            raise ValueError("Invalid node gene inventory")
        nodes[key] = {**row, "members": set(members)}
        children[row["parent_node_id"]].append(key)
    if len(children[""]) != 1 or any(key and key not in nodes for key in children):
        raise ValueError("Node table does not have one valid root")
    root = children[""][0]
    active, visited = set(), set()
    root_duplications, propagated_splits, high_pairs = [], [], set()
    def visit(key):
        if key in active or key in visited:
            raise ValueError("Cyclic or repeated node traversal")
        active.add(key)
        row = nodes[key]
        offspring = children.get(key, [])
        species = {owners[g] for g in row["members"]}
        if set(row["species"].split(",")) != species:
            raise ValueError("Recorded node species differ from gene ownership")
        if not offspring:
            if row["event"] != "leaf" or len(row["members"]) != 1 or row["pair_event"] != "leaf":
                raise ValueError("Invalid leaf node")
            contains_duplication, groups = False, [row["members"]]
        else:
            descendants = [visit(child) for child in offspring]
            child_genes = [nodes[child]["members"] for child in offspring]
            if (len(offspring) < 2 or set().union(*child_genes) != row["members"]
                    or sum(map(len, child_genes)) != len(row["members"])):
                raise ValueError("Children do not partition parent genes")
            child_species = [{owners[g] for g in group} for group in child_genes]
            overlaps = set().union(*(a & b for a, b in combinations(child_species, 2)))
            if int(row["species_overlap_count"]) != len(overlaps) or row["mapping_conflict"] not in {"true", "false"}:
                raise ValueError("Recorded overlap or mapping status invalid")
            conflict = row["mapping_conflict"] == "true"
            expected_pair = "duplication" if overlaps else "uncertain" if conflict else "speciation"
            expected_event = "duplication" if overlaps or conflict else "speciation"
            if row["pair_event"] != expected_pair or row["event"] != expected_event:
                raise ValueError("Node calls differ from fixed positive-paralogy rules")
            if expected_pair == "speciation" and row["event_confidence"] != "high":
                raise ValueError("Speciation confidence differs from fixed rules")
            if expected_pair == "uncertain" and row["event_confidence"] != "medium":
                raise ValueError("Uncertain confidence differs from fixed rules")
            if row["pair_event"] == "speciation":
                for a, b in combinations(child_genes, 2):
                    high_pairs.update(tuple(sorted((left, right))) for left, right in product(a, b) if owners[left] != owners[right])
            is_root_duplication = bool(overlaps) and row["species_tree_node"] == "S0000"
            contains_duplication = is_root_duplication or any(value[0] for value in descendants)
            if is_root_duplication:
                root_duplications.append(key)
            elif contains_duplication:
                propagated_splits.append(key)
            groups = [group for _, parts in descendants for group in parts] if contains_duplication else [row["members"]]
        active.remove(key)
        visited.add(key)
        return contains_duplication, groups
    _, groups = visit(root)
    if visited != set(nodes) or nodes[root]["members"] != genes:
        raise ValueError("Disconnected nodes or incomplete candidate genes")
    return {"root_groups": groups, "high_confidence_pairs": high_pairs,
            "root_duplication_nodes": sorted(root_duplications), "propagated_split_nodes": sorted(propagated_splits)}


def reconstruct_bypass(genes, owners):
    counts = Counter(owners[g] for g in genes)
    if len(genes) >= 3 and len(counts) >= 2 and max(counts.values()) > 1:
        raise ValueError("Ambiguous candidate is missing reconciliation nodes")
    return {"root_groups": [genes], "high_confidence_pairs": {tuple(sorted(pair)) for pair in combinations(genes, 2) if owners[pair[0]] != owners[pair[1]]},
            "root_duplication_nodes": [], "propagated_split_nodes": []}


def apply_logged_constraints(reconstruction, constraints, genes):
    detached, evidence = {}, []
    for event_index, row in constraints:
        source, target = (set(row[key]) for key in ("source_genes", "target_genes"))
        if not source or not target or source & target or not (source | target) <= genes:
            raise ValueError("Invalid constraint sides or candidate boundary")
        supporting_pair = next((pair for pair in sorted(reconstruction["high_confidence_pairs"])
                                if (pair[0] in source and pair[1] in target) or (pair[1] in source and pair[0] in target)), None)
        evidence.append({"event_index": event_index, "supported": supporting_pair is not None,
                         "supporting_pair": supporting_pair, "source_genes": sorted(source), "target_genes": sorted(target)})
        if supporting_pair is None:
            for gene in source:
                if gene in detached:
                    raise ValueError("Gene occurs in multiple detached sources")
                detached[gene] = event_index
    groups = []
    for group in reconstruction["root_groups"]:
        parts = defaultdict(set)
        for gene in group:
            parts[detached.get(gene, -1)].add(gene)
        groups.extend(parts.values())
    return sorted(groups, key=lambda group: tuple(sorted(group))), evidence
