"""Reconstruct a deterministic sample of native pairs from saved node events."""

import argparse
from collections import defaultdict
import csv
import hashlib
from itertools import combinations, product
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


def reconstruct(rows, groups, constrained):
    nodes = {r["node_id"]: r for r in rows}
    if not nodes or len(nodes) != len(rows):
        raise ValueError("Empty or duplicate node inventory")
    children = defaultdict(list)
    for row in rows:
        parent = row["parent_node_id"]
        if parent and parent not in nodes:
            raise ValueError("Unknown parent")
        children[parent].append(row["node_id"])
    if len(children[""]) != 1:
        raise ValueError("Require one root")
    pending = [(children[""][0], False)]
    visited, descendants, owners, pairs = set(), {}, {}, set()
    while pending:
        name, processed = pending.pop()
        row = nodes[name]
        if not processed:
            if name in visited:
                raise ValueError("Cycle or reused node")
            visited.add(name)
            pending.append((name, True))
            pending.extend((child, False) for child in children[name])
            continue
        genes = row["genes"].split(",")
        if not all(genes) or len(genes) != len(set(genes)):
            raise ValueError("Empty or repeated descendant gene")
        if not children[name]:
            if row["event"] != "leaf" or row["pair_event"] != "leaf" or len(genes) != 1 or "," in row["species"] or not row["species"]:
                raise ValueError("Invalid leaf")
            if genes[0] in owners:
                raise ValueError("Duplicate leaf gene")
            owners[genes[0]] = row["species"]
        else:
            blocks = [descendants[c] for c in children[name]]
            flat = [g for block in blocks for g in block]
            if len(blocks) < 2 or len(flat) != len(set(flat)) or set(flat) != set(genes):
                raise ValueError("Children do not partition descendants")
            taxa = [{owners[g] for g in block} for block in blocks]
            overlap = set().union(*(a & b for a, b in combinations(taxa, 2)))
            if row["mapping_conflict"] not in {"true", "false"} or int(row["species_overlap_count"]) != len(overlap):
                raise ValueError("Invalid recorded conflict or overlap")
            conflict = row["mapping_conflict"] == "true"
            event = "duplication" if overlap or conflict else "speciation"
            pair_event = "duplication" if overlap else "uncertain" if conflict else "speciation"
            if row["event"] != event or row["pair_event"] != pair_event:
                raise ValueError("Node labels disagree with positive-paralogy rule")
            if pair_event != "duplication":
                for a, b in combinations(blocks, 2):
                    pairs.update(tuple(sorted((x, y))) for x, y in product(a, b) if owners[x] != owners[y])
        if set(row["species"].split(",")) != {owners[g] for g in genes}:
            raise ValueError("Recorded descendant species disagree")
        descendants[name] = genes
    if visited != set(nodes) or set(groups) != set(owners):
        raise ValueError("Disconnected nodes or incomplete root-HOG coverage")
    unfiltered = len(pairs)
    if constrained:
        pairs = {pair for pair in pairs if groups[pair[0]] == groups[pair[1]]}
    return pairs, {"nodes": len(nodes), "genes": len(owners), "event_pairs": unfiltered,
                   "retained_pairs": len(pairs), "group_filtered_pairs": unfiltered - len(pairs)}


def audit(admission_path, expected_sha, count):
    admission_record = record(admission_path)
    if admission_record["sha256"] != expected_sha or not 1 <= count <= 256:
        raise ValueError("Unreviewed admission or invalid sample size")
    admission = json.loads(admission_path.read_text())
    if admission["status"] != "qfo_native_pair_output_verified" or admission["cell"] not in {"p0_c0_r1", "p0_c1_r1", "p1_c0_r1", "p1_c1_r1"}:
        raise ValueError("Unsupported admission")
    directory = Path(admission["native_pairs"]["path"]).parent
    check(admission["native_pairs"])
    manifest_record = admission["native_group_integrity"]["native_manifest"]
    check(manifest_record)
    execution_record = admission["native_group_integrity"]["integrity"]["execution_status"]
    check(execution_record)
    execution = json.loads(Path(execution_record["path"]).read_text())
    artifacts = execution["methods"][admission["cell"]]["outputs"]
    by_path = {r["absolute_path"]: r for r in artifacts}
    if len(by_path) != len(artifacts):
        raise ValueError("Duplicate execution artifact path")
    paths = [directory / name for name in ("orthohmm_reconciliation_nodes.tsv", "orthohmm_root_hogs.tsv")]
    identities = [record(p) for p in paths]
    for identity in identities:
        original = by_path.get(identity["path"])
        if original is None or any(original[k] != identity[k] for k in ("bytes", "sha256")):
            raise ValueError("Node or group table differs from admitted execution")
    # Family IDs are selected without scores, reference labels or observed errors.
    available = {Path(p).name.removesuffix(".reconciled.nwk") for p in by_path
                 if Path(p).parent == directory / "gene_trees" and p.endswith(".reconciled.nwk")}
    selected = sorted(available, key=lambda s: (hashlib.sha256(("20260923:" + s).encode()).hexdigest(), s))[:count]
    if len(selected) != count:
        raise ValueError("Insufficient reconciled families")
    rows, groups = {s: [] for s in selected}, {s: {} for s in selected}
    with paths[0].open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if row["source_family"] in rows:
                rows[row["source_family"]].append(row)
    with paths[1].open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if row["source_family"] in groups:
                group = groups[row["source_family"]]
                for gene in row["genes"].split(","):
                    if gene in group:
                        raise ValueError("Repeated root-HOG gene")
                    group[gene] = row["root_hog"]
    gene_family = {}
    expected, summaries = {}, {}
    for family in selected:
        expected[family], summaries[family] = reconstruct(rows[family], groups[family], "_c1_" in admission["cell"])
        for gene in groups[family]:
            if gene in gene_family:
                raise ValueError("Shared sampled gene")
            gene_family[gene] = family
    observed = {s: set() for s in selected}
    with Path(admission["native_pairs"]["path"]).open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            a, b = row["gene_a"], row["gene_b"]
            family = gene_family.get(a, gene_family.get(b))
            if family is None:
                continue
            if gene_family.get(a) != family or gene_family.get(b) != family:
                raise ValueError("Native pair leaves sampled family")
            pair = tuple(sorted((a, b)))
            if pair in observed[family]:
                raise ValueError("Duplicate sampled native pair")
            observed[family].add(pair)
    for family in selected:
        if observed[family] != expected[family]:
            raise ValueError(f"Reconstructed pair mismatch for {family}: missing={len(expected[family] - observed[family])}, extra={len(observed[family] - expected[family])}")
    for identity in [admission_record, admission["native_pairs"], manifest_record, execution_record, *identities]:
        check(identity)
    return {"status": "sampled_recorded_event_pair_conversion_verified", "cell": admission["cell"],
            "source": record(__file__), "admission": admission_record, "inputs": identities,
            "execution_status": execution_record,
            "native_pairs": admission["native_pairs"], "selection_seed": 20260923,
            "eligible_reconciled_families": len(available), "selected_families": selected,
            "families": summaries, "accuracy_evaluated": False,
            "limitations": ["Sampled reconciled families only; bypass families and unsampled pairs not reconstructed.",
                            "Validates recorded node-event to pair conversion, not independent tree inference or biological orthology truth.",
                            "Mapping-conflict labels and final root-HOG membership are conditioned on, not independently reconciled or regenerated.",
                            "Confidence labels and correctness of membership-constraint selection are not tested."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("admission", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--admission-sha256", required=True)
    parser.add_argument("--families", type=int, default=64)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.admission, args.admission_sha256, args.families)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
