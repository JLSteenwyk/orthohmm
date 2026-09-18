"""Validate OrthoMCL 1.4 final groups against the native MCL partition."""

import argparse
from array import array
from collections import Counter
import json
from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.orthomcl_matrix_to_pairwise import load_species
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


INTEGER = re.compile(r"(?:0|[1-9][0-9]*)\Z")
HEADER = re.compile(r"ORTHOMCL(0|[1-9][0-9]*)\(([1-9][0-9]*) genes,([1-9][0-9]*) taxa\):\s*(.*)\Z")
MEMBER = re.compile(r"([^\s()]+)\(([^\s()]+)\)\Z")


def raw_index(path, owners):
    genes, seen = [], set()
    with path.open() as stream:
        for line in stream:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 2 or fields[0] != str(len(genes)):
                raise ValueError("Invalid/noncontiguous native index")
            gene = fields[1]
            if gene not in owners or gene in seen:
                raise ValueError("Unknown/duplicate native index gene")
            genes.append(gene)
            seen.add(gene)
    if not genes:
        raise ValueError("Empty native index")
    return genes


def partition(path, gene_count):
    """Read the pinned MCL matrix format, including wrapped cluster rows."""
    assigned = array("q", [-1]) * gene_count
    sizes, pending = [], []
    state, clusters = 0, None
    with path.open() as stream:
        for line in stream:
            text = line.strip()
            if not text:
                continue
            if state < 6:
                if state == 2:
                    match = re.fullmatch(r"dimensions ([0-9]+)x([0-9]+)", text)
                    if not match or int(match[1]) != gene_count or not 0 < int(match[2]) <= gene_count:
                        raise ValueError("Invalid native partition dimensions")
                    clusters = int(match[2])
                elif text != {0: "(mclheader", 1: "mcltype matrix", 3: ")", 4: "(mclmatrix", 5: "begin"}[state]:
                    raise ValueError("Invalid native partition header")
                state += 1
                continue
            if state == 7:
                raise ValueError("Trailing native partition content")
            if text == ")":
                if pending or len(sizes) != clusters:
                    raise ValueError("Truncated native partition")
                state = 7
                continue
            tokens = text.split()
            ended = tokens[-1] == "$"
            if ended:
                tokens.pop()
            if any(not INTEGER.fullmatch(token) for token in tokens):
                raise ValueError("Invalid native partition token")
            pending.extend(map(int, tokens))
            if ended:
                if len(pending) < 2 or pending[0] != len(sizes) or len(sizes) >= clusters:
                    raise ValueError("Invalid native cluster ID/size")
                cid = pending[0]
                for index in pending[1:]:
                    if not 0 <= index < gene_count or assigned[index] != -1:
                        raise ValueError("Out-of-range/duplicate partition member")
                    assigned[index] = cid
                sizes.append(len(pending) - 1)
                pending.clear()
    if state != 7 or any(cid < 0 for cid in assigned):
        raise ValueError("Incomplete native partition")
    return assigned, sizes


def validate(groups, mcl, index, gg):
    owners = load_species(gg)
    genes = raw_index(index, owners)
    assigned, sizes = partition(mcl, len(genes))
    positions = {gene: position for position, gene in enumerate(genes)}
    seen_genes, seen_clusters = bytearray(len(genes)), set()
    grouped = cross_pairs = same_species_groups = 0
    with groups.open() as stream:
        for line in stream:
            match = HEADER.fullmatch(line.rstrip("\r\n"))
            if not match:
                raise ValueError("Malformed native final-group header")
            cid, count, taxa = map(int, match.group(1, 2, 3))
            if cid >= len(sizes) or cid in seen_clusters or sizes[cid] < 2:
                raise ValueError("Unknown/duplicate/singleton final-group ID")
            tokens = match[4].split()
            if count != sizes[cid] or len(tokens) != count:
                raise ValueError("Final-group count differs from native partition")
            species = Counter()
            for token in tokens:
                member = MEMBER.fullmatch(token)
                if not member or member[1] not in positions:
                    raise ValueError("Malformed/unknown final-group member")
                gene, taxon = member.groups()
                position = positions[gene]
                if seen_genes[position] or assigned[position] != cid:
                    raise ValueError("Duplicate/wrong-cluster final-group member")
                if owners[gene] != taxon:
                    raise ValueError("Final-group taxon differs from GG")
                seen_genes[position] = 1
                species[taxon] += 1
            if taxa != len(species):
                raise ValueError("Final-group taxon count differs from GG")
            grouped += count
            same_species_groups += len(species) == 1
            cross_pairs += (count * count - sum(n * n for n in species.values())) // 2
            seen_clusters.add(cid)
    if len(seen_clusters) != sum(size >= 2 for size in sizes):
        raise ValueError("Missing nonsingleton final groups")
    return {
        "input_proteins": len(owners), "input_species": len(set(owners.values())),
        "indexed_proteins": len(genes), "mcl_clusters": len(sizes),
        "mcl_singleton_clusters": sizes.count(1), "final_groups": len(seen_clusters),
        "grouped_proteins": grouped, "ungrouped_input_proteins": len(owners) - grouped,
        "input_proteins_absent_from_index": len(owners) - len(genes),
        "single_species_final_groups": same_species_groups,
        "cross_species_clique_pairs": cross_pairs,
    }


def audit(groups, mcl, index, gg, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    checked = [record(path) for path in (groups, mcl, index, gg, Path(__file__),
        Path(__file__).with_name("orthomcl_matrix_to_pairwise.py"),
        Path(__file__).with_name("prepare_ob_candidate_neighborhood.py"))]
    content = validate(groups, mcl, index, gg)
    for item in checked:
        check(item)
    result = {"status": "native_final_groups_match_mcl_partition", "content": content,
              "checked_records": checked, "accuracy_admitted": False, "publication_ready": False,
              "limitations": [
                  "Structural final-output validation, not execution provenance or biological accuracy.",
                  "Singleton clusters are omitted by native OrthoMCL 1.4 mcl_backindex; no singletons are added.",
                  "Pair count describes cross-species final-group cliques, not native pairwise orthologs or graph edges.",
                  "Source search failures and full input coverage require separate audits."]}
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("groups", "mcl", "index", "gg", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.groups.resolve(), args.mcl.resolve(), args.index.resolve(),
                   args.gg.resolve(), args.output.absolute())
    print(json.dumps({"status": result["status"], "content": result["content"]}))
