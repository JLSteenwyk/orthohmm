"""Validate native duplication/loss simulation truth; reject transfer histories."""

import argparse
import csv
import itertools
import json
import math
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

from Bio import Phylo, SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.benchmark_production import file_record


def event_graph(path):
    graph = {}
    previous = -math.inf
    origins = []
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["TIME", "EVENT", "NODES"]:
            raise ValueError("Unexpected event table header")
        for row in reader:
            time = float(row["TIME"])
            if not math.isfinite(time) or time < 0 or time < previous:
                raise ValueError("Invalid event chronology")
            previous = time
            event, parts = row["EVENT"], row["NODES"].split(";")
            if event == "O":
                origins.append(parts)
                continue
            if event not in {"S", "D", "F", "L"}:
                raise ValueError(f"Unsupported event: {event}")
            if len(parts) != (6 if event in {"S", "D"} else 2) or any(not p for p in parts):
                raise ValueError("Malformed event nodes")
            names = [parts[i] + "_" + parts[i + 1] for i in range(0, len(parts), 2)]
            if names[0] in graph:
                raise ValueError("Duplicate terminal event for a gene lineage")
            graph[names[0]] = (event, tuple(names[1:]))
    if origins != [["Root"]]:
        raise ValueError("Only one root-origin family is supported")
    return graph


def xml_graph(path):
    document = ET.parse(path).getroot()
    roots = document.findall("./phylogeny/clade")
    if document.tag != "recGeneTree" or len(roots) != 1:
        raise ValueError("Expected one reconciled gene tree")
    graph = {}
    translate = {"speciation": "S", "duplication": "D", "P": "F", "loss": "L"}
    for clade in roots[0].iter("clade"):
        name = clade.findtext("name")
        events = clade.findall("./eventsRec/*")
        if not name or name in graph or len(events) != 1 or events[0].tag not in translate:
            raise ValueError("Unsupported or ambiguous reconciled event")
        if events[0].get("speciesLocation") != name.rsplit("_", 1)[0]:
            raise ValueError("Reconciled species location mismatch")
        graph[name] = (translate[events[0].tag], tuple(c.findtext("name") for c in clade.findall("clade")))
    return graph


def ortholog_truth(graph, root="Root_1"):
    """Cross-species pairs whose event-labeled common ancestor is speciation."""
    visited, pairs = set(), set()
    def visit(node):
        if node not in graph or node in visited:
            raise ValueError("Missing, cyclic, or multiply parented lineage")
        visited.add(node)
        event, children = graph[node]
        if event not in {"S", "D", "F", "L"} or len(children) != (2 if event in {"S", "D"} else 0):
            raise ValueError("Invalid lineage event/arity")
        if event == "F":
            return {node}
        if event == "L":
            return set()
        left, right = map(visit, children)
        if event == "S":
            for a, b in itertools.product(left, right):
                if a.rsplit("_", 1)[0] == b.rsplit("_", 1)[0]:
                    raise ValueError("Speciation descendants overlap in extant species")
                pairs.add(tuple(sorted((a, b))))
        return left | right
    leaves = visit(root)
    if visited != set(graph):
        raise ValueError("Disconnected lineage events")
    return leaves, pairs


def read_fasta(path):
    records = {}
    for record in SeqIO.parse(path, "fasta"):
        sequence = str(record.seq)
        if record.id in records or not sequence or set(sequence) - set("ACDEFGHIKLMNPQRSTVWY"):
            raise ValueError("Duplicate ID or invalid protein sequence")
        records[record.id] = sequence
    return records


def validate_run(run):
    extant_tree = run / "T/ExtantTree.nwk"
    leaves = [n.name for n in Phylo.read(extant_tree, "newick").get_terminals()]
    if not leaves or len(set(leaves)) != len(leaves):
        raise ValueError("Nonempty unique extant species required")
    genomes = {}
    sources = [extant_tree]
    for species in leaves:
        path = run / "G/Genomes" / f"{species}_GENOME.tsv"
        sources.append(path)
        with path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != ["POSITION", "GENE_FAMILY", "ORIENTATION", "GENE_ID"]:
                raise ValueError("Unexpected genome header")
            for row in reader:
                key = (row["GENE_FAMILY"], species + "_" + row["GENE_ID"])
                if key in genomes:
                    raise ValueError("Duplicate extant genome gene")
                genomes[key] = species
    sequences, families, pairs, checked = {}, {}, [], set()
    event_files = sorted((run / "G/Gene_families").glob("*_events.tsv"))
    if not event_files:
        raise ValueError("No family histories")
    for events in event_files:
        family = events.name.removesuffix("_events.tsv")
        xml = run / "G/Gene_trees" / f"{family}_rec.xml"
        fasta = run / "S" / f"{family}_complete.fasta"
        sources.extend([events, xml, fasta])
        graph = event_graph(events)
        reconciliation = xml_graph(xml)
        # Child ordering is not a biological property.
        normalized = lambda g: {n: (e, tuple(sorted(c))) for n, (e, c) in g.items()}
        if normalized(graph) != normalized(reconciliation):
            raise ValueError("Event history and reconciled tree disagree")
        terminal, orthologs = ortholog_truth(graph)
        expected = {gene for f, gene in genomes if f == family}
        if terminal != expected:
            raise ValueError("Terminal history genes disagree with extant genomes")
        proteins = read_fasta(fasta)
        if not terminal <= proteins.keys():
            raise ValueError("Missing extant sequences")
        native = run / "G/Gene_trees" / f"{family}_prunedtree.nwk"
        sources.append(native)
        text = native.read_text().strip()
        tips = [] if text == ";" else [n.name for n in Phylo.read(native, "newick").get_terminals()]
        if len(tips) != len(set(tips)) or set(tips) != terminal:
            raise ValueError("Pruned gene tree disagrees with terminal events")
        names = {gene: f"F{family}__{gene}" for gene in terminal}
        families[family] = sorted(names.values())
        for gene, name in names.items():
            if name in sequences:
                raise ValueError("Global simulation ID collision")
            sequences[name] = (genomes[(family, gene)], proteins[gene])
            checked.add((family, gene))
        pairs.extend(tuple(sorted((names[a], names[b]))) for a, b in orthologs)
    if checked != set(genomes):
        raise ValueError("Unaccounted extant genome genes")
    return {"schema_version": 1, "scope": "root-origin duplication/loss histories; no transfer",
            "species": sorted(leaves), "extant_genes": len(sequences), "families": families,
            "ortholog_pairs": sorted(pairs), "ortholog_pair_count": len(pairs),
            "inputs": [file_record(p, run) for p in sources]}, sequences


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("Refusing to overwrite extracted truth")
    report, sequences = validate_run(args.run)
    args.output.mkdir(parents=True)
    fasta_dir = args.output / "input"
    fasta_dir.mkdir()
    for species in report["species"]:
        with (fasta_dir / f"{species}.fasta").open("w") as handle:
            for name, (owner, sequence) in sorted(sequences.items()):
                if owner == species:
                    handle.write(f">{name}\n{sequence}\n")
    report["source"] = file_record(Path(__file__).resolve(), Path(__file__).resolve().parent)
    report["prepared_inputs"] = [file_record(p, args.output) for p in sorted(fasta_dir.glob("*.fasta"))]
    (args.output / "truth.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(f"Validated {report['extant_genes']} extant genes and {report['ortholog_pair_count']} ortholog pairs")


if __name__ == "__main__":
    main()
