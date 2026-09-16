#!/usr/bin/env python3
"""Run the frozen descriptive YGOB/development homology screen, not accuracy scoring."""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import re
import subprocess
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.orthobench_stage_diagnostics import file_provenance


def search_command(binary, query, database, output, threads):
    return [str(binary), "blastp", "--query", str(query), "--db", str(database),
            "--out", str(output), "--threads", str(threads), "--very-sensitive",
            "--evalue", "1e-5", "--id", "30", "--query-cover", "50", "--subject-cover", "50",
            "--max-target-seqs", "1", "--max-hsps", "1", "--outfmt", "6",
            "qseqid", "sseqid", "pident", "qcovhsp", "scovhsp", "evalue", "bitscore"]


def summarize_hits(path, query_species, references, target_counts, labels):
    gene_to_reference = {}
    for name, genes in references.items():
        for gene in genes:
            if gene in gene_to_reference or gene not in query_species:
                raise ValueError(f"Invalid reference membership: {gene}")
            gene_to_reference[gene] = name
    matched = set()
    best_sources = Counter()
    with path.open() as handle:
        for number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 7:
                raise ValueError(f"Expected seven columns at {path}:{number}")
            query, target = fields[:2]
            if query not in query_species or query in matched:
                raise ValueError(f"Unknown or repeated query: {query}")
            match = re.fullmatch(r"D(\d+)G(\d+)", target)
            if not match:
                raise ValueError(f"Invalid target identifier: {target}")
            dataset, ordinal = map(int, match.groups())
            if dataset >= len(target_counts) or not 1 <= ordinal <= target_counts[dataset]:
                raise ValueError(f"Target outside indexed range: {target}")
            identity, qcover, scover, evalue, bitscore = map(float, fields[2:])
            if not all(math.isfinite(v) for v in (identity, qcover, scover, evalue, bitscore)):
                raise ValueError("Nonfinite DIAMOND output")
            if not (30 <= identity <= 100 and 50 <= qcover <= 100 and 50 <= scover <= 100
                    and 0 <= evalue <= 1e-5 and bitscore >= 0):
                raise ValueError(f"Hit violates frozen thresholds at {path}:{number}")
            matched.add(query)
            best_sources[labels[dataset]] += 1
    matched_pillars = {gene_to_reference[g] for g in matched if g in gene_to_reference}
    totals = Counter(query_species.values())
    by_species = Counter(query_species[g] for g in matched)
    return {
        "query_proteins": len(query_species), "matching_query_proteins": len(matched),
        "matching_query_fraction": len(matched) / len(query_species),
        "reference_pillars": len(references), "reference_pillars_with_hit": len(matched_pillars),
        "reference_pillar_fraction_with_hit": len(matched_pillars) / len(references),
        "best_hit_dataset_counts": dict(best_sources),
        "by_species": {s: {"proteins": totals[s], "with_hit": by_species[s],
                           "fraction_with_hit": by_species[s] / totals[s]} for s in sorted(totals)},
        "reference_pillars_with_hit_ids": sorted(matched_pillars),
        "limitations": [
            "A qualifying hit is a homology-overlap screen, not proof of orthology or reference-label reuse.",
            "No qualifying hit does not establish absence of remote homology or independence.",
            "Only the best qualifying hit is retained; dataset counts are not exhaustive per-dataset overlap counts.",
            "No method accuracy results or outcome-based family selection are used by this screen.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--development", action="append", required=True, help="NAME=DIRECTORY")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--diamond", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=32)
    args = parser.parse_args()
    if args.threads < 1:
        raise ValueError("Threads must be positive")
    version = subprocess.check_output([str(args.diamond), "version"], text=True).strip()
    if version != "diamond version 2.1.11":
        raise ValueError(f"Unexpected DIAMOND version: {version}")
    if args.output.exists():
        raise ValueError(f"Refusing to overwrite screen: {args.output}")
    manifest_path = args.prepared / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    for item in manifest["inputs"] + [manifest["reference"]]:
        if file_provenance(Path(item["path"]))["sha256"] != item["sha256"]:
            raise ValueError(f"Changed prepared artifact: {item['path']}")
    development = {}
    for entry in args.development:
        label, path = entry.split("=", 1)
        if not label or label in development:
            raise ValueError("Development labels must be unique and nonempty")
        development[label] = Path(path)
    args.output.mkdir(parents=True)
    query = args.output / "queries.fasta"
    query_species = {}
    with query.open("w") as handle:
        for item in manifest["inputs"]:
            path = Path(item["path"])
            for record in SeqIO.parse(path, "fasta"):
                if record.id in query_species:
                    raise ValueError(f"Duplicate query ID: {record.id}")
                query_species[record.id] = path.stem
                SeqIO.write(record, handle, "fasta")
    target = args.output / "development.fasta"
    mapping = args.output / "target_ids.tsv"
    counts, sources = [], {}
    print("Preparing development sequences", flush=True)
    with target.open("w") as handle, mapping.open("w") as table:
        table.write("screen_id\tdataset\tfile\toriginal_id\n")
        for index, (label, directory) in enumerate(development.items()):
            paths = sorted(p for p in directory.iterdir() if p.suffix in {".fa", ".faa", ".fasta", ".fsa"})
            if not paths:
                raise ValueError(f"No FASTAs in {directory}")
            count = 0
            for path in paths:
                for record in SeqIO.parse(path, "fasta"):
                    count += 1
                    identifier = f"D{index}G{count}"
                    sequence = str(record.seq).upper().removesuffix("*")
                    SeqIO.write(SeqRecord(Seq(sequence), id=identifier, description=""), handle, "fasta")
                    table.write(f"{identifier}\t{label}\t{path.name}\t{record.id}\n")
            counts.append(count)
            sources[label] = [file_provenance(p) for p in paths]
    database = args.output / "development"
    hits = args.output / "hits.tsv"
    commands = [[str(args.diamond), "makedb", "--in", str(target), "--db", str(database), "--threads", str(args.threads)],
                search_command(args.diamond, query, database, hits, args.threads)]
    provenance = {"schema_version": 1, "started_at": datetime.now(timezone.utc).isoformat(),
                  "source": file_provenance(Path(__file__)), "manifest": file_provenance(manifest_path),
                  "development_inputs": sources, "commands": commands, "diamond_version": version,
                  "diamond_binary": file_provenance(args.diamond), "target_counts": counts,
                  "command": [sys.executable, *sys.argv]}
    (args.output / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    for phase, command in zip(("makedb", "search"), commands):
        print(f"Running {phase}", flush=True)
        with (args.output / (phase + ".log")).open("w") as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    references = json.loads(Path(manifest["reference"]["path"]).read_text())
    result = summarize_hits(hits, query_species, references, counts, list(development))
    result.update(status="complete", generated_at=datetime.now(timezone.utc).isoformat(),
                  provenance=provenance, hits=file_provenance(hits), target_mapping=file_provenance(mapping))
    (args.output / "summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(f"Complete: {result['matching_query_proteins']} of {result['query_proteins']} queries have qualifying hits", flush=True)


if __name__ == "__main__":
    main()
