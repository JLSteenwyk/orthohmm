"""Prepare label-blind DIAMOND searches for the OrthoBench search-engine control."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_historical_profile_ablation import verify_file
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.run_simulation_methods import read_frozen

PREPARED_SHA = "5c325f4d77865e0c7571fe4bb4df0be0959977a49e1d22169f192518700f9382"


def search_command(binary, queries, database, output):
    return [str(binary), "blastp", "--query", str(queries), "--db", str(database),
            "--out", str(output), "--threads", "32", "--very-sensitive", "--evalue", "0.0001",
            "--matrix", "BLOSUM62", "--gapopen", "11", "--gapextend", "1",
            "--comp-based-stats", "1", "--masking", "1", "--max-target-seqs", "0", "--max-hsps", "1",
            "--outfmt", "6", "qseqid", "sseqid", "qlen", "slen", "score", "bitscore", "evalue"]


def write_queries(records, output):
    seen, owners, lengths = set(), {}, {}
    with output.open("x") as handle:
        for item in records:
            path = Path(item["path"])
            verify_file(path, item)
            for record in SeqIO.parse(path, "fasta"):
                if record.id in seen or not len(record.seq):
                    raise ValueError("Duplicate or empty input sequence")
                seen.add(record.id)
                owners[record.id] = path.name
                lengths[record.id] = len(record.seq)
                SeqIO.write(record, handle, "fasta")
    return owners, lengths


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--diamond", required=True, type=Path)
    args = parser.parse_args()
    root, output, binary = args.root.resolve(), args.output.resolve(), args.diamond.resolve()
    if output.exists():
        raise FileExistsError(output)
    manifest = root / "benchmark_tools/results/orthobench_factorial_prepared_20260916.json"
    prepared = read_frozen(manifest, PREPARED_SHA)
    version = subprocess.check_output([str(binary), "version"], text=True).strip()
    if version != "diamond version 2.1.11":
        raise ValueError("Unexpected DIAMOND version")
    records = sorted(prepared["fasta_inputs"], key=lambda r: r["path"])
    if len(records) != 12:
        raise ValueError("Expected frozen 12-species OrthoBench input panel")
    output.mkdir(parents=True)
    queries = output / "queries.fasta"
    owners, lengths = write_queries(records, queries)
    if len(owners) != 251378:
        raise ValueError("Unexpected input universe")
    metadata = output / "gene_metadata.json"
    metadata.write_text(json.dumps({g: {"species": owners[g], "length": lengths[g]} for g in sorted(owners)}, sort_keys=True) + "\n")
    searches = []
    for index, item in enumerate(records):
        directory = output / f"target_{index:02d}"
        directory.mkdir()
        database, hits = directory / "target", directory / "hits.tsv"
        searches.append({"index": index, "target_fasta": item,
                         "makedb": [str(binary), "makedb", "--in", item["path"], "--db", str(database), "--threads", "32"],
                         "search": search_command(binary, queries, database, hits), "output": str(hits)})
    report = {"schema_version": 1, "status": "prepared_not_searched", "accuracy_evaluated": False,
              "source": file_provenance(Path(__file__)), "frozen_factorial": file_provenance(manifest),
              "inputs": records, "queries": file_provenance(queries), "gene_metadata": file_provenance(metadata),
              "diamond": file_provenance(binary), "diamond_version": version, "searches": searches,
              "normalization": "DIAMOND raw score / sqrt(query_length * target_length), once only",
              "limits": "All reported targets; one HSP per query-target. Preserve self hits and directed asymmetry.",
              "scope": "Search-engine diagnostic with frozen downstream graph settings; sensitivity equivalence is not assumed.",
              "runtime_kind": "Search and database construction measured separately; shared-machine timings."}
    (output / "manifest.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(output / "manifest.json")


if __name__ == "__main__":
    main()
