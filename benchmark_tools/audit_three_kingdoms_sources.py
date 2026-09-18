"""Trace retained Three Kingdoms files without downloading or changing inputs."""

import argparse
import configparser
import csv
import gzip
import hashlib
import json
from pathlib import Path
import shlex
import sys

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

DOWNLOAD_SHA = "e27bc332a3709547f9298bbd908c6b376921a40723305d210bd40932b1475b25"
REFERENCE_SHA = "a5f3447056ecfa305442caff0d898d524eed350f13587d3247d3e1ba4c19757d"
CODES = {"Amph", "Anol", "Dani", "Xeno", "Arab", "Oryz", "Sela", "Zeam", "Asfu", "Scer", "Scom", "Spom"}


def source_rows(text):
    rows = {}
    # Read quoted pipe-delimited array records with a shell lexer, never execute the script.
    for token in shlex.split(text, comments=True):
        fields = token.split("|")
        if len(fields) != 3 or fields[1] not in CODES:
            continue
        species, code, url = fields
        for prefix, host in (("${EBI}", "https://ftp.ebi.ac.uk"), ("${ENS}", "https://ftp.ensembl.org")):
            if url.startswith(prefix):
                url = host + url[len(prefix):]
        if code in rows or not url.startswith("https://") or "$" in url:
            raise ValueError("Ambiguous or unsupported source record")
        rows[code] = {"species": species, "code": code, "source_url": url,
                      "moving_release_url": "/current_release/" in url}
    if set(rows) != CODES:
        raise ValueError("Incomplete source inventory")
    return rows


def sequence_inventory(path):
    genes = {}
    for sequence in SeqIO.parse(path, "fasta"):
        if not sequence.id or sequence.id in genes or not sequence.seq:
            raise ValueError("Duplicate or empty FASTA record")
        data = str(sequence.seq).encode("ascii")
        genes[sequence.id] = {"length": len(data), "sha256": hashlib.sha256(data).hexdigest()}
    if not genes:
        raise ValueError("Empty proteome")
    return genes


def decompressed_identity(path):
    digest, size = hashlib.sha256(), 0
    with gzip.open(path, "rb") as stream:
        while block := stream.read(1024 * 1024):
            digest.update(block)
            size += len(block)
    return {"bytes": size, "sha256": digest.hexdigest()}


def change_details(raw_path, staged_path, changed, reference_genes):
    if not changed:
        return []
    raw = {r.id: str(r.seq) for r in SeqIO.parse(raw_path, "fasta") if r.id in changed}
    staged = {r.id: str(r.seq) for r in SeqIO.parse(staged_path, "fasta") if r.id in changed}
    return [{"gene": gene, "raw_length": len(raw[gene]), "staged_length": len(staged[gene]),
             "raw_stop_markers": raw[gene].count("*"), "staged_stop_markers": staged[gene].count("*"),
             "exactly_explained_by_removing_stop_markers": raw[gene].replace("*", "") == staged[gene],
             "in_scored_reference": gene in reference_genes} for gene in sorted(changed)]


def complete_hits(path, genes):
    result = {}
    lines = path.read_text().splitlines()
    headers = [line for line in lines if line.startswith("#")]
    expected = ["# BUSCO version is: 5.8.2", "# The lineage dataset is: eukaryota_odb10 (Creation date: 2024-01-08, number of genomes: 70, number of BUSCOs: 255)"]
    if [line.strip() for line in headers[:2]] != expected:
        raise ValueError("BUSCO version or lineage differs")
    for row in csv.reader((line for line in lines if line and not line.startswith("#")), delimiter="\t"):
        if len(row) < 2 or row[1] not in {"Complete", "Duplicated", "Fragmented", "Missing"}:
            raise ValueError("Malformed BUSCO row")
        if row[1] != "Complete":
            continue
        if len(row) < 3:
            raise ValueError("Malformed Complete BUSCO row")
        bid, _, gene = row[:3]
        if bid in result or gene not in genes:
            raise ValueError("Duplicate Complete family or unknown gene")
        result[bid] = gene
    return result


def audit(root):
    script = root / "download_proteomes.sh"
    source = record(script)
    if source["sha256"] != DOWNLOAD_SHA:
        raise ValueError("Download-intent script changed")
    sources = source_rows(script.read_text())
    reference = record(root / "busco/reference_orthogroups.txt")
    if reference["sha256"] != REFERENCE_SHA:
        raise ValueError("Retained reference differs from scored reference")
    reference_genes = set((root / "busco/reference_orthogroups.txt").read_text().split())
    cfg = root / "busco_downloads/lineages/eukaryota_odb10/dataset.cfg"
    config_record = record(cfg)
    config = configparser.ConfigParser()
    config.read_string("[dataset]\n" + cfg.read_text())
    values = dict(config["dataset"])
    required = {"name": "eukaryota_odb10", "creation_date": "2024-01-08",
                "number_of_buscos": "255", "number_of_species": "70", "orthodb_version": "10.1"}
    if any(values.get(k) != v for k, v in required.items()):
        raise ValueError("Unexpected lineage configuration")
    records = [source, reference, config_record, record(root / "run_busco.sh"), record(root / "build_busco_reference_ogs.py")]
    rows, families, universe = [], {}, set()
    for code, metadata in sorted(sources.items()):
        paths = {"compressed": root / f"proteomes/{code}.fa.gz", "raw": root / f"input.raw/{code}.fasta",
                 "staged": root / f"input/{code}.fasta", "busco_table": root / f"busco/{code}/run_eukaryota_odb10/full_table.tsv"}
        files = {key: record(path) for key, path in paths.items()}
        records.extend(files.values())
        decompressed = decompressed_identity(paths["compressed"])
        if any(decompressed[k] != files["raw"][k] for k in ("bytes", "sha256")):
            raise ValueError("Compressed download does not reproduce retained raw FASTA")
        raw, staged = sequence_inventory(paths["raw"]), sequence_inventory(paths["staged"])
        if universe.intersection(staged):
            raise ValueError("Gene identifiers collide across species")
        universe.update(staged)
        complete = complete_hits(paths["busco_table"], staged)
        for bid, gene in complete.items():
            families.setdefault(bid, []).append(gene)
        shared = set(raw) & set(staged)
        changed = {g for g in shared if raw[g] != staged[g]}
        rows.append({**metadata, "files": files, "decompressed": decompressed,
                     "raw_proteins": len(raw), "staged_proteins": len(staged),
                     "added_ids": sorted(set(staged) - set(raw)), "removed_ids": sorted(set(raw) - set(staged)),
                     "changed_sequence_ids": sorted(changed),
                     "sequence_changes": change_details(paths["raw"], paths["staged"], changed, reference_genes),
                     "raw_staged_byte_match": files["raw"]["sha256"] == files["staged"]["sha256"],
                     "complete_busco_families": len(complete)})
    groups = [sorted(genes) for _, genes in sorted(families.items()) if len(genes) >= 2]
    reconstructed = "".join(" ".join(group) + "\n" for group in groups).encode()
    if reconstructed != (root / "busco/reference_orthogroups.txt").read_bytes():
        raise ValueError("Complete-hit tables do not reproduce scored reference")
    counts = {"proteomes": len(rows), "proteins": len(universe), "reference_groups": len(groups),
              "reference_genes": len({g for group in groups for g in group}),
              "reference_pairs": sum(len(g) * (len(g) - 1) // 2 for g in groups)}
    if counts != {"proteomes": 12, "proteins": 443217, "reference_groups": 255, "reference_genes": 2035, "reference_pairs": 7352}:
        raise ValueError("Dataset counts differ from retained publication panel")
    for item in records:
        check(item)
    return {"status": "retained_three_kingdoms_lineage_and_reference_verified", "publication_ready": False,
            "source": record(__file__), "counts": counts, "lineage": values, "inputs": rows,
            "supporting_files": records[:5], "reference_reconstruction_byte_match": True,
            "limitations": ["URLs are retained download intent, not independently authenticated upstream archive identity.",
                "No upstream archive was downloaded again, and a moving current_release URL is not reproducible by itself.",
                "Sequence comparisons preserve case and content; header/line-wrap differences are separate from byte identity.",
                "Does not prove original BUSCO execution environment/input provenance or validate every BUSCO classification.",
                "Does not grant redistribution rights; provider, lineage and release-specific terms require separate review.",
                "Reference evaluates conserved single-copy complete hits, not proteome-wide false positives."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.root.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps(result["counts"], sort_keys=True))
