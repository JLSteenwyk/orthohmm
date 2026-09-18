"""Check native SonicParanoid group tables against inputs and normalized groups."""

import argparse
import csv
import json
from pathlib import Path
import sys

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_three_kingdoms_pair_counts import partition, record
from benchmark_tools.audit_three_kingdoms_method_inputs import snapshot_rows


def validate(table, inputs, snapshot, normalized):
    inputs = [Path(p) for p in inputs]
    evidence = [record(p) for p in [table, snapshot, normalized, *inputs]]
    if len(inputs) < 2 or len({p.name for p in inputs}) != len(inputs):
        raise ValueError("Require distinct input species")
    owners, counts = {}, {}
    for path in inputs:
        count = 0
        for sequence in SeqIO.parse(path, "fasta"):
            if not sequence.id or not sequence.seq or sequence.id in owners:
                raise ValueError("Empty sequence or duplicate input identifier")
            owners[sequence.id] = path.name
            count += 1
        if not count:
            raise ValueError("Empty species input")
        counts[path.name] = count
    snap = snapshot_rows(Path(snapshot).read_text())
    if set(snap) != set(counts):
        raise ValueError("Snapshot species mismatch")
    for item in evidence[3:]:
        row = snap[Path(item["path"]).name]
        if row["sha256"] != item["sha256"] or row["proteins"] != counts[Path(item["path"]).name]:
            raise ValueError("Snapshot input identity/count mismatch")
    norm, sizes = partition(normalized)
    seen, group_ids, normalized_groups = set(), set(), set()
    with Path(table).open(newline="") as stream:
        rows = csv.reader(stream, delimiter="\t")
        header = next(rows, [])
        if (header[:4] != ["group_id", "group_size", "sp_in_grp", "seed_ortholog_cnt"]
                or len(header[4:]) != len(counts) or set(header[4:]) != set(counts)):
            raise ValueError("Unexpected native species columns")
        for row in rows:
            if len(row) != len(header) or not row[0] or row[0] in group_ids:
                raise ValueError("Malformed row or duplicate native group")
            group_ids.add(row[0])
            genes, occupied = [], 0
            for species, cell in zip(header[4:], row[4:]):
                if cell in ("", "*"):
                    continue
                occupied += 1
                members = cell.split(",")
                if any(not gene or gene.strip() != gene or owners.get(gene) != species for gene in members):
                    raise ValueError("Invalid member or species ownership")
                genes.extend(members)
            if not genes or len(set(genes)) != len(genes) or seen.intersection(genes):
                raise ValueError("Empty group or duplicate gene membership")
            seen.update(genes)
            try:
                size, species_count, seeds = map(int, row[1:4])
            except ValueError as exc:
                raise ValueError("Invalid native counts") from exc
            if size != len(genes) or species_count != occupied or not 0 <= seeds <= size:
                raise ValueError("Native count mismatch")
            indices = {norm.get(gene) for gene in genes}
            if len(indices) != 1 or None in indices:
                raise ValueError("Normalized membership differs")
            index = next(iter(indices))
            if index in normalized_groups or sizes[index] != len(genes):
                raise ValueError("Normalized partition differs")
            normalized_groups.add(index)
    if not group_ids or len(normalized_groups) != len(sizes) or seen != set(norm):
        raise ValueError("Normalized/native coverage mismatch")
    for item in evidence:
        if record(item["path"]) != item:
            raise ValueError("Input changed during validation")
    return {"status": "sonic_native_group_conversion_verified", "source": record(__file__),
            "input_species": len(counts), "input_genes": len(owners), "groups": len(group_ids),
            "grouped_genes": len(seen), "unassigned_input_genes": len(owners) - len(seen),
            "evidence": evidence, "accuracy_admitted": False,
            "limitations": ["Does not establish terminal execution or runtime provenance.",
                            "Snapshot residue totals retained by source, not independently validated here.",
                            "Unassigned input genes are counted, not assumed to be missing output errors.",
                            "Group co-membership is not native pair-table orthology."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", required=True, type=Path)
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--snapshot", required=True, type=Path)
    parser.add_argument("--normalized", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = validate(args.table, sorted(args.input_dir.glob("*.fasta")), args.snapshot, args.normalized)
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps({key: report[key] for key in ("status", "groups", "grouped_genes", "unassigned_input_genes")}))
