"""Strictly stream native FastOMA pairs; preserve order and multiplicity."""

import argparse
import gzip
from pathlib import Path
import sys

from Bio import SeqIO


def input_owners(paths):
    owners = {}
    species = set()
    for path in paths:
        path = Path(path)
        if path.stem in species:
            raise ValueError("Duplicate species input")
        species.add(path.stem)
        count = 0
        for entry in SeqIO.parse(path, "fasta"):
            fields = entry.id.split("|")
            identifier = fields[1] if len(fields) > 1 else entry.id
            if not identifier or identifier in owners or not entry.seq:
                raise ValueError("Empty sequence/accession or duplicate accession: " + identifier)
            owners[identifier] = path.stem
            count += 1
        if not count:
            raise ValueError("Empty species input: " + str(path))
    if len(species) < 2:
        raise ValueError("Require at least two species inputs")
    return owners


def iter_pairs(path, owners):
    """Reject malformed/foreign/intraspecies rows; do not silently drop any row."""
    opener = gzip.open if Path(path).suffix == ".gz" else open
    with opener(path, "rt", encoding="ascii", newline="") as stream:
        for number, line in enumerate(stream, 1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 2 or any(not value or any(c.isspace() for c in value) for value in fields):
                raise ValueError(f"Malformed FastOMA pair at {path}:{number}")
            first, second = fields
            if first not in owners or second not in owners:
                raise ValueError(f"Unknown FastOMA accession at {path}:{number}")
            if first == second or owners[first] == owners[second]:
                raise ValueError(f"Self/intraspecies FastOMA pair at {path}:{number}")
            yield tuple(sorted((first, second)))


def write_pairs(path, owners, output):
    count = 0
    for first, second in iter_pairs(path, owners):
        output.write(f"{first}\t{second}\n")
        count += 1
    if not count:
        raise ValueError("Empty FastOMA pair output")
    return count


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pairs", type=Path)
    parser.add_argument("input_directory", type=Path)
    args = parser.parse_args()
    owners = input_owners(sorted(args.input_directory.glob("*.fasta")))
    count = write_pairs(args.pairs, owners, sys.stdout)
    print(f"Validated {count} native FastOMA pair rows; order/multiplicity retained, uniqueness not asserted", file=sys.stderr)
