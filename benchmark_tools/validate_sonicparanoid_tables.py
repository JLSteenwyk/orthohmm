"""Validate native SonicParanoid table coverage and accession ownership."""

from itertools import combinations
import math
from pathlib import Path

from benchmark_tools.fastoma_to_pairwise import input_owners
from benchmark_tools.sonicparanoid_to_pairwise import EXPECTED_HEADER, _accession, _members, iter_pairs


def validate_tables(directory, inputs):
    directory = Path(directory)
    inputs = [Path(p) for p in inputs]
    owners = input_owners(inputs)
    species = {p.name: p.stem for p in inputs}
    expected = {tuple(sorted(pair)) for pair in combinations(species, 2)}
    seen = set()
    rows = relations = unique = 0
    stats = {"duplicates": 0}
    files = sorted(p for p in directory.rglob("*") if p.is_file())
    for path in files:
        left = path.parent.name
        prefix = left + "-"
        if path.parent.parent != directory or left not in species or not path.name.startswith(prefix):
            raise ValueError("Unexpected species-pair path: " + str(path))
        right = path.name[len(prefix):]
        key = tuple(sorted((left, right)))
        if right not in species or key not in expected or key in seen:
            raise ValueError("Invalid/duplicate species-pair table: " + str(path))
        seen.add(key)
        with path.open(encoding="ascii") as stream:
            if next(stream, "").rstrip("\n").split("\t") != EXPECTED_HEADER:
                raise ValueError("Unexpected SonicParanoid header")
            for number, line in enumerate(stream, 2):
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 4:
                    raise ValueError("Malformed SonicParanoid row")
                sizes = []
                for field, name in zip(fields[2:], (left, right)):
                    members = _members(field, path, number)
                    accessions = [_accession(identifier) for identifier in members]
                    if len(set(accessions)) != len(accessions):
                        raise ValueError("Duplicate accession aliases in one member list")
                    if any(owners.get(identifier) != species[name] for identifier in accessions):
                        raise ValueError("Accession ownership mismatch")
                    if any(not math.isfinite(float(value)) for value in field.split()[1::2]):
                        raise ValueError("Nonfinite confidence")
                    sizes.append(len(members))
                if int(fields[0]) != sum(sizes) or int(fields[1]) != sizes[0] * sizes[1]:
                    raise ValueError("Native row counts differ")
                rows += 1
                relations += int(fields[1])
        unique += sum(1 for _ in iter_pairs(path, stats))
    if seen != expected:
        raise ValueError("Incomplete input-species matrix")
    if relations != unique + stats["duplicates"]:
        raise ValueError("Relation accounting differs")
    return {"input_species": len(species), "input_accessions": len(owners),
            "species_pair_tables": len(files), "native_rows": rows,
            "raw_relations": relations, "distinct_pairs": unique,
            "duplicate_relations": stats["duplicates"]}
