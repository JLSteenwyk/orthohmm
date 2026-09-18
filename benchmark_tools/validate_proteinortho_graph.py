"""Validate native post-clustering graph ownership against complete input FASTAs."""

import math
from itertools import combinations
from pathlib import Path

from benchmark_tools.fastoma_to_pairwise import input_owners
from benchmark_tools.proteinortho_to_pairwise import _accession, iter_pairs


def validate_graph(path, inputs):
    """Require all input-species sections, valid ownership, and finite scores.

    Empty sections are allowed: no ortholog edge is required for a species pair.
    This checks representation, not biological accuracy or workflow completion.
    """
    path = Path(path)
    if not path.name.endswith(".proteinortho-graph"):
        raise ValueError("Require post-clustering .proteinortho-graph, not search graph")
    inputs = [Path(p) for p in inputs]
    owners = input_owners(inputs)
    names = {p.name: p.stem for p in inputs}
    expected = {tuple(sorted(pair)) for pair in combinations(names, 2)}
    seen = set()
    section = None
    count = 0
    with path.open(encoding="ascii") as stream:
        for number, line in enumerate(stream, 1):
            if line.startswith("#"):
                fields = line[1:].strip().split("\t")
                if len(fields) == 2 and fields != ["file_a", "file_b"]:
                    if any(name not in names for name in fields):
                        raise ValueError(f"Unknown species header at row {number}")
                    key = tuple(sorted(fields))
                    if key not in expected or key in seen:
                        raise ValueError(f"Invalid/duplicate species section at row {number}")
                    seen.add(key)
                    section = tuple(names[name] for name in fields)
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 6 or section is None:
                raise ValueError(f"Malformed relation or missing section at row {number}")
            for identifier, species in zip(fields[:2], section):
                if not identifier or any(c.isspace() for c in identifier):
                    raise ValueError(f"Malformed identifier at row {number}")
                if owners.get(_accession(identifier)) != species:
                    raise ValueError(f"Accession ownership mismatch at row {number}")
            if any(not math.isfinite(float(value)) for value in fields[2:]):
                raise ValueError(f"Nonfinite score at row {number}")
            count += 1
    if seen != expected:
        raise ValueError("Incomplete input-species matrix")
    # Preserve the existing converter's duplicate and structural checks as well.
    if sum(1 for _ in iter_pairs(path)) != count:
        raise ValueError("Converter row count differs")
    return {"input_species": len(names), "input_accessions": len(owners),
            "species_sections": len(seen), "validated_pair_rows": count}
