"""Audit native OrthoFinder tables with memory bounded by one species pair."""

import csv
from itertools import combinations, product
from pathlib import Path

from benchmark_tools.orthofinder_to_pairwise import _accession
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


def read_table(path, a, b, owners, membership=None):
    pairs, rows, expanded = set(), 0, 0
    with path.open(newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        if next(reader, None) != ["Orthogroup", a, b]:
            raise ValueError("Wrong native table header or orientation")
        for row in reader:
            if len(row) != 3 or not row[0]:
                raise ValueError("Malformed native orthologue row")
            left, right = ([g.strip() for g in cell.split(",")] for cell in row[1:])
            if any(not g or owners.get(g) != a for g in left) or any(not g or owners.get(g) != b for g in right):
                raise ValueError("Unknown, empty or incorrectly owned native gene")
            if membership is not None:
                groups = {membership.get(g) for g in left + right}
                if None in groups or len(groups) != 1:
                    raise ValueError("Native relation crosses or escapes MCL checkpoint groups")
            rows += 1
            expanded += len(left) * len(right)
            pairs.update(tuple(sorted(pair)) for pair in product(left, right))
    return pairs, {"rows": rows, "expanded_relations": expanded,
                   "distinct_relations": len(pairs), "duplicate_relations": expanded - len(pairs)}


def audit_tables(results, owners, species, membership=None):
    species = sorted(species)
    if (len(species) < 2 or len(species) != len(set(species))
            or set(owners.values()) != set(species)
            or any(not s or Path(s).name != s for s in species)):
        raise ValueError("Invalid native species universe")
    if membership is not None and set(membership) != set(owners):
        raise ValueError("MCL checkpoint does not cover exact input universe")
    accessions = [_accession(g) for g in owners]
    if any(not a for a in accessions) or len(set(accessions)) != len(owners):
        raise ValueError("Accession normalization is not injective")
    directory = Path(results) / "Orthologues"
    expected = {directory / f"Orthologues_{a}" / f"{a}__v__{b}.tsv"
                for a in species for b in species if a != b}
    summary_paths = {directory / f"{s}.tsv" for s in species}
    actual = {p for p in directory.rglob("*.tsv") if p.is_file()}
    if actual - summary_paths != expected:
        raise ValueError("Incomplete or unexpected native species-pair table inventory")
    rows, checked = [], []
    for a, b in combinations(species, 2):
        forward = directory / f"Orthologues_{a}" / f"{a}__v__{b}.tsv"
        reverse = directory / f"Orthologues_{b}" / f"{b}__v__{a}.tsv"
        identities = [record(forward), record(reverse)]
        left, left_stats = read_table(forward, a, b, owners, membership)
        right, right_stats = read_table(reverse, b, a, owners, membership)
        if left != right:
            raise ValueError("Native table orientations disagree on ortholog relations")
        for item in identities:
            check(item)
        rows.append({"species": [a, b], "forward": left_stats, "reverse": right_stats})
        checked.extend(identities)
        del left, right
    for item in checked:
        check(item)
    if {p for p in directory.rglob("*.tsv") if p.is_file()} != actual:
        raise ValueError("Native table inventory changed during audit")
    return {"status": "native_orthofinder_pair_tables_verified", "species": len(species),
            "input_genes": len(owners), "directed_tables": len(checked), "species_pairs": rows,
            "distinct_pairs": sum(r["forward"]["distinct_relations"] for r in rows),
            "converter_emitted_pairs": sum(r["forward"]["expanded_relations"] for r in rows),
            "converter_duplicate_pairs": sum(r["forward"]["duplicate_relations"] for r in rows),
            "mcl_membership_checked": membership is not None, "checked_files": checked,
            "source": record(__file__), "accuracy_admitted": False,
            "limitations": ["Native table integrity only; caller must admit execution, inputs and species tree separately.",
                            "Duplicate native relations are counted, not silently removed or interpreted as independent observations.",
                            "Top-level per-species summary tables are not pair-converter inputs and are not semantically audited here.",
                            "Memory is bounded by two orientations of the largest species-pair relation set, not constant."]}
