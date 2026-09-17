"""Strict native group readers preserving species columns and missing assignments."""

import csv

from benchmark_tools.score_wgd_application import membership
from benchmark_tools.score_ygob_groups import read_predictions

HEADERS = {
    "orthofinder_root_hogs": ("HOG", "OG", "Gene Tree Parent Clade"),
    "sonicparanoid": ("group_id", "group_size", "sp_in_grp", "seed_ortholog_cnt"),
}


def read_species_table(path, format, owners, column_species):
    if format not in HEADERS:
        raise ValueError("Unknown species-table format")
    metadata = HEADERS[format]
    groups = {}
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        columns = header[len(metadata):]
        if (tuple(header[:len(metadata)]) != metadata or not columns
                or len(set(header)) != len(header) or set(columns) != set(column_species)
                or len(set(column_species.values())) != len(column_species)
                or set(column_species.values()) != set(owners.values())):
            raise ValueError("Unexpected native species columns or metadata")
        for row in reader:
            if len(row) != len(header):
                raise ValueError("Truncated or extra native table fields")
            name = row[0].strip()
            if not name or name in groups:
                raise ValueError("Empty or duplicate native group ID")
            genes, occupied = [], 0
            for column, cell in zip(columns, row[len(metadata):]):
                cell = cell.strip()
                if not cell or (format == "sonicparanoid" and cell == "*"):
                    continue
                members = [g.strip() for g in cell.split(",")]
                if any(not g or owners.get(g) != column_species[column] for g in members):
                    raise ValueError("Foreign gene, empty member or wrong species column")
                genes.extend(members)
                occupied += 1
            if format == "sonicparanoid":
                try:
                    size, species_count, seeds = map(int, row[1:4])
                except ValueError as error:
                    raise ValueError("Invalid SonicParanoid counts") from error
                if size != len(genes) or species_count != occupied or not 0 <= seeds <= size:
                    raise ValueError("SonicParanoid metadata counts differ")
            groups[name] = genes
    membership(groups, owners)
    if not groups:
        raise ValueError("Empty native group table")
    return groups


def read_orthohmm(path, format, owners):
    groups = read_predictions(path, format)
    membership(groups, owners)
    if not groups:
        raise ValueError("Empty native group table")
    return groups
