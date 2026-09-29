"""Inventory a historical derived TreeFam archive without assigning QfO labels."""

import argparse
import csv
import json
from pathlib import Path

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def rows(path):
    with path.open(newline="") as handle:
        return [row for row in csv.reader(handle, delimiter="\t")
                if row and not row[0].startswith("#")]


def inspect(directory):
    download_path = directory / "download.json"
    download = json.loads(download_path.read_text())
    for artifact in download["artifacts"]:
        check(artifact["file"])
    genes = {}
    for label in ("human", "fly", "worm"):
        values = rows(directory / f"tf.{label}.genes.raw")
        if any(len(r) != 8 or r[-1] != "" or not r[0].isdigit() or not r[1] for r in values):
            raise ValueError("Unexpected historical gene-table format")
        if len({r[1] for r in values}) != len(values):
            raise ValueError("Duplicate source gene ID")
        genes[label] = values
    memberships = {}
    catalogs = {label: {r[1] for r in values} for label, values in genes.items()}
    all_ids = set.union(*catalogs.values())
    for suffix in ("raw", "all"):
        groups = {}
        for row in rows(directory / f"tf.human_fly_worm.fam_genes.{suffix}"):
            if len(row) < 2 or not row[0].startswith("TF") or row[0] in groups or any(not v for v in row[1:]):
                raise ValueError("Invalid or repeated family record")
            groups[row[0]] = row[1:]
        ids = {v for values in groups.values() for v in values}
        memberships[suffix] = dict(families=len(groups), membership_records=sum(map(len, groups.values())),
            unique_ids=len(ids), ids_absent_from_downloaded_gene_tables=len(ids - all_ids),
            duplicate_ids_within_families=sum(len(v) - len(set(v)) for v in groups.values()))
    result = dict(download=record(download_path), source=record(Path(__file__)),
        inputs=[a["file"] for a in download["artifacts"]],
        gene_table_rows={label: len(values) for label, values in genes.items()},
        fly_and_worm_data_rows_identical=genes["fly"] == genes["worm"],
        fly_and_worm_id_intersection=len(catalogs["fly"] & catalogs["worm"]),
        memberships=memberships, original_trees_recovered=False, qfo_mapping_recovered=False,
        scientific_use_admitted=False,
        limitations=["Derived archive labels do not establish the species or original family membership.",
                     "HTTP download hashes provide byte identity, not authenticated source provenance.",
                     "No topology, duplication events or QfO mapping was reconstructed."])
    for artifact in download["artifacts"]:
        check(artifact["file"])
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = inspect(args.directory)
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps({k: result[k] for k in ("gene_table_rows", "fly_and_worm_data_rows_identical", "memberships")}))
