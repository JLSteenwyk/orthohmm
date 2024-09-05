import os
from typing import (
    Dict, List, Tuple
)

import numpy as np


def write_clusters_file(
    output_directory: str,
    clustering_res: List[List[str]],
) -> None:
    with open(f"{output_directory}/orthohmm_orthogroups.txt", 'w') as file:
        for cluster in clustering_res:
            genes_in_cluster = cluster[1:]
            genes_in_cluster.sort()
            file.write(
                cluster[0] + " " + " ".join(map(str, genes_in_cluster)) + "\n"
            )


def write_copy_number_file(
    output_directory: str,
    og_cn: Dict[str, List[str]]
) -> None:

    with open(f"{output_directory}/orthohmm_gene_count.txt", 'w') as file:
        for key, value in og_cn.items():
            file.write(f"{key} {' '.join(value)}\n")


def write_file_of_single_copy_ortholog_names(
    output_directory: str,
    og_cn: Dict[str, List[str]]
) -> None:
    with open(
        f"{output_directory}/orthohmm_single_copy_orthogroups.txt", 'w'
    ) as file:
        for key in og_cn.keys():
            if key != "files:":
                file.write(f"{key[:-1]}\n")


def write_fasta_files_for_all_ogs(
    output_directory: str,
    ogs_dat: Dict[str, List[str]],
) -> None:
    output_dir = f"{output_directory}/orthohmm_orthogroups"
    os.makedirs(output_dir, exist_ok=True)

    for og_id, fasta_dat in ogs_dat.items():
        with open(
            f"{output_directory}/orthohmm_orthogroups/{og_id}.fa", "w"
        ) as file:
            file.write("\n".join(fasta_dat)+"\n")


def write_fasta_files_for_single_copy_orthologs(
    output_directory: str,
    ogs_dat: Dict[str, List[str]],
    gene_lengths: np.ndarray,
    single_copy_ogs: List[str],
    extensions: Tuple
) -> None:
    output_dir = f"{output_directory}/orthohmm_single_copy_orthogroups"
    os.makedirs(output_dir, exist_ok=True)

    name_to_species = {row["name"]: row["spp"] for row in gene_lengths}

    for single_copy_og in single_copy_ogs:
        updated_entries = []
        for entry in ogs_dat[single_copy_og]:
            if entry.startswith(">"):
                gene_name = entry[1:]
                # Find the species, handling the absence case
                species = name_to_species.get(gene_name, "UnknownSpecies")
                taxon_name = next((species[:-len(ext)] for ext in extensions if species.endswith(ext)), species)
                new_header = f">{taxon_name}|{gene_name}"
                updated_entries.append(new_header)
            else:
                updated_entries.append(entry)

        filepath = os.path.join(output_dir, f"{single_copy_og}.fa")
        with open(filepath, "w") as file:
            file.write("\n".join(updated_entries) + "\n")
