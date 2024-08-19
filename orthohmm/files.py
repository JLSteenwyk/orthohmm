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
    if not os.path.isdir(f"{output_directory}/orthohmm_orthogroups"):
        os.mkdir(f"{output_directory}/orthohmm_orthogroups")
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
    extensions: Tuple,
) -> None:
    # write single-copy og files
    if not os.path.isdir(f"{output_directory}/orthohmm_single_copy_orthogroups"):
        os.mkdir(f"{output_directory}/orthohmm_single_copy_orthogroups")
    for single_copy_og in single_copy_ogs:
        for idx in range(len(ogs_dat[single_copy_og])):
            if ogs_dat[single_copy_og][idx][0] == ">":
                taxon_name = \
                    gene_lengths[gene_lengths["name"] == ogs_dat[single_copy_og][idx][1:]]["spp"][0]
                # remove extension
                taxon_name = next((taxon_name[:-len(ext)] for ext in extensions if taxon_name.endswith(ext)), taxon_name)
                ogs_dat[single_copy_og][idx] = f">{taxon_name}|{ogs_dat[single_copy_og][idx][1:]}"
        with open(f"{output_directory}/orthohmm_single_copy_orthogroups/{single_copy_og}.fa", "w") as file:
            file.write("\n".join(ogs_dat[single_copy_og]) + "\n")
