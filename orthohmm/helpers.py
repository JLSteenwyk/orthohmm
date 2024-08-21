from enum import Enum
import itertools
from typing import Tuple, List, Dict

import numpy as np

from .files import (
    write_fasta_files_for_all_ogs,
    write_fasta_files_for_single_copy_orthologs,
    write_file_of_single_copy_ortholog_names,
    write_clusters_file,
    write_copy_number_file,
)


class StopStep(Enum):
    prepare = "prepare"
    infer = "infer"
    write = "write"


class StartStep(Enum):
    search_res = "search_res"


def generate_phmmer_cmds(
    files: List[str],
    phmmer: str,
    output_directory: str,
    fasta_directory: str,
    cpu: int,
    stop: str,
):
    pairwise_combos = list(itertools.product(files, repeat=2))
    phmmer_cmds = []
    for combo in pairwise_combos:
        if stop == "prepare":
            phmmer_cmds.append(f"{phmmer} --noali --notextw --cpu {cpu} --tblout {output_directory}/orthohmm_working_res/{combo[0]}_2_{combo[1]}.phmmerout.txt {fasta_directory}/{combo[0]} {fasta_directory}/{combo[1]}")
        else:
            phmmer_cmds.append(f"{phmmer} --noali --notextw --tblout {output_directory}/orthohmm_working_res/{combo[0]}_2_{combo[1]}.phmmerout.txt {fasta_directory}/{combo[0]} {fasta_directory}/{combo[1]}")

    return phmmer_cmds


def get_sequence_lengths(
    fasta_directory: str,
    files: List[str],
) -> np.ndarray:
    # Get sequence lengths
    gene_lengths = []
    for file in files:
        with open(f"{fasta_directory}/{file}", 'r') as fasta_file:
            sequence_id = None
            sequence_length = 0

            for line in fasta_file:
                line = line.strip()
                if line.startswith('>'):
                    if sequence_id:
                        gene_lengths.append([file, sequence_id, sequence_length])
                    sequence_id = line[1:]
                    sequence_length = 0
                else:
                    sequence_length += len(line)

            if sequence_id:
                gene_lengths.append([file, sequence_id, sequence_length])
    # make sure types are correct
    gene_lengths = [(str(row[0]), str(row[1]), int(row[2])) for row in gene_lengths]
    dtype = [("spp", "U50"), ("name", "U50"), ("length", int)]

    return np.array(gene_lengths, dtype=dtype)


def merge_with_gene_lengths(
    res: np.ndarray,
    gene_lengths: np.ndarray,
) -> np.ndarray:
    res_merged = np.empty([len(res), 6], dtype=object)
    for idx in range(len(res)):
        temp = [None, None, None, None, None, None]
        for spp_gene_length in gene_lengths:
            if spp_gene_length[1] == res[idx][0]:
                temp[0] = res[idx][0]
                temp[1] = res[idx][1]
                temp[2] = res[idx][2]
                temp[3] = res[idx][3]
                temp[4] = spp_gene_length[2]

            if spp_gene_length[1] == res[idx][1]:
                temp[5] = spp_gene_length[2]

            res_merged[idx] = temp
    return res_merged


def read_and_filter_phmmer_output(
    taxon_a: str,
    taxon_b: str,
    output_directory: str,
) -> np.ndarray:
    res_path = f"{output_directory}/orthohmm_working_res/{taxon_a}_2_{taxon_b}.phmmerout.txt"

    dtype_res = [
        ("target_name", "U50"),
        ("query_name", "U50"),
        ("evalue", float),
        ("score", float)
    ]

    res = np.genfromtxt(res_path, comments="#", dtype=dtype_res, usecols=[0, 2, 4, 5], encoding='utf-8')

    res = res[res["evalue"] < 0.0001]

    return res


def normalize_by_gene_length(res_merged: np.ndarray) -> np.ndarray:
    res_merged[:, 3] = res_merged[:, 3] / (res_merged[:, 4] + res_merged[:, 5])

    return res_merged


def correct_by_phylogenetic_distance(
    best_hits_A_to_B: Dict[np.str_, Dict[np.str_, np.float64]],
    best_hits_B_to_A: Dict[np.str_, Dict[np.str_, np.float64]],
    pair: Tuple[str, str],
    pairwise_rbh_corr: Dict[frozenset, np.float64],
) -> Tuple[
        Dict[np.str_, np.float64],
        Dict[np.str_, np.float64],
        Dict[frozenset, np.float64],
]:
    # get rbh scores
    rbh_scores = []
    rbh_pairs_identified = 0
    for query, best_hit in best_hits_A_to_B.items():
        target = best_hit['target']
        # Check if the reciprocal hit is also the best hit
        if target in best_hits_B_to_A and best_hits_B_to_A[target]["target"] == query:
            score_between_rbh_pair = (
                (best_hit["score"] + best_hits_B_to_A[target]["score"]) / 2
            )
            rbh_scores.append(score_between_rbh_pair)
            rbh_pairs_identified += 1

    if frozenset(pair) not in pairwise_rbh_corr:
        pairwise_rbh_corr[frozenset(pair)] = np.mean(rbh_scores)
    else:
        pairwise_rbh_corr[frozenset(pair)] = (pairwise_rbh_corr[frozenset(pair)] + np.mean(rbh_scores)) / 2

    # Phylogenetic correction
    best_hit_scores_A_to_B = {key: value["score"] / pairwise_rbh_corr[frozenset(pair)] for key, value in best_hits_A_to_B.items()}
    best_hit_scores_B_to_A = {key: value["score"] / pairwise_rbh_corr[frozenset(pair)] for key, value in best_hits_B_to_A.items()}

    return best_hit_scores_A_to_B, best_hit_scores_B_to_A, pairwise_rbh_corr


def get_best_hits_and_scores(
    res_merged: np.ndarray
) -> Dict[np.str_, Dict[np.str_, np.float64]]:
    """
    get dictionaries of scores for best hit and best hit
    """
    best_hits = dict()

    for record in res_merged:
        query = record[1]
        target = record[0]
        score = record[3]
        if query not in best_hits or best_hits[query]["score"] < score:
            best_hits[query] = {"target": target, "score": score}

    return best_hits


def get_threshold_per_gene(
    best_hits_A_to_B: Dict[np.str_, Dict[np.str_, np.float64]],
    best_hits_B_to_A: Dict[np.str_, Dict[np.str_, np.float64]],
    best_hit_scores_A_to_B: Dict[np.str_, np.float64],
    best_hit_scores_B_to_A: Dict[np.str_, np.float64],
    reciprocal_best_hit_thresholds: Dict[np.str_, np.float64],
):
    for geneA, geneB in best_hits_A_to_B.items():
        geneB = geneB["target"]
        if best_hits_B_to_A.get(geneB)["target"] == geneA:
            score = ((best_hit_scores_A_to_B[geneA] + best_hit_scores_B_to_A[geneB]) / 2)
            if geneA in reciprocal_best_hit_thresholds:
                if score < reciprocal_best_hit_thresholds[geneA]:
                    reciprocal_best_hit_thresholds[geneA] = score
            else:
                reciprocal_best_hit_thresholds[geneA] = score

    return reciprocal_best_hit_thresholds


def determine_edge_thresholds(
    files: List[str],
    fasta_directory: str,
    output_directory: str,
) -> Tuple[
        np.ndarray,
        Dict[np.str_, np.float64],
        Dict[frozenset, np.float64],
]:
    pairwise_rbh_corr = dict()
    reciprocal_best_hit_thresholds = dict()

    gene_lengths = get_sequence_lengths(fasta_directory, files)

    for file in files:
        file_pairs = [(file, i) for i in files]

        for pair in file_pairs:
            fwd_res = read_and_filter_phmmer_output(
                pair[0], pair[1],
                output_directory
            )
            rev_res = read_and_filter_phmmer_output(
                pair[1], pair[0],
                output_directory
            )

            fwd_res_merged = merge_with_gene_lengths(fwd_res, gene_lengths)
            rev_res_merged = merge_with_gene_lengths(rev_res, gene_lengths)

            fwd_res_merged = normalize_by_gene_length(fwd_res_merged)
            rev_res_merged = normalize_by_gene_length(rev_res_merged)

            best_hits_A_to_B = get_best_hits_and_scores(fwd_res_merged)
            best_hits_B_to_A = get_best_hits_and_scores(rev_res_merged)

            best_hit_scores_A_to_B, best_hit_scores_B_to_A, pairwise_rbh_corr = \
                correct_by_phylogenetic_distance(
                    best_hits_A_to_B,
                    best_hits_B_to_A,
                    pair,
                    pairwise_rbh_corr
                )

            reciprocal_best_hit_thresholds = \
                get_threshold_per_gene(
                    best_hits_A_to_B, best_hits_B_to_A,
                    best_hit_scores_A_to_B, best_hit_scores_B_to_A,
                    reciprocal_best_hit_thresholds
                )

    return gene_lengths, reciprocal_best_hit_thresholds, pairwise_rbh_corr


def determine_network_edges(
    files: List[str],
    output_directory: str,
    gene_lengths: np.ndarray,
    pairwise_rbh_corr: Dict[frozenset, np.float64],
    reciprocal_best_hit_thresholds: Dict[np.str_, np.float64],
) -> Dict[frozenset, np.float64]:
    edges = dict()

    gene_lengths = {
        str(row["name"]): int(row["length"]) for row in gene_lengths
    }

    for file in files:

        file_pairs = [(file, i) for i in files]

        for pair in file_pairs:
            res = read_and_filter_phmmer_output(
                pair[0], pair[1],
                output_directory
            )

            for hit in res:
                query_length = gene_lengths[hit["query_name"]]
                target_length = gene_lengths[hit["target_name"]]

                norm_score = (
                    hit['score'] / (query_length + target_length)
                ) / pairwise_rbh_corr[frozenset(pair)]

                try:
                    if norm_score >= reciprocal_best_hit_thresholds[hit["query_name"]]:
                        genes = frozenset(
                            [hit["query_name"], hit["target_name"]]
                        )

                        if len(genes) == 2:
                            if genes in edges:
                                if edges[genes] <= norm_score:
                                    edges[genes] = norm_score
                            else:
                                edges[genes] = norm_score

                # except reciprocal_best_hit_thresholds[hit["query_name"]] doesn't exist
                except KeyError:
                    continue

    with open(f"{output_directory}/orthohmm_working_res/orthohmm_edges.txt", "w") as file:
        for key, value in edges.items():
            key_str = "\t".join(map(str, key))
            file.write(f"{key_str}\t{value}\n")

    return edges


def get_singletons(
    gene_lengths: np.ndarray,
    clustering_res: List[List[str]],
) -> Tuple[List[List[str]], List[List[str]]]:
    singletons = list(
        set(gene_lengths["name"]) - set([j for i in clustering_res for j in i])
    )
    singletons = [[str(i)] for i in singletons]
    clustering_res.extend(singletons)

    return clustering_res, singletons


def get_all_fasta_entries(
    fasta_directory: str,
    files: List[str]
) -> Dict[str, str]:
    entries = {}
    for fasta_file in files:
        fasta_file_entries = {}

        header = None
        sequence = []

        with open(f"{fasta_directory}/{fasta_file}", "r") as file:
            for line in file:
                line = line.strip()
                if line.startswith(">"):
                    if header:
                        fasta_file_entries[header] = ''.join(sequence)
                    header = line[1:]  # Remove the '>' character
                    sequence = []
                else:
                    sequence.append(line)

            # Don't forget to add the last entry to the dictionary
            if header:
                fasta_file_entries[header] = ''.join(sequence)

        entries[fasta_file] = fasta_file_entries

    return entries


def get_orthogroup_information(
    files: List[str],
    gene_lengths: np.ndarray,
    clustering_res: List[List[str]],
    single_copy_threshold: float,
    entries: Dict[str, str],
) -> Tuple[
    List[List[str]],
    Dict[str, List[str]],
    Dict[str, List[str]],
    List[str],
]:
    ogs_dat = dict()
    og_cn = dict()
    single_copy_ogs = list()

    total_ogs = len(clustering_res)
    width = len(str(total_ogs))
    total_number_of_taxa = len(np.unique(gene_lengths["spp"]))

    og_cn["files:"] = files

    for i in range(total_ogs):
        # Format the OG number with leading zeros
        og_id = f"OG{i:0{width}}:"
        og_rows = gene_lengths[
            np.isin(gene_lengths["name"], clustering_res[i])
        ]
        # test if single-copy
        if len(np.unique(og_rows["spp"])) == len(og_rows["spp"]):
            # test if sufficient occupancy
            if len(np.unique(og_rows["spp"])) / total_number_of_taxa > single_copy_threshold:
                single_copy_ogs.append(f"OG{i}")
        og_dat = list()
        for row in og_rows:
            lines = [entries[row["spp"]][row["name"]][i:i+70] for i in range(0, len(entries[row["spp"]][row["name"]]), 70)]
            og_dat.append(f">{row['name']}")
            og_dat.extend(lines)
        ogs_dat[f"OG{i}"] = og_dat

        spp_values, counts = np.unique(og_rows["spp"], return_counts=True)
        spp_counts = dict(zip(spp_values, counts))

        cnts = []
        for file in files:
            try:
                cnts.append(str(spp_counts.get(file, 0)))
            except IndexError:
                cnts.append("0")
        og_cn[og_id] = cnts

        clustering_res[i].insert(0, og_id)

    return clustering_res, og_cn, ogs_dat, single_copy_ogs


def generate_orthogroup_clusters_file(
    output_directory: str,
    gene_lengths: np.ndarray,
    files: List[str],
    single_copy_threshold: float,
    fasta_directory: str,
) -> Tuple[
    List[List[str]],
    Dict[str, List[str]],
    Dict[str, List[str]],
    List[str],
]:
    clustering_res = list()

    entries = get_all_fasta_entries(fasta_directory, files)

    with open(
        f"{output_directory}/orthohmm_working_res/orthohmm_edges_clustered.txt",
        "r",
    ) as file:
        for line in file:
            line = line.strip()
            if line:
                clustering_res.append(line.split())

    # get singletons - i.e., genes that aren't in groups with other genes
    clustering_res, singletons = get_singletons(
        gene_lengths,
        clustering_res,
    )

    clustering_res, og_cn, ogs_dat, single_copy_ogs = \
        get_orthogroup_information(
            files,
            gene_lengths,
            clustering_res,
            single_copy_threshold,
            entries,
        )

    write_clusters_file(output_directory, clustering_res)

    return singletons, og_cn, ogs_dat, single_copy_ogs


def generate_orthogroup_files(
    output_directory: str,
    gene_lengths: np.ndarray,
    extensions: Tuple[str],
    og_cn: Dict[str, List[str]],
    ogs_dat: Dict[str, List[str]],
    single_copy_ogs: List[str],
) -> None:
    write_copy_number_file(output_directory, og_cn)
    write_file_of_single_copy_ortholog_names(output_directory, og_cn)
    write_fasta_files_for_all_ogs(output_directory, ogs_dat)
    write_fasta_files_for_single_copy_orthologs(
        output_directory, ogs_dat, gene_lengths, single_copy_ogs, extensions
    )
