import os
from typing import Tuple, List, Dict

from Bio import SeqIO
import numpy as np


def get_sequence_lengths(
    fasta_directory: str,
    files: List[str]
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
    temporary_directory: str,
) -> np.ndarray:
    res_path = f"{temporary_directory}/{taxon_a}_2_{taxon_b}.phmmerout.txt"

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


def determine_edge_thresholds(
    files: List[str],
    fasta_directory: str,
    temporary_directory: str,
) -> Tuple[np.ndarray, Dict[np.str_, np.float64], Dict[frozenset, np.float64]]:
    gene_lengths = get_sequence_lengths(fasta_directory, files)

    reciprocal_best_hit_thresholds = {}
    pairwise_rbh_corr = {}

    for file in files:
        file_pairs = [(file, i) for i in files]

        for pair in file_pairs:
            fwd_res = read_and_filter_phmmer_output(pair[0], pair[1], temporary_directory)
            rev_res = read_and_filter_phmmer_output(pair[1], pair[0], temporary_directory)

            fwd_res_merged = merge_with_gene_lengths(fwd_res, gene_lengths)
            rev_res_merged = merge_with_gene_lengths(rev_res, gene_lengths)

            fwd_res_merged = normalize_by_gene_length(fwd_res_merged)
            rev_res_merged = normalize_by_gene_length(rev_res_merged)

            # get dictionaries of scores for best hit and best hit
            best_hits_A_to_B = get_best_hits_and_scores(fwd_res_merged)
            best_hits_B_to_A = get_best_hits_and_scores(rev_res_merged)

            # TODO: continue refactoring here
            # get rbh scores
            rbh_scores = []
            rbh_pairs_identified = 0
            for query, best_hit in best_hits_A_to_B.items():
                target = best_hit['target']
                # Check if the reciprocal hit is also the best hit
                if target in best_hits_B_to_A and best_hits_B_to_A[target]['target'] == query:
                    score_between_rbh_pair = ((best_hit['score'] + best_hits_B_to_A[target]['score']) / 2)
                    rbh_scores.append(score_between_rbh_pair)
                    rbh_pairs_identified += 1

            if frozenset(pair) not in pairwise_rbh_corr:
                pairwise_rbh_corr[frozenset(pair)] = np.mean(rbh_scores)
            else:
                pairwise_rbh_corr[frozenset(pair)] = (pairwise_rbh_corr[frozenset(pair)] + np.mean(rbh_scores)) / 2

            # Phylogenetic correction
            best_hit_scores_A_to_B = {key: value["score"] / pairwise_rbh_corr[frozenset(pair)] for key, value in best_hits_A_to_B.items()}
            best_hit_scores_B_to_A = {key: value["score"] / pairwise_rbh_corr[frozenset(pair)] for key, value in best_hits_B_to_A.items()}

            # get lower bound threshold, which is sequence length and distance corrected
            for geneA, geneB in best_hits_A_to_B.items():
                geneB = geneB["target"]
                if best_hits_B_to_A.get(geneB)["target"] == geneA:
                    score = ((best_hit_scores_A_to_B[geneA] + best_hit_scores_B_to_A[geneB]) / 2)
                    if geneA in reciprocal_best_hit_thresholds:
                        if score < reciprocal_best_hit_thresholds[geneA]:
                            reciprocal_best_hit_thresholds[geneA] = score
                    else:
                        reciprocal_best_hit_thresholds[geneA] = score

    return gene_lengths, reciprocal_best_hit_thresholds, pairwise_rbh_corr


def determine_network_edges(
    files: List[str],
    temporary_directory: str,
    gene_lengths: np.ndarray,
    pairwise_rbh_corr: Dict[frozenset, np.float64],
    reciprocal_best_hit_thresholds: Dict[np.str_, np.float64],
) -> Dict[frozenset, np.float64]:
    edges = dict()
    gene_lengths = {str(row['name']): int(row['length']) for row in gene_lengths}
    for file in files:

        file_pairs = [(file, i) for i in files]

        for pair in file_pairs:
            res = f"{temporary_directory}/{pair[0]}_2_{pair[1]}.phmmerout.txt"

            dtype_res = [("target_name", "U50"), ("query_name", "U50"), ("evalue", float), ("score", float)]
            res = np.genfromtxt(res, comments="#", dtype=dtype_res, usecols=[0, 2, 4, 5], encoding='utf-8')

            res = res[res["evalue"] < 0.0001]

            for hit in res:
                query_length = gene_lengths[hit['query_name']]
                target_length = gene_lengths[hit['target_name']]
                norm_score = (hit['score'] / (query_length + target_length)) / pairwise_rbh_corr[frozenset(pair)]
                try:
                    if norm_score >= reciprocal_best_hit_thresholds[hit["query_name"]]:
                        genes = frozenset([hit['query_name'], hit['target_name']])
                        if len(genes) == 2:
                            if genes in edges:
                                if edges[genes] <= norm_score:
                                    edges[genes] = norm_score
                            else:
                                edges[genes] = norm_score
                # except reciprocal_best_hit_thresholds[hit["query_name"]] doesn't exist
                except KeyError:
                    continue

    with open(f"{temporary_directory}/orthohmm_edges.txt", "w") as file:
        for key, value in edges.items():
            key_str = "\t".join(map(str, key))
            file.write(f"{key_str}\t{value}\n")

    return edges


def generate_orthogroup_files(
    output_directory: str,
    gene_lengths: np.ndarray,
    files: list,
    fasta_directory: str,
    single_copy_threshold: float,
    extensions: Tuple[str],
    temporary_directory: str,
) -> Tuple[list, list, dict]:
    clustering_res = []
    with open(
        f"{temporary_directory}/orthohmm_edges_clustered.txt",
        "r"
    ) as file:
        for line in file:
            line = line.strip()
            if line:
                clustering_res.append(line.split())

    # get singletons - genes that aren't in groups with other genes
    singletons = list(set(gene_lengths["name"]) - set([j for i in clustering_res for j in i]))
    singletons = [[x] for x in singletons]
    clustering_res.extend(singletons)

    total_ogs = len(clustering_res)
    width = len(str(total_ogs))  # Determine the width needed for the largest number
    total_number_of_taxa = len(np.unique(gene_lengths["spp"]))

    # get all FASTA entries
    entries = {}
    for fasta_file in files:
        fasta_file_entries = {}

        for record in SeqIO.parse(f"{fasta_directory}/{fasta_file}", "fasta"):
            # Store gene ID and sequence in the dictionary
            fasta_file_entries[record.id] = str(record.seq)

        # Add the file_entries dictionary to the all_entries dictionary
        entries[fasta_file] = fasta_file_entries

    ogs_dat = dict()
    og_cn = dict()
    single_copy_ogs = list()
    og_cn["files:"] = files
    for i in range(total_ogs):
        # Format the OG number with leading zeros
        og_id = f"OG{i:0{width}}:"
        og_rows = gene_lengths[np.isin(gene_lengths['name'], clustering_res[i])]
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
                # cnts.append(str(spp_counts[spp_counts["spp"] == file]["cnt"].values[0]))
                cnts.append(str(spp_counts.get(file, 0)))
            except IndexError:
                cnts.append("0")
        og_cn[og_id] = cnts

        clustering_res[i].insert(0, og_id)

    # write clusters file
    with open(f"{output_directory}/orthohmm_orthogroups.txt", 'w') as file:
        for cluster in clustering_res:
            # Join the elements of the sublist with spaces and write to the file
            genes_in_cluster = cluster[1:]
            genes_in_cluster.sort()
            file.write(cluster[0] + " " + " ".join(map(str, genes_in_cluster)) + "\n")

    # write copy number per species
    with open(f"{output_directory}/orthohmm_gene_count.txt", 'w') as file:
        for key, value in og_cn.items():
            # Write each key-value pair as a space-delimited string
            file.write(f"{key} {' '.join(value)}\n")

    # write list of single-copy orthologs
    with open(f"{output_directory}/orthohmm_single_copy_orthogroups.txt", 'w') as file:
        for key in og_cn.keys():
            if key != "files:":
                file.write(f"{key[:-1]}\n")

    # write all og files
    if not os.path.isdir(f"{output_directory}/orthohmm_orthogroups"):
        os.mkdir(f"{output_directory}/orthohmm_orthogroups")
    for og_id, fasta_dat in ogs_dat.items():
        with open(f"{output_directory}/orthohmm_orthogroups/{og_id}.fa", "w") as file:
            file.write("\n".join(fasta_dat)+"\n")

    # write single-copy og files
    if not os.path.isdir(f"{output_directory}/orthohmm_single_copy_orthogroups"):
        os.mkdir(f"{output_directory}/orthohmm_single_copy_orthogroups")
    for single_copy_og in single_copy_ogs:
        for idx in range(len(ogs_dat[single_copy_og])):
            if ogs_dat[single_copy_og][idx][0] == ">":
                taxon_name = gene_lengths[gene_lengths["name"] == ogs_dat[single_copy_og][idx][1:]]["spp"][0]
                # remove extension
                taxon_name = next((taxon_name[:-len(ext)] for ext in extensions if taxon_name.endswith(ext)), taxon_name)
                ogs_dat[single_copy_og][idx] = f">{taxon_name}|{ogs_dat[single_copy_og][idx][1:]}"
        with open(f"{output_directory}/orthohmm_single_copy_orthogroups/{single_copy_og}.fa", "w") as file:
            file.write("\n".join(ogs_dat[single_copy_og]) + "\n")

    return single_copy_ogs, singletons, ogs_dat
