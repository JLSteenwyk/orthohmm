import itertools
import multiprocessing
import os
import subprocess
import sys
import time
from typing import Tuple

from Bio import SeqIO
import pandas as pd


def run_bash_commands(
    command: list
) -> None:
    subprocess.run(command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def run_phmmer(
    phmmer_cmds: list,
    cpu: int,
) -> None:
    # Create a pool of workers
    pool = multiprocessing.Pool(processes=cpu)

    # Map the commands to the worker pool
    pool.map(run_bash_commands, phmmer_cmds)

    # Close the pool and wait for the work to finish
    pool.close()
    pool.join()


def calculate_sequence_lengths(
    fasta_directory: str,
    file: str,
) -> list:
    sequence_lengths = []
    with open(f"{fasta_directory}/{file}", 'r') as fasta_file:
        sequence_id = None
        sequence_length = 0

        for line in fasta_file:
            line = line.strip()
            if line.startswith('>'):
                if sequence_id:
                    sequence_lengths.append([file, sequence_id, sequence_length])
                sequence_id = line[1:]
                sequence_length = 0
            else:
                sequence_length += len(line)

        if sequence_id:
            sequence_lengths.append([file, sequence_id, sequence_length])

    return sequence_lengths


def get_best_hits(
    df: pd.DataFrame
) -> pd.DataFrame:
    best_hits = df.loc[df.groupby("query name")["norm_score"].idxmax()]
    return best_hits.set_index("query name")["target name"].to_dict(), \
        best_hits.set_index("query name")["norm_score"].to_dict()


def execute_phmmer_search(
        files: list,
        cpu: int,
        output_directory: str,
        fasta_directory: str,
        phmmer: str,
) -> None:
    # create phmmer cmds
    pairwise_combos = list(itertools.product(files, repeat=2))
    phmmer_cmds = []
    for combo in pairwise_combos:
        phmmer_cmds.append(f"{phmmer} --noali --notextw --cpu {cpu} --tblout {output_directory}/working_dir/{combo[0]}_2_{combo[1]}.phmmerout.txt {fasta_directory}/{combo[0]} {fasta_directory}/{combo[1]}")

    # run phmmer cmds
    run_phmmer(phmmer_cmds, cpu)


def determine_edge_thresholds(
        files: list,
        fasta_directory: str,
        output_directory: str,
        col_names: list,
        columns_to_drop: list,
) -> Tuple[pd.DataFrame, dict, dict]:
    # get sequence lengths
    gene_lengths = []
    for file in files:
        gene_lengths.extend(calculate_sequence_lengths(fasta_directory, file))
    gene_lengths = pd.DataFrame(gene_lengths, columns=["spp", "query name", "length"]) 

    reciprocal_best_hit_thresholds = {}  # for saving thresholds
    pairwise_rbh_corr = {}  # for phylogenetic correction
    for file in files:
        file_pairs = [(file, i) for i in files]

        for pair in file_pairs:
            fwd_res = f"{output_directory}/working_dir/{pair[0]}_2_{pair[1]}.phmmerout.txt"
            rev_res = f"{output_directory}/working_dir/{pair[1]}_2_{pair[0]}.phmmerout.txt"

            fwd_res = pd.read_csv(fwd_res, comment="#", sep="\s+", header=None)
            rev_res = pd.read_csv(rev_res, comment="#", sep="\s+", header=None)

            fwd_res.columns = col_names
            rev_res.columns = col_names

            fwd_res = fwd_res[fwd_res["E-value"] < 0.0001]
            rev_res = rev_res[rev_res["E-value"] < 0.0001]

            gene_lengths.columns = ["spp", "query name", "query length"]
            fwd_res = pd.merge(fwd_res, gene_lengths, on="query name")
            rev_res = pd.merge(rev_res, gene_lengths, on="query name")
            gene_lengths.columns = ["spp", "target name", "target length"]
            fwd_res = pd.merge(fwd_res, gene_lengths, on="target name")
            rev_res = pd.merge(rev_res, gene_lengths, on="target name")

            fwd_res = fwd_res.drop(columns=columns_to_drop)
            rev_res = rev_res.drop(columns=columns_to_drop)

            fwd_res["norm_score"] = fwd_res["score"]/(fwd_res["query length"]+fwd_res["target length"])
            rev_res["norm_score"] = rev_res["score"]/(rev_res["query length"]+rev_res["target length"])

            best_hits_A_to_B, best_hit_scores_A_to_B = get_best_hits(fwd_res)
            best_hits_B_to_A, best_hit_scores_B_to_A = get_best_hits(rev_res)

            # get lower threshold per gene
            # first, get lower threshold and normalize by sequence length and phylo dist
            rbh_scores = []
            rbh_pairs_identified = 0
            for geneA, geneB in best_hits_A_to_B.items():
                if best_hits_B_to_A.get(geneB) == geneA:
                    score_between_rbh_pair = ((best_hit_scores_A_to_B[geneA] + best_hit_scores_B_to_A[geneB]) / 2)
                    rbh_scores.append(score_between_rbh_pair)
                    rbh_pairs_identified += 1
            if frozenset(pair) not in pairwise_rbh_corr:
                pairwise_rbh_corr[frozenset(pair)] = sum(rbh_scores) / len(rbh_scores) 
            else:
                pairwise_rbh_corr[frozenset(pair)] = (pairwise_rbh_corr[frozenset(pair)] + (sum(rbh_scores) / len(rbh_scores))) / 2

            # phylo correction
            best_hit_scores_A_to_B = {key: value / pairwise_rbh_corr[frozenset(pair)] for key, value in best_hit_scores_A_to_B.items()}
            best_hit_scores_B_to_A = {key: value / pairwise_rbh_corr[frozenset(pair)] for key, value in best_hit_scores_B_to_A.items()}

            # get lower bound threshold, which is sequence length and distance corrected
            for geneA, geneB in best_hits_A_to_B.items():
                if best_hits_B_to_A.get(geneB) == geneA:
                    score = ((best_hit_scores_A_to_B[geneA] + best_hit_scores_B_to_A[geneB]) / 2)
                    if geneA in reciprocal_best_hit_thresholds:
                        if score < reciprocal_best_hit_thresholds[geneA]:
                            reciprocal_best_hit_thresholds[geneA] = score
                    else:
                        reciprocal_best_hit_thresholds[geneA] = score

    return gene_lengths, reciprocal_best_hit_thresholds, pairwise_rbh_corr


def determine_network_edges(
    files: list,
    output_directory: str,
    col_names: list,
    columns_to_drop: list,
    gene_lengths: list,
    pairwise_rbh_corr: dict,
    reciprocal_best_hit_thresholds: dict,
) -> dict:
    edges = dict()
    for file in files:

        file_pairs = [(file, i) for i in files]

        for pair in file_pairs:
            res = f"{output_directory}/working_dir/{pair[0]}_2_{pair[1]}.phmmerout.txt"

            res = pd.read_csv(res, comment="#", sep="\s+", header=None)

            res.columns = col_names

            gene_lengths.columns = ["spp", "query name", "query length"]
            res = pd.merge(res, gene_lengths, on="query name")
            gene_lengths.columns = ["spp", "target name", "target length"]
            res = pd.merge(res, gene_lengths, on="target name")

            res = res[res["E-value"] < 0.0001]

            res = res.drop(columns=columns_to_drop)

            res["norm_score"] = res["score"] / (res["query length"] + res["target length"])
            res["norm_score"] = res["norm_score"] / pairwise_rbh_corr[frozenset(pair)]

            res = res.values.tolist()

            for hit in res:
                try:
                    if hit[7] >= reciprocal_best_hit_thresholds[hit[1]]:
                        genes = frozenset(hit[0:2])
                        if len(genes) == 2:
                            if genes in edges:
                                if edges[genes] <= hit[7]:
                                    edges[genes] = hit[7]
                            else:
                                edges[genes] = hit[7]
                # except reciprocal_best_hit_thresholds[hit[1]] doesn't exist
                except KeyError:
                    continue

    with open(f"{output_directory}/working_dir/orthohmm_edges.txt", "w") as file:
        for key, value in edges.items():
            key_str = "\t".join(map(str, key))
            file.write(f"{key_str}\t{value}\n")

    return edges


def execute_mcl(
    mcl: str,
    inflation_value: float,
    output_directory: str,
) -> None:
    cmd = f"{mcl} {output_directory}/working_dir/orthohmm_edges.txt --abc -I {inflation_value} -o {output_directory}/working_dir/orthohmm_edges_clustered.txt"
    subprocess.run(
        cmd,
        shell=True,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )


def generate_orthogroup_files(
    output_directory: str,
    gene_lengths: pd.DataFrame,
    files: list,
    fasta_directory: str,
    single_copy_threshold: float,
    extensions: tuple,
) -> Tuple[list, list, dict]:
    clustering_res = []
    with open(
        f"{output_directory}/working_dir/orthohmm_edges_clustered.txt",
        "r"
    ) as file:
        for line in file:
            line = line.strip()
            if line:
                clustering_res.append(line.split())

    # get singletons - genes that aren't in groups with other genes
    singletons = list(set(gene_lengths["target name"]) - set([j for i in clustering_res for j in i]))
    singletons = [[x] for x in singletons]
    clustering_res.extend(singletons)

    total_ogs = len(clustering_res)
    width = len(str(total_ogs))  # Determine the width needed for the largest number
    total_number_of_taxa = gene_lengths["spp"].nunique()

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
        og_rows = gene_lengths[gene_lengths["target name"].isin(clustering_res[i])]
        # test if single-copy
        if og_rows["spp"].nunique() == len(og_rows["spp"]):
            # test if sufficient occupancy
            if og_rows["spp"].nunique() / total_number_of_taxa > single_copy_threshold:
                single_copy_ogs.append(f"OG{i}")
        og_dat = list()
        for row in og_rows.itertuples(index=True):
            lines = [entries[row.spp][row._2][i:i+70] for i in range(0, len(entries[row.spp][row._2]), 70)]
            og_dat.append(f">{row._2}")
            og_dat.extend(lines)
        ogs_dat[f"OG{i}"] = og_dat

        counts = og_rows["spp"].value_counts()
        counts_df = counts.reset_index()
        counts_df.columns = ["spp", "cnt"]

        cnts = []
        for file in files:
            try:
                cnts.append(str(counts_df[counts_df["spp"] == file]["cnt"].values[0]))
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
                taxon_name = gene_lengths.loc[gene_lengths["target name"] == ogs_dat[single_copy_og][idx][1:], "spp"].values[0]
                # remove extension
                taxon_name = next((taxon_name[:-len(ext)] for ext in extensions if taxon_name.endswith(ext)), taxon_name)
                ogs_dat[single_copy_og][idx] = f">{taxon_name}|{ogs_dat[single_copy_og][idx][1:]}"
        with open(f"{output_directory}/orthohmm_single_copy_orthogroups/{single_copy_og}.fa", "w") as file:
            file.write("\n".join(ogs_dat[single_copy_og]) + "\n")

    return single_copy_ogs, singletons, ogs_dat
