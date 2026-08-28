#!/usr/bin/env python

import logging
import os
import sys
import time
from typing import Union

from .args_processing import process_args
from .externals import (
    execute_leiden,
    execute_mcl,
    execute_phmmer_search,
)
from .files import fetch_fasta_files
from .helpers import (
    determine_edge_thresholds,
    determine_network_edges,
    generate_orthogroup_clusters_file,
    generate_orthogroup_files,
    generate_phmmer_cmds,
    StartStep,
    StopStep,
    SubstitutionMatrix,
)
from .parser import create_parser
from .refinement import (
    DEFAULT_COPY_SPLIT_MIN_DATASET_SPECIES,
    resolve_refinement_profile,
    refine_cluster_indices,
)
from .writer import (
    write_user_args,
    write_output_stats
)

logger = logging.getLogger(__name__)


def _collect_search_hit_arrays(search_results, evalue_threshold, gene_to_id):
    hit_queries = []
    hit_targets = []
    hit_scores = []
    for result in search_results.values():
        for query, target, evalue, score in zip(
            result.query_names,
            result.target_names,
            result.evalues,
            result.scores,
        ):
            if evalue < evalue_threshold:
                query_id = gene_to_id.get(str(query))
                target_id = gene_to_id.get(str(target))
                if query_id is None or target_id is None:
                    continue
                hit_queries.append(query_id)
                hit_targets.append(target_id)
                hit_scores.append(float(score))
    return hit_queries, hit_targets, hit_scores


def _collect_edge_arrays(edges, gene_to_id):
    edge_queries = []
    edge_targets = []
    edge_scores = []
    for edge, score in edges.items():
        if len(edge) != 2:
            continue
        gene_a, gene_b = tuple(edge)
        gene_a_id = gene_to_id.get(str(gene_a))
        gene_b_id = gene_to_id.get(str(gene_b))
        if gene_a_id is None or gene_b_id is None:
            continue
        edge_queries.append(gene_a_id)
        edge_targets.append(gene_b_id)
        edge_scores.append(float(score))
    return edge_queries, edge_targets, edge_scores


def _refine_cluster_file(output_directory, gene_lengths, search_results, edges,
                         evalue_threshold, refinement_profile="default", **kwargs):
    if search_results is None:
        return

    cluster_path = f"{output_directory}/orthohmm_working_res/orthohmm_edges_clustered.txt"
    gene_names = [str(row["name"]) for row in gene_lengths]
    gene_to_id = {gene: idx for idx, gene in enumerate(gene_names)}
    species_to_id = {}
    gene_to_species = []
    for row in gene_lengths:
        species = str(row["spp"])
        species_to_id.setdefault(species, len(species_to_id))
        gene_to_species.append(species_to_id[species])
    broad_copy_only = (
        len(species_to_id) >= DEFAULT_COPY_SPLIT_MIN_DATASET_SPECIES
    )

    clusters = []
    with open(cluster_path, "r") as handle:
        for line in handle:
            genes = line.strip().split()
            if genes:
                cluster_ids = [
                    gene_to_id[gene]
                    for gene in genes
                    if gene in gene_to_id
                ]
                if cluster_ids:
                    clusters.append(cluster_ids)

    if broad_copy_only:
        hit_queries = []
        hit_targets = []
        hit_scores = []
        edge_queries, edge_targets, edge_scores = _collect_edge_arrays(
            edges,
            gene_to_id,
        )
    else:
        hit_queries, hit_targets, hit_scores = _collect_search_hit_arrays(
            search_results,
            evalue_threshold,
            gene_to_id,
        )
        edge_queries, edge_targets, edge_scores = _collect_edge_arrays(edges, gene_to_id)
    refinement_kwargs = resolve_refinement_profile(refinement_profile)
    refinement_kwargs.update(kwargs)
    refined = refine_cluster_indices(
        clusters,
        hit_queries,
        hit_targets,
        hit_scores,
        edge_queries,
        edge_targets,
        gene_to_species,
        rbnh_scores=edge_scores,
        **refinement_kwargs,
    )

    with open(cluster_path, "w") as handle:
        for cluster in refined:
            if cluster:
                handle.write(" ".join(sorted(gene_names[i] for i in cluster)) + "\n")


def execute(
    fasta_directory: str,
    output_directory: str,
    phmmer: str,
    cpu: int,
    single_copy_threshold: float,
    mcl: str,
    inflation_value: float,
    start: Union[StartStep, None],
    stop: Union[StopStep, None],
    substitution_matrix: SubstitutionMatrix,
    evalue_threshold: float,
    search_mode: str = "builtin",
    clustering: str = "leiden",
    cpm_resolution=0.1,
    refinement_profile: str = "default",
    **kwargs,
) -> None:
    # for reporting runtime duration to user
    start_time = time.time()

    # make working dir
    working_dir = f"{output_directory}/orthohmm_working_res/"
    os.makedirs(working_dir, exist_ok=True)

    files = fetch_fasta_files(fasta_directory)

    search_results = None

    if start != StartStep.search_res and search_mode == "phmmer":
        phmmer_cmds = generate_phmmer_cmds(
            files,
            phmmer,
            output_directory,
            fasta_directory,
            cpu,
            substitution_matrix,
        )

    # print phmmer cmds and exit is users only want to prepare phmmer cmds
    if stop == StopStep.prepare:
        if search_mode == "builtin":
            print("--stop prepare is only supported with --search_mode phmmer")
            sys.exit(1)
        for cmd in phmmer_cmds:
            print(cmd)
        sys.exit()

    # display to user what args are being used in stdout
    write_user_args(
        fasta_directory,
        output_directory,
        phmmer,
        mcl,
        cpu,
        single_copy_threshold,
        files,
        start,
        stop,
        substitution_matrix,
        evalue_threshold,
        inflation_value,
        search_mode,
        clustering,
        cpm_resolution,
    )

    # set current step and determine the total number of
    # steps that will be used in the run
    current_step = 1
    if stop == StopStep.infer:
        total_steps = 4
    else:
        total_steps = 5

    if start == StartStep.search_res:
        total_steps -= 1
    elif search_mode == "builtin":
        print(f"Step {current_step}/{total_steps}: Conducting all-to-all comparisons (built-in search).")
        from .search.engine import execute_builtin_search
        search_results = execute_builtin_search(
            files,
            fasta_directory,
            output_directory,
            cpu,
            substitution_matrix,
            evalue_threshold=evalue_threshold,
        )
        print("\r          Completed!      \n")
        current_step += 1
    else:
        print(f"Step {current_step}/{total_steps}: Conducting all-to-all comparisons.")
        execute_phmmer_search(
            phmmer_cmds,
            cpu,
        )
        print("\r          Completed!      \n")
        current_step += 1

    print(f"Step {current_step}/{total_steps}: Determining edge thresholds")
    gene_lengths, reciprocal_best_hit_thresholds, pairwise_rbh_corr = \
        determine_edge_thresholds(
            files,
            fasta_directory,
            output_directory,
            cpu,
            evalue_threshold,
            search_results=search_results,
        )
    print("\r          Completed!      \n")
    current_step += 1

    print(f"Step {current_step}/{total_steps}: Identifying network edges")
    edges = determine_network_edges(
        files,
        output_directory,
        gene_lengths,
        pairwise_rbh_corr,
        reciprocal_best_hit_thresholds,
        evalue_threshold,
        cpu,
        search_results=search_results,
    )
    print("\r          Completed!      \n")
    current_step += 1

    print(f"Step {current_step}/{total_steps}: Conducting clustering")
    if clustering == "mcl":
        execute_mcl(
            mcl,
            inflation_value,
            cpu,
            output_directory,
        )
    else:
        execute_leiden(
            cpm_resolution,
            output_directory,
        )
        _refine_cluster_file(
            output_directory,
            gene_lengths,
            search_results,
            edges,
            evalue_threshold,
            refinement_profile=refinement_profile,
        )
    singletons, og_cn, ogs_dat, single_copy_ogs = \
        generate_orthogroup_clusters_file(
            output_directory,
            gene_lengths,
            files,
            single_copy_threshold,
            fasta_directory,
        )
    print("          Completed!\n")
    current_step += 1

    # exit if users only want orthogroups to be inferred
    if stop == StopStep.infer:
        sys.exit()

    print(f"Step {current_step}/{total_steps}: Writing orthogroup information")
    generate_orthogroup_files(
        output_directory,
        gene_lengths,
        og_cn,
        ogs_dat,
        single_copy_ogs,
    )
    print("          Completed!\n")

    write_output_stats(
        start_time,
        single_copy_ogs,
        singletons,
        ogs_dat,
        edges,
        gene_lengths,
    )


def main(argv=None):
    """
    Function that parses and collects arguments
    """
    parser = create_parser()
    args = parser.parse_args()

    execute(**process_args(args))


if __name__ == "__main__":
    main(sys.argv[1:])
