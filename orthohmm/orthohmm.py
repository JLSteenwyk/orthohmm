#!/usr/bin/env python

import logging
import os
import sys
import time
from typing import Union

import numpy as np

from .accuracy import (
    build_rbnh_edges,
    build_singleton_assignment_edges,
    combine_edges,
    read_index_clusters,
    resolve_accuracy_profile,
    write_accuracy_checkpoint,
)
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
    get_sequence_lengths,
    IndexedEdges,
    StartStep,
    StopStep,
    SubstitutionMatrix,
    write_indexed_edges,
)
from .parser import create_parser
from .metrics import PipelineMetrics
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
    from .search.engine import IndexedSearchResults

    if isinstance(search_results, IndexedSearchResults):
        query_chunks = []
        target_chunks = []
        score_chunks = []
        species_global_ids = {
            species: np.fromiter(
                (gene_to_id[gene] for gene in names),
                dtype=np.int32,
                count=len(names),
            )
            for species, names in search_results.species_ids.items()
        }
        for (query_species, target_species), result in search_results.items():
            mask = result.evalues < evalue_threshold
            if not mask.any():
                continue
            query_chunks.append(
                species_global_ids[query_species][result.query_indices[mask]]
            )
            target_chunks.append(
                species_global_ids[target_species][result.target_indices[mask]]
            )
            score_chunks.append(result.scores[mask])
        if not query_chunks:
            return (
                np.empty(0, dtype=np.int32),
                np.empty(0, dtype=np.int32),
                np.empty(0, dtype=np.float64),
            )
        return (
            np.concatenate(query_chunks),
            np.concatenate(target_chunks),
            np.concatenate(score_chunks),
        )

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
    if isinstance(edges, IndexedEdges):
        return edges.sources, edges.targets, edges.weights

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
        return {}

    started = time.perf_counter()

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
    setup_finished = time.perf_counter()

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
    collection_finished = time.perf_counter()
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
    refinement_finished = time.perf_counter()

    with open(cluster_path, "w") as handle:
        for cluster in refined:
            if cluster:
                handle.write(" ".join(sorted(gene_names[i] for i in cluster)) + "\n")
    write_finished = time.perf_counter()
    return {
        "setup_and_cluster_read": round(setup_finished - started, 6),
        "hit_and_edge_collection": round(
            collection_finished - setup_finished, 6
        ),
        "numeric_refinement": round(
            refinement_finished - collection_finished, 6
        ),
        "cluster_write": round(write_finished - refinement_finished, 6),
    }


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
    accuracy_profile: str = "standard",
    **kwargs,
) -> None:
    metrics_json = kwargs.pop("metrics_json", None)
    threads_per_worker = kwargs.pop("threads_per_worker", 8)
    with PipelineMetrics(metrics_json) as metrics:
        return _execute(
            fasta_directory=fasta_directory,
            output_directory=output_directory,
            phmmer=phmmer,
            cpu=cpu,
            single_copy_threshold=single_copy_threshold,
            mcl=mcl,
            inflation_value=inflation_value,
            start=start,
            stop=stop,
            substitution_matrix=substitution_matrix,
            evalue_threshold=evalue_threshold,
            search_mode=search_mode,
            clustering=clustering,
            cpm_resolution=cpm_resolution,
            refinement_profile=refinement_profile,
            accuracy_profile=accuracy_profile,
            metrics=metrics,
            threads_per_worker=threads_per_worker,
            **kwargs,
        )


def _execute(
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
    search_mode: str,
    clustering: str,
    cpm_resolution,
    refinement_profile: str,
    accuracy_profile: str,
    metrics: PipelineMetrics,
    threads_per_worker: int,
    **kwargs,
) -> None:
    # for reporting runtime duration to user
    start_time = time.time()

    # make working dir
    working_dir = f"{output_directory}/orthohmm_working_res/"
    os.makedirs(working_dir, exist_ok=True)

    files = fetch_fasta_files(fasta_directory)
    phylogeny = kwargs.pop("phylogeny", "off")
    species_tree_mode = kwargs.pop("species_tree_mode", "supplied")
    species_tree = kwargs.pop("species_tree", None)
    aligner = kwargs.pop("aligner", "mafft")
    tree_builder = kwargs.pop("tree_builder", "FastTree")
    phylogeny_root_rule = kwargs.pop(
        "phylogeny_root_rule", "supported_children"
    )
    phylogeny_pair_rule = kwargs.pop("phylogeny_pair_rule", "lca")
    accuracy_config = resolve_accuracy_profile(accuracy_profile)
    if accuracy_config.multipass_graph and (
        search_mode != "builtin" or clustering != "leiden"
    ):
        raise ValueError(
            "high_sensitivity accuracy requires built-in search and Leiden clustering"
        )
    metrics.add_metadata(
        fasta_directory=os.path.abspath(fasta_directory),
        output_directory=os.path.abspath(output_directory),
        search_mode=search_mode,
        clustering=clustering,
        cpm_resolution=cpm_resolution,
        substitution_matrix=substitution_matrix.value,
        evalue_threshold=evalue_threshold,
        accuracy_profile=accuracy_config.name,
        search_kmer_k=accuracy_config.kmer_k,
        search_max_candidates_per_query=(
            accuracy_config.max_candidates_per_query
        ),
        leiden_seed=accuracy_config.leiden_seed,
        cpu_budget=cpu,
    )
    metrics.add_counts(species=len(files), species_pairs=len(files) ** 2)
    if phylogeny != "off":
        metrics.add_metadata(
            phylogeny=phylogeny,
            species_tree_mode=species_tree_mode,
            species_tree=species_tree,
            aligner=aligner,
            tree_builder=tree_builder,
            phylogeny_root_rule=phylogeny_root_rule,
            phylogeny_pair_rule=phylogeny_pair_rule,
        )

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
        accuracy_config.name,
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
        from .search.engine import resolve_parallelism
        search_workers, worker_threads = resolve_parallelism(
            cpu, len(files) ** 2, threads_per_worker
        )
        metrics.add_metadata(
            search_workers=search_workers,
            search_threads_per_worker=worker_threads,
            search_total_threads=search_workers * worker_threads,
        )
        with metrics.stage("search"):
            search_results = execute_builtin_search(
                files,
                fasta_directory,
                output_directory,
                cpu,
                substitution_matrix,
                kmer_k=accuracy_config.kmer_k,
                max_candidates_per_query=(
                    accuracy_config.max_candidates_per_query
                ),
                evalue_threshold=evalue_threshold,
                threads_per_worker=threads_per_worker,
            )
        metrics.add_counts(
            search_candidates=sum(
                result.candidate_count for result in search_results.values()
            ),
            significant_hits=sum(
                len(result.scores) for result in search_results.values()
            ),
        )
        print("\r          Completed!      \n")
        current_step += 1
    else:
        print(f"Step {current_step}/{total_steps}: Conducting all-to-all comparisons.")
        with metrics.stage("search"):
            execute_phmmer_search(
                phmmer_cmds,
                cpu,
            )
        print("\r          Completed!      \n")
        current_step += 1

    print(f"Step {current_step}/{total_steps}: Determining edge thresholds")
    with metrics.stage("edge_thresholds"):
        if accuracy_config.multipass_graph:
            gene_lengths = get_sequence_lengths(fasta_directory, files)
            gene_order = np.argsort(
                np.asarray([str(row["name"]) for row in gene_lengths]),
                kind="stable",
            )
            gene_lengths = gene_lengths[gene_order]
            reciprocal_best_hit_thresholds = None
            pairwise_rbh_corr = None
        else:
            gene_lengths, reciprocal_best_hit_thresholds, pairwise_rbh_corr = \
                determine_edge_thresholds(
                files,
                fasta_directory,
                output_directory,
                cpu,
                evalue_threshold,
                search_results=search_results,
            )
    metrics.add_counts(genes=len(gene_lengths))
    print("\r          Completed!      \n")
    current_step += 1

    print(f"Step {current_step}/{total_steps}: Identifying network edges")
    accuracy_hits = None
    with metrics.stage("network_edges"):
        if accuracy_config.multipass_graph:
            if search_results is None:
                raise ValueError(
                    "high_sensitivity accuracy requires in-memory built-in search results"
                )
            gene_names = [str(row["name"]) for row in gene_lengths]
            gene_to_id = {gene: idx for idx, gene in enumerate(gene_names)}
            species_to_id = {}
            gene_to_species = []
            for row in gene_lengths:
                species = str(row["spp"])
                species_to_id.setdefault(species, len(species_to_id))
                gene_to_species.append(species_to_id[species])
            accuracy_hits = _collect_search_hit_arrays(
                search_results,
                evalue_threshold,
                gene_to_id,
            )
            checkpoint_path = write_accuracy_checkpoint(
                output_directory,
                gene_names,
                gene_to_species,
                *accuracy_hits,
            )
            metrics.add_metadata(
                high_sensitivity_checkpoint=str(checkpoint_path.resolve())
            )
            edges = build_rbnh_edges(
                gene_names,
                gene_to_species,
                *accuracy_hits,
            )
            rbnh_edges = edges
            write_indexed_edges(edges, output_directory)
        else:
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
    metrics.add_counts(network_edges=len(edges))
    print("\r          Completed!      \n")
    current_step += 1

    print(f"Step {current_step}/{total_steps}: Conducting clustering")
    with metrics.stage("clustering"):
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
                edges=edges,
                include_isolates=accuracy_config.multipass_graph,
                seed=accuracy_config.leiden_seed,
            )
            if accuracy_config.multipass_graph:
                cluster_path = (
                    f"{output_directory}/orthohmm_working_res/"
                    "orthohmm_edges_clustered.txt"
                )
                initial_clusters = read_index_clusters(
                    cluster_path,
                    edges.gene_names,
                )
                singleton_edges = build_singleton_assignment_edges(
                    edges.gene_names,
                    initial_clusters,
                    *accuracy_hits,
                )
                metrics.add_counts(
                    high_sensitivity_rbnh_edges=len(edges),
                    high_sensitivity_singleton_edges=len(singleton_edges),
                )
                edges = combine_edges(edges, singleton_edges)
                write_indexed_edges(edges, output_directory)
                execute_leiden(
                    cpm_resolution,
                    output_directory,
                    edges=edges,
                    include_isolates=True,
                    seed=accuracy_config.leiden_seed,
                )
                metrics.add_counts(network_edges=len(edges))
    if accuracy_config.profile_expansion:
        from .search.profile_expansion import expand_profiles

        with metrics.stage("profile_expansion"):
            seed_clusters = read_index_clusters(
                cluster_path,
                edges.gene_names,
            )
            profile_result = expand_profiles(
                seed_clusters,
                edges.gene_names,
                fasta_directory,
                files,
                substitution_matrix.value,
                cpu,
                evalue_threshold,
                *accuracy_hits,
            )
            profile_base_edges = combine_edges(
                rbnh_edges,
                profile_result.edges,
            )
            execute_leiden(
                cpm_resolution,
                output_directory,
                edges=profile_base_edges,
                include_isolates=True,
                seed=accuracy_config.leiden_seed,
            )
            profile_base_clusters = read_index_clusters(
                cluster_path,
                edges.gene_names,
            )
            final_singleton_edges = build_singleton_assignment_edges(
                edges.gene_names,
                profile_base_clusters,
                *accuracy_hits,
            )
            edges = combine_edges(
                profile_base_edges,
                final_singleton_edges,
            )
            write_indexed_edges(edges, output_directory)
            execute_leiden(
                cpm_resolution,
                output_directory,
                edges=edges,
                include_isolates=True,
                seed=accuracy_config.leiden_seed,
            )
        metrics.add_counts(
            high_sensitivity_profiles=profile_result.profiles_built,
            high_sensitivity_profile_candidates=(
                profile_result.profile_candidates
            ),
            high_sensitivity_profile_hits=(
                profile_result.significant_profile_hits
            ),
            high_sensitivity_profile_edges=len(profile_result.edges),
            high_sensitivity_final_singleton_edges=(
                len(final_singleton_edges)
            ),
            network_edges=len(edges),
        )
    with metrics.stage("refinement"):
        refinement_substages = _refine_cluster_file(
            output_directory,
            gene_lengths,
            search_results,
            edges,
            evalue_threshold,
            refinement_profile=refinement_profile,
        )
    metrics.add_metadata(refinement_substages_s=refinement_substages)
    if phylogeny == "reconcile":
        from .phylogeny import PhylogenyConfig
        from .phylogeny_pipeline import run_phylogeny_stage

        print("          Reconciling ambiguous families with the species tree")
        with metrics.stage("phylogeny"):
            phylogeny_result = run_phylogeny_stage(
                fasta_directory=fasta_directory,
                output_directory=output_directory,
                files=files,
                config=PhylogenyConfig(
                    mode=phylogeny,
                    species_tree_mode=species_tree_mode,
                    species_tree=species_tree,
                    aligner=aligner,
                    tree_builder=tree_builder,
                    root_duplication_rule=phylogeny_root_rule,
                    pair_orthology_rule=phylogeny_pair_rule,
                ),
                cpu=cpu,
            )
        metrics.add_counts(
            phylogeny_candidate_families=phylogeny_result.candidate_families,
            phylogeny_reconciled_families=phylogeny_result.reconciled_families,
            phylogeny_bypassed_families=phylogeny_result.bypassed_families,
            phylogeny_checkpoint_hits=phylogeny_result.checkpoint_hits,
            phylogeny_root_hogs=phylogeny_result.root_hogs,
            phylogeny_ortholog_pairs=phylogeny_result.ortholog_pairs,
            phylogeny_duplications=phylogeny_result.duplications,
            phylogeny_speciations=phylogeny_result.speciations,
            phylogeny_uncertain_events=phylogeny_result.uncertain_events,
            phylogeny_species_tree_families=(
                phylogeny_result.species_tree_families
            ),
            phylogeny_species_tree_checkpoint_hit=(
                phylogeny_result.species_tree_checkpoint_hit
            ),
        )
    with metrics.stage("orthogroup_materialization"):
        singletons, og_cn, ogs_dat, single_copy_ogs = \
            generate_orthogroup_clusters_file(
                output_directory,
                gene_lengths,
                files,
                single_copy_threshold,
                fasta_directory,
            )
    metrics.add_counts(
        orthogroups=len(ogs_dat),
        singletons=len(singletons),
        single_copy_orthogroups=len(single_copy_ogs),
    )
    print("          Completed!\n")
    current_step += 1

    # exit if users only want orthogroups to be inferred
    if stop == StopStep.infer:
        sys.exit()

    print(f"Step {current_step}/{total_steps}: Writing orthogroup information")
    with metrics.stage("output"):
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
