#!/usr/bin/env python

import logging
import os
import sys
import time
from typing import Union

from .args_processing import process_args
from .externals import execute_mcl
from .files import fetch_fasta_files
from .helpers import (
    determine_edge_thresholds_pyhmmer,
    determine_network_edges_pyhmmer,
    generate_orthogroup_clusters_file,
    generate_orthogroup_files,
    StartStep,
    StopStep,
    SubstitutionMatrix,
)
from .parser import create_parser
from .pyhmmer_externals import execute_pyhmmer_search_all_pairs
from .writer import (
    write_user_args,
    write_output_stats
)

logger = logging.getLogger(__name__)


def execute_pyhmmer(
    fasta_directory: str,
    output_directory: str,
    phmmer: str,  # Not used in pyhmmer version but kept for compatibility
    cpu: int,
    single_copy_threshold: float,
    mcl: str,
    inflation_value: float,
    start: Union[StartStep, None],
    stop: Union[StopStep, None],
    substitution_matrix: SubstitutionMatrix,
    evalue_threshold: float,
    **kwargs,
) -> None:
    # for reporting runtime duration to user
    start_time = time.time()

    # make working dir
    working_dir = f"{output_directory}/orthohmm_working_res/"
    os.makedirs(working_dir, exist_ok=True)

    files = fetch_fasta_files(fasta_directory)

    # display to user what args are being used in stdout
    write_user_args(
        fasta_directory,
        output_directory,
        "pyhmmer",  # Indicate we're using pyhmmer
        mcl,
        cpu,
        single_copy_threshold,
        files,
        start,
        stop,
        substitution_matrix,
        evalue_threshold,
        inflation_value,
    )

    # set current step and determine the total number of
    # steps that will be used in the run
    current_step = 1
    if stop == StopStep.infer:
        total_steps = 4
    else:
        total_steps = 5

    pyhmmer_results = None
    
    if start == StartStep.search_res:
        total_steps -= 1
        # Would need to load previous results from files
        print("Starting from search results not yet implemented for pyhmmer version")
        sys.exit(1)
    else:
        print(f"Step {current_step}/{total_steps}: Conducting all-to-all comparisons with PyHMMER.")
        pyhmmer_results = execute_pyhmmer_search_all_pairs(
            files,
            fasta_directory,
            substitution_matrix,
            evalue_threshold,
            cpu,
        )
        print("\r          Completed!      \n")
        current_step += 1

    print(f"Step {current_step}/{total_steps}: Determining edge thresholds")
    gene_lengths, reciprocal_best_hit_thresholds, pairwise_rbh_corr = \
        determine_edge_thresholds_pyhmmer(
            files,
            fasta_directory,
            substitution_matrix,
            cpu,
            evalue_threshold,
            pyhmmer_results,
        )
    print("\r          Completed!      \n")
    current_step += 1

    print(f"Step {current_step}/{total_steps}: Identifying network edges")
    edges = determine_network_edges_pyhmmer(
        files,
        gene_lengths,
        pairwise_rbh_corr,
        reciprocal_best_hit_thresholds,
        evalue_threshold,
        pyhmmer_results,
        cpu,
        output_directory,
    )
    print("\r          Completed!      \n")
    current_step += 1

    print(f"Step {current_step}/{total_steps}: Conducting clustering")
    execute_mcl(
        mcl,
        inflation_value,
        cpu,
        output_directory,
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


def main_pyhmmer(argv=None):
    """
    PyHMMER version of main function
    """
    parser = create_parser()
    args = parser.parse_args()

    execute_pyhmmer(**process_args(args))


if __name__ == "__main__":
    main_pyhmmer(sys.argv[1:])