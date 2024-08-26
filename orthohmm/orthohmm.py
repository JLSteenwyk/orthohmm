#!/usr/bin/env python

import logging
import os
import sys
import time
from typing import Union, List

from .args_processing import process_args
from .externals import (
    execute_mcl,
    execute_phmmer_search,
)
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
from .writer import (
    write_user_args,
    write_output_stats
)

logger = logging.getLogger(__name__)


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
    **kwargs,
) -> None:
    # for reporting runtime duration to user
    start_time = time.time()

    # make working dir
    if not os.path.exists(f"{output_directory}/orthohmm_working_res/"):
        os.makedirs(f"{output_directory}/orthohmm_working_res/")

    # get FASTA files to identify orthologs from
    extensions = (".fa", ".faa", ".fas", ".fasta")
    files: List[str] = [file for file in os.listdir(fasta_directory) if os.path.splitext(file)[1].lower() in extensions]

    if start != StartStep.search_res:
        phmmer_cmds = generate_phmmer_cmds(
            files,
            phmmer,
            output_directory,
            fasta_directory,
            cpu,
            stop,
            substitution_matrix,
        )

    # print phmmer cmds and exit is users only want to prepare phmmer cmds
    if stop == StopStep.prepare:
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
    else:
        print(f"Step {current_step}/{total_steps}: Conducting all-to-all comparisons.")
        print("          This is typically the longest step.")
        execute_phmmer_search(
            phmmer_cmds,
            cpu,
        )
        print("\r          Completed!   \n")
        current_step += 1

    print(f"Step {current_step}/{total_steps}: Determining edge thresholds")
    gene_lengths, reciprocal_best_hit_thresholds, pairwise_rbh_corr = \
        determine_edge_thresholds(
            files,
            fasta_directory,
            output_directory,
        )
    print("          Completed!\n")
    current_step += 1

    print(f"Step {current_step}/{total_steps}: Identifying network edges")
    edges = determine_network_edges(
        files,
        output_directory,
        gene_lengths,
        pairwise_rbh_corr,
        reciprocal_best_hit_thresholds,
    )
    print("          Completed!\n")
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
        extensions,
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
