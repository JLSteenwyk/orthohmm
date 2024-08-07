#!/usr/bin/env python

import logging
import os
import shutil
import sys
import time
import textwrap

from .args_processing import process_args
from .helpers import (
    determine_edge_thresholds,
    determine_network_edges,
    execute_mcl,
    execute_phmmer_search,
    generate_orthogroup_files,
)
from .parser import create_parser
from .version import __version__
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
    **kwargs,
):
    print(textwrap.dedent(
            f"""\
          ____       _   _           _    _ __  __ __  __ 
         / __ \     | | | |         | |  | |  \/  |  \/  |
        | |  | |_ __| |_| |__   ___ | |__| | \  / | \  / |
        | |  | | '__| __| '_ \ / _ \|  __  | |\/| | |\/| |
        | |__| | |  | |_| | | | (_) | |  | | |  | | |  | |
         \____/|_|   \__|_| |_|\___/|_|  |_|_|  |_|_|  |_|

        Version: {__version__}
        """  # noqa
    ))

    # for reporting runtime duration to user
    start_time = time.time()

    # column names and columns to drop
    col_names = [
        "target name", "accession", "query name", "accession", "E-value", "score", "bias", "domain E-value", "domain score",
        "domain bias", "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc", "description of target"
    ]

    columns_to_drop = [
        "accession", "E-value", "bias", "domain E-value", "domain score", "domain bias", "exp",
        "reg", "clu", "ov", "env", "dom", "rep", "inc", "description of target"
    ]

    # create working directory
    if not os.path.isdir(f"{output_directory}/working_dir"):
        os.mkdir(f"{output_directory}/working_dir")

    # get FASTA files to identify orthologs from
    extensions = (".fa", ".faa", ".fas", ".fasta")
    files = [file for file in os.listdir(fasta_directory) if os.path.splitext(file)[1].lower() in extensions]

    # display to user what args are being used in stdout
    write_user_args(
        fasta_directory,
        output_directory,
        phmmer,
        mcl,
        cpu,
        single_copy_threshold,
        files
    )

    # Step 1: all-to-all comparisons
    print("Step 1/6: Conducting all-to-all comparisons.")
    print("          This is typically the longest step.")
    execute_phmmer_search(
        files,
        cpu,
        output_directory,
        fasta_directory,
        phmmer,
    )
    print("          Completed!\n")

    # Step 2: Determining edge thresholds
    print("Step 2/6: Determining edge thresholds")
    gene_lengths, reciprocal_best_hit_thresholds, pairwise_rbh_corr = \
        determine_edge_thresholds(
            files,
            fasta_directory,
            output_directory,
            col_names,
            columns_to_drop,
        )
    print("          Completed!\n")

    # Step 3: Determining network edges
    print("Step 3/6: Identifying network edges")
    edges = determine_network_edges(
        files,
        output_directory,
        col_names,
        columns_to_drop,
        gene_lengths,
        pairwise_rbh_corr,
        reciprocal_best_hit_thresholds,
    )
    print("          Completed!\n")

    # Step 4: Conduct mcl clustering
    print("Step 4/6: Conducting clustering")
    execute_mcl(mcl, inflation_value, output_directory)
    print("          Completed!\n")

    # Step 5: Write out orthogroup files
    print("Step 5/6: Writing orthogroup information")
    single_copy_ogs, singletons, ogs_dat = generate_orthogroup_files(
        output_directory,
        gene_lengths,
        files,
        fasta_directory,
        single_copy_threshold,
    )
    print("          Completed!\n")

    # Step 6: Clean
    print("Step 6/6: Cleaning up workspace")
    shutil.rmtree(f"{output_directory}/working_dir")
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
