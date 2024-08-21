import textwrap
import time
from typing import List, Union

from .helpers import StartStep, StopStep
from .version import __version__


def write_user_args(
    fasta_directory: str,
    output_directory: str,
    phmmer: str,
    mcl: str,
    cpu: int,
    single_copy_threshold: float,
    files: List[str],
    start: Union[StartStep, None],
    stop: Union[StopStep, None],
) -> None:

    try:
        if start.value:
            start_print = start.value
    except AttributeError:
        start_print = "NA"

    try:
        if stop.value:
            stop_print = stop.value
    except AttributeError:
        stop_print = "NA"

    print(
        textwrap.dedent(
            f"""\
          ____       _   _           _    _ __  __ __  __ 
         / __ \     | | | |         | |  | |  \/  |  \/  |
        | |  | |_ __| |_| |__   ___ | |__| | \  / | \  / |
        | |  | | '__| __| '_ \ / _ \|  __  | |\/| | |\/| |
        | |__| | |  | |_| | | | (_) | |  | | |  | | |  | |
         \____/|_|   \__|_| |_|\___/|_|  |_|_|  |_|_|  |_|

        Version: {__version__}

    -------------
    | Arguments |
    -------------
    Directory of FASTA files: {fasta_directory}
    Number of FASTA files: {len(files)}
    Directory for output files: {output_directory}
    Step to start analysis: {start_print}
    Step to stop analysis: {stop_print}
    Path to phmmer: {phmmer}
    Path to mcl: {mcl}
    Single-copy threshold: {single_copy_threshold}
    CPUs: {cpu}

    """  # noqa
        )
    )


def write_output_stats(
    start_time: float,
    single_copy_ogs: list,
    singletons: list,
    ogs_dat: dict,
    edges: dict,
    gene_lengths: list,
) -> None:
    """
    Function to print out output statistics
    """
    print(
        textwrap.dedent(
            f"""\

        ---------------------
        | Output Statistics |
        ---------------------
        Number of genes processed: {len(gene_lengths)}
        Number of orthogroups: {len(ogs_dat)}
        Number of edges in network: {len(edges)}
        Number of single-copy orthogroups: {len(single_copy_ogs)}
        Number of singletons: {len(singletons)}

        Execution time: {round(time.time() - start_time, 3)}s
    """  # noqa
        )
    )
