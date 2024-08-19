import textwrap
import time


def write_user_args(
    fasta_directory: str,
    output_directory: str,
    phmmer: str,
    mcl: str,
    cpu: int,
    single_copy_threshold: float,
    files: list,
    temporary_directory: str,
) -> None:

    print(
        textwrap.dedent(
            f"""\

    -------------
    | Arguments |
    -------------
    Directory of FASTA files: {fasta_directory}
    Number of FASTA files: {len(files)}
    Directory for output files: {output_directory}
    Temporary directory: {temporary_directory}
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
