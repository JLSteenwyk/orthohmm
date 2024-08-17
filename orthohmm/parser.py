import sys
import textwrap

from argparse import (
    ArgumentParser,
    SUPPRESS,
    RawDescriptionHelpFormatter,
)

from .version import __version__


def create_parser() -> ArgumentParser:
    parser = ArgumentParser(
        add_help=False,
        formatter_class=RawDescriptionHelpFormatter,
        usage=SUPPRESS,
        description=textwrap.dedent(
            f"""\
          ____       _   _           _    _ __  __ __  __ 
         / __ \     | | | |         | |  | |  \/  |  \/  |
        | |  | |_ __| |_| |__   ___ | |__| | \  / | \  / |
        | |  | | '__| __| '_ \ / _ \|  __  | |\/| | |\/| |
        | |__| | |  | |_| | | | (_) | |  | | |  | | |  | |
         \____/|_|   \__|_| |_|\___/|_|  |_|_|  |_|_|  |_|

        
        Version: {__version__}
        Citation: Steenwyk et al. YEAR, JOURNAL. doi: DOI
        LINK

        HMM-based inference of orthologous groups.

        Usage: orthohmm <input> [optional arguments]
        """  # noqa
        ),
    )

    # if no arguments are given, print help and exit
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit()

    # required arguments
    required = parser.add_argument_group(
        "required argument",
        description=textwrap.dedent(
            """\
        <fasta_directory>                           Directory of FASTA files ending in
                                                    .fa, .faa, .fas, or .fasta
                                                    (must be the first argument)
        """
        ),
    )

    required.add_argument("fasta_directory", type=str, help=SUPPRESS)

    # optional arguments
    optional = parser.add_argument_group(
        "optional arguments",
        description=textwrap.dedent(
            """\
        -o, --output_directory <path>               output directory name 
                                                    (default: same directory as
                                                    directory of FASTA files)

        -p, --phmmer <path>                         path to phmmer from HMMER suite
                                                    (default: phmmer)

        -c, --cpu <integer>                         number of parallel CPU workers
                                                    to use for multithreading
                                                    (default: auto detect)
        
        -s, --single_copy_threshold <float>         taxon occupancy threshold
                                                    for single-copy orthologs
                                                    (default: 0.5)
        
        -m, --mcl <path>                            path to mcl software
                                                    (default: mcl)

        -i, --inflation_value <float>               mcl inflation parameter
                                                    (default: 1.5)

        -t, --temporary_directory <path>            specify path for temporary
                                                    directory of intermediate
                                                    files
                                                    (default: /tmp/)

        -------------------------------------
        | Detailed explanation of arguments | 
        -------------------------------------
        Output directory (-o, --output_directory)
            Output directory name to store OrthoHMM results. This directory
            should already exist. By default, results files will be written
            to the same directory as the input directory of FASTA files.

        Phmmer (-p, --phmmer) 
            Path to phmmer executable from HMMER suite. By default, phmmer
            is assumed to be in the PATH variable; in other words, phmmer
            can be evoked by typing `phmmer`.

        CPU (-c, --cpu) 
            Number of CPU workers for multithreading during sequence search.
            This argument is used by phmmer during all-vs-all comparisons.
            By default, the number of CPUs available will be auto-detected.
        
        Single-Copy Threshold (-s, --single_copy_threshold)
            Taxon occupancy threshold when identifying single-copy orthologs.
            By default, the threshold is 50% taxon occupancy, which is specified
            as a fraction - that is, 0.5.
        
        MCL (-m, --mcl)
            Path to mcl executable from MCL software. By default, mcl
            is assumed to be in the PATH variable; in other words,
            mcl can be evoked by typing `mcl`.

        Inflation Value (-i, --inflation_value)
            MCL inflation parameter for clustering genes into orthologous groups.
            Lower values are more permissive resulting in larger orthogroups.
            Higher values are stricter resulting in smaller orthogroups.
            The default value is 1.5.
        
        Temporary Directory (-t, --temporary_directory)
            Path for temporary directory location. This temporary directory will
            have the prefix "orthohmm" and a randomly generated suffix. Files
            like the output of phmmer are kept here while OrthoHMM processes them.
            The default path is /tmp/.

        -------------------
        | OrthoHMM output | 
        -------------------
        All OrthoHMM outputs have the prefix `orthohmm` so that they are easy to find.

        orthohmm_gene_count.txt
            A gene count matrix per taxa for each orthogroup. Space delimited.
        
        orthohmm_orthogroups.txt
            Genes present in each orthogroup. Space delimited.
        
        orthohmm_single_copy_orthogroups.txt
            A single-column list of single-copy orthologs.
        
        orthohmm_orthogroups
            A directory of FASTA files wherein each file is an orthogroup.
        
        orthohmm_single_copy_orthogroups
            A directory of FASTA files wherein each file is a single-copy ortholog.
            Headers are modified to have taxon names come before the gene identifier.
            Taxon names are the file name excluding the extension. Taxon name and gene
            identifier are separated by a pipe symbol "|". This aims to help streamline
            phylogenomic workflows wherein sequences will be concatenated downstream
            based on taxon names.
        """  # noqa
        ),
    )

    optional.add_argument(
        "-o",
        "--output_directory",
        help=SUPPRESS,
        metavar="output_directory"
    )

    optional.add_argument(
        "-t",
        "--temporary_directory",
        help=SUPPRESS,
        metavar="temporary_directory"
    )

    optional.add_argument(
        "-s",
        "--single_copy_threshold",
        type=float,
        required=False,
        help=SUPPRESS,
        metavar="single_copy_threshold",
    )

    optional.add_argument(
        "-c",
        "--cpu",
        help=SUPPRESS,
        metavar="cpu"
    )

    optional.add_argument(
        "-p",
        "--phmmer",
        help=SUPPRESS,
        metavar="phmmer"
    )

    optional.add_argument(
        "-m",
        "--mcl",
        help=SUPPRESS,
        metavar="mcl"
    )

    optional.add_argument(
        "-i",
        "--inflation_value",
        type=float,
        required=False,
        help=SUPPRESS,
        metavar="inflation_value",
    )

    optional.add_argument(
        "-h",
        "--help",
        action="help",
        help=SUPPRESS,
    )

    optional.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"orthohmm v{__version__}",
        help=SUPPRESS,
    )

    return parser
