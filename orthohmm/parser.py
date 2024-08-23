import sys
import textwrap

from argparse import (
    ArgumentParser,
    SUPPRESS,
    RawDescriptionHelpFormatter,
)

from .helpers import (
    StartStep,
    StopStep,
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

        -x, --substitution_matrix <subs. matrix>    substitution matrix to use for
                                                    residue probabilities
                                                    (default: BLOSUM62)

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
    
        --stop <prepare, infer, write>              options for stopping
                                                    an analysis at a specific
                                                    intermediate step

        --start <search_res>                        start analysis from
                                                    completed all-vs-all
                                                    search results

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

        Substitution matrix (-x, --substitution_matrix)
            Residue alignment probabilities will be determined from the
            specified substitution matrix. Supported substitution matrices
            include: BLOSUM45, BLOSUM50, BLOSUM62, BLOSUM80, BLOSUM90,
            PAM30, PAM70, PAM120, and PAM240.

        CPU (-c, --cpu) 
            Number of CPU workers for multithreading during sequence search.
            This argument is used by phmmer during all-by-all comparisons.
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
        
        Stop (--stop)
            Similar to other ortholog calling algorithms, different steps in the
            OrthoHMM workflow can be cpu or memory intensive. Thus, users may
            want to stop OrthoHMM at certain steps, to faciltiate more
            practical resource allocation. There are three choices for when to
            stop the analysis: prepare, infer, and write.
            - prepare: Stop after preparing input files for the all-by-all search
            - infer: Stop after inferring the orthogroups
            - write: Stop after writing sequence files for the orthogroups.
                    Currently, this is synonymous with not specifying a step
                    to stop the analysis at.

        Start (--start)
            Start analysis from a specific intermediate step. Currently, this
            can only be applied to the results from the all-by-all search.
            - search_res: Start analysis from all-by-all search results.
                          
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
        
        orthohmm_working_res
            Various intermediate results files that help OrthoHMM start analyses
            from an intermediate step in the analysis. This includes outputs
            from phmmer searches, initial edges inputted to MCL, and the 
            output from MCL clustering.
        """  # noqa
        ),
    )

    optional.add_argument(
        "-o",
        "--output_directory",
        help=SUPPRESS,
        metavar="output_directory"
    )

    substitution_matrix_choices = [step.value for step in StartStep]
    optional.add_argument(
        "-x",
        "--substitution_matrix",
        type=str,
        required=False,
        help=SUPPRESS,
        metavar="substitution_model",
        choices=substitution_matrix_choices,
    )

    optional.add_argument(
        "-s",
        "--single_copy_threshold",
        type=float,
        required=False,
        help=SUPPRESS,
        metavar="single_copy_threshold",
    )

    stop_choices = [step.value for step in StopStep]
    optional.add_argument(
        "--stop",
        type=str,
        required=False,
        help=SUPPRESS,
        metavar="stop step",
        choices=stop_choices
    )

    start_choices = [step.value for step in StartStep]
    optional.add_argument(
        "--start",
        type=str,
        required=False,
        help=SUPPRESS,
        metavar="start step",
        choices=start_choices
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
