import shutil
import logging
import multiprocessing
import os.path
import sys

# distutils.spawn was removed in Python 3.12; shutil.which is the
# standard-library replacement that's been available since Python 3.3.
find_executable = shutil.which

from .helpers import (
    StartStep,
    StopStep,
    SubstitutionMatrix,
)


logger = logging.getLogger(__name__)


def process_args(args) -> dict:
    """
    Process args from argparser and set defaults
    """
    # required argument
    fasta_directory = args.fasta_directory

    if not os.path.isdir(fasta_directory):
        logger.warning("Input directory does not exist")
        sys.exit()

    # assign optional arguments
    output_directory = args.output_directory or args.fasta_directory
    if not os.path.isdir(output_directory):
        logger.warning("Output directory does not exist")
        sys.exit()

    search_mode = getattr(args, "search_mode", None) or "builtin"

    # Only require phmmer when using phmmer search mode
    phmmer = None
    if search_mode == "phmmer":
        if args.phmmer:
            phmmer = args.phmmer
            if not shutil.which(phmmer):
                logger.warning(f"phmmer can't be found at {phmmer}.")
                sys.exit()
        else:
            if find_executable("phmmer"):
                phmmer = "phmmer"
            else:
                logger.warning("phmmer can't be found. Provide path with the -p argument or add path to PATH variable")
                sys.exit()
    else:
        # Built-in mode: phmmer path is optional
        if args.phmmer:
            phmmer = args.phmmer
        elif find_executable("phmmer"):
            phmmer = "phmmer"
        else:
            phmmer = "builtin"

    if args.cpu:
        cpu = int(args.cpu)
        if cpu > multiprocessing.cpu_count():
            logger.warning(f"{cpu} CPUs requested exceeds {multiprocessing.cpu_count()} CPUs available.")
            logger.warning(f"Changing CPUs to {multiprocessing.cpu_count()}.")
            cpu = multiprocessing.cpu_count()
    else:
        cpu = multiprocessing.cpu_count()

    single_copy_threshold = float(args.single_copy_threshold) if args.single_copy_threshold is not None else 0.5

    clustering = getattr(args, "clustering", None) or "leiden"
    cpm_resolution = getattr(args, "cpm_resolution", None)
    if cpm_resolution is None:
        cpm_resolution = "0.1"
    # Pass-through string. The library accepts a float or "auto"; we
    # validate/coerce here so downstream sees a clean value.
    if isinstance(cpm_resolution, str) and cpm_resolution.lower() != "auto":
        try:
            cpm_resolution = float(cpm_resolution)
        except ValueError:
            logger.warning(f"--cpm_resolution must be a float or 'auto', got {cpm_resolution!r}")
            sys.exit()

    # MCL is only required when the user explicitly opts into it. The default
    # `leiden` clustering uses Python igraph + leidenalg (pip-installed).
    if clustering == "mcl":
        if args.mcl:
            mcl = args.mcl
            if not shutil.which(mcl):
                logger.warning(f"mcl can't be found at {mcl}.")
                sys.exit()
        else:
            if find_executable("mcl"):
                mcl = "mcl"
            else:
                logger.warning("mcl can't be found. Provide path with the -m argument or add path to PATH variable")
                sys.exit()
    else:
        # Leiden mode: mcl is optional. Resolve it if available so legacy
        # callers passing --mcl still get a sensible value back.
        if args.mcl:
            mcl = args.mcl
        elif find_executable("mcl"):
            mcl = "mcl"
        else:
            mcl = "mcl"  # placeholder; never invoked in leiden mode

    inflation_value = args.inflation_value or 1.5

    start = StartStep(args.start) if args.start else None
    stop = StopStep(args.stop) if args.stop else None

    substitution_matrix = SubstitutionMatrix(args.substitution_matrix) if args.substitution_matrix else SubstitutionMatrix.blosum62

    evalue_threshold = args.evalue or 0.0001

    return dict(
        fasta_directory=fasta_directory,
        output_directory=output_directory,
        phmmer=phmmer,
        cpu=int(cpu),
        single_copy_threshold=single_copy_threshold,
        mcl=mcl,
        inflation_value=float(inflation_value),
        start=start,
        stop=stop,
        substitution_matrix=substitution_matrix,
        evalue_threshold=evalue_threshold,
        search_mode=search_mode,
        clustering=clustering,
        cpm_resolution=cpm_resolution,
    )
