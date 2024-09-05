from distutils.spawn import find_executable
import shutil
import logging
import multiprocessing
import os.path
import sys

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

    if args.cpu:
        cpu = int(args.cpu)
        if cpu > multiprocessing.cpu_count():
            logger.warning(f"{cpu} CPUs requested exceeds {multiprocessing.cpu_count()} CPUs available.")
            logger.warning(f"Changing CPUs to {multiprocessing.cpu_count()}.")
            cpu = multiprocessing.cpu_count()
    else:
        cpu = multiprocessing.cpu_count()

    single_copy_threshold = float(args.single_copy_threshold) if args.single_copy_threshold is not None else 0.5

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
    )
