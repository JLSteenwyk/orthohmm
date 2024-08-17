import itertools
import multiprocessing
import subprocess
from typing import List


def run_bash_commands(
    command: List[str]
) -> None:
    subprocess.run(
        command,
        shell=True,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def run_phmmer(
    phmmer_cmds: List[str],
    cpu: int,
) -> None:
    # Create a pool of workers
    pool = multiprocessing.Pool(processes=cpu)

    # # Map the commands to the worker pool
    # pool.map(run_bash_commands, phmmer_cmds)
    for command in phmmer_cmds:
        pool.apply_async(run_bash_commands, args=(command,))

    # Close the pool and wait for the work to finish
    pool.close()
    pool.join()


def execute_phmmer_search(
        files: List[str],
        cpu: int,
        fasta_directory: str,
        phmmer: str,
        temporary_directory: str,
) -> None:
    # create phmmer cmds
    pairwise_combos = list(itertools.product(files, repeat=2))
    phmmer_cmds = []
    for combo in pairwise_combos:
        # TODO: ask Thomas is the setting CPUs here and in the multi threader is okay
        # phmmer_cmds.append(f"{phmmer} --noali --notextw --cpu {cpu} --tblout {temporary_directory}/{combo[0]}_2_{combo[1]}.phmmerout.txt {fasta_directory}/{combo[0]} {fasta_directory}/{combo[1]}")
        phmmer_cmds.append(f"{phmmer} --noali --notextw --tblout {temporary_directory}/{combo[0]}_2_{combo[1]}.phmmerout.txt {fasta_directory}/{combo[0]} {fasta_directory}/{combo[1]}")

    # run phmmer cmds
    run_phmmer(phmmer_cmds, cpu)


def execute_mcl(
    mcl: str,
    inflation_value: float,
    cpu: int,
    temporary_directory: str,
) -> None:
    cmd = f"{mcl} {temporary_directory}/orthohmm_edges.txt -te {cpu} --abc -I {inflation_value} -o {temporary_directory}/orthohmm_edges_clustered.txt"
    subprocess.run(
        cmd,
        shell=True,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
