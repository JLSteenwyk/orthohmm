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


def execute_phmmer_search(
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


def execute_mcl(
    mcl: str,
    inflation_value: float,
    cpu: int,
    output_directory: str,
) -> None:
    cmd = f"{mcl} {output_directory}/orthohmm_working_res/orthohmm_edges.txt -te {cpu} --abc -I {inflation_value} -o {output_directory}/orthohmm_working_res/orthohmm_edges_clustered.txt"
    subprocess.run(
        cmd,
        shell=True,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
