import math
import multiprocessing
from multiprocessing.synchronize import Lock
from multiprocessing.sharedctypes import Synchronized
import subprocess
import sys
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


def update_progress(
    lock: Lock,
    completed_tasks: Synchronized,
    total_tasks: int,
) -> None:
    with lock:
        completed_tasks.value += 1
        progress = (completed_tasks.value / total_tasks) * 100
        sys.stdout.write(f"\r          {math.floor(progress)}% complete")
        sys.stdout.flush()


def execute_phmmer_search(
    phmmer_cmds: List[str],
    cpu: int,
) -> None:

    # Create a pool of workers
    pool = multiprocessing.Pool(processes=cpu)

    # Create a counter and lock for tracking progress
    completed_tasks = multiprocessing.Value("i", 0)
    total_tasks = len(phmmer_cmds)
    lock = multiprocessing.Lock()

    # Apply async with a callback to update progress
    for command in phmmer_cmds:
        pool.apply_async(
            run_bash_commands,
            args=(command,),
            callback=lambda _: update_progress(
                lock, completed_tasks, total_tasks
            )
        )

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
