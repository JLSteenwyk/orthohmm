import os
import subprocess
import sys
from typing import List
import multiprocessing
from multiprocessing.synchronize import Lock
from multiprocessing.sharedctypes import Synchronized


def run_bash_command(command: str) -> None:
    subprocess.run(
        command.split(),
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
        sys.stdout.write(f"\r          {progress:.1f}% complete")
        sys.stdout.flush()


def execute_phmmer_search(
    phmmer_cmds: List[str],
    cpu: int,
) -> None:

    # create a pool of workers
    pool = multiprocessing.Pool(processes=cpu)

    # create a counter and lock for tracking progress
    completed_tasks = multiprocessing.Value('i', 0)
    total_tasks = len(phmmer_cmds)
    lock = multiprocessing.Lock()

    # apply async with a callback to update progress
    for command in phmmer_cmds:
        pool.apply_async(
            run_bash_command,
            args=(command,),
            callback=lambda _: update_progress(
                lock, completed_tasks, total_tasks
            )
        )

    # close the pool and wait for the work to finish
    pool.close()
    pool.join()


def check_if_phmmer_command_completed(
    file_to_check: str
) -> bool:
    if not os.path.isfile(file_to_check):
        return False

    with open(file_to_check, "r") as file:
        lines = file.readlines()
        if lines and lines[-1].strip() == "# [ok]":
            return True
    return False


def execute_mcl(
    mcl: str,
    inflation_value: float,
    cpu: int,
    output_directory: str,
) -> None:
    if not check_if_mcl_command_completed(f"{output_directory}/orthohmm_working_res/orthohmm_edges_clustered.txt"):
        cmd = f"{mcl} {output_directory}/orthohmm_working_res/orthohmm_edges.txt -te {cpu} --abc -I {inflation_value} -o {output_directory}/orthohmm_working_res/orthohmm_edges_clustered.txt"
        subprocess.run(
            cmd,
            shell=True,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )


def check_if_mcl_command_completed(
    file_to_check: str
) -> bool:
    if not os.path.isfile(file_to_check):
        return False

    with open(file_to_check, "r") as file:
        lines = file.readlines()
        if lines and lines[-1].strip() == "    ( http://link.aip.org/link/?SJMAEL/30/121/1 )":
            return True
    return False
