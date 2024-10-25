import math
import os
import subprocess
import sys
from typing import List


def run_bash_command(command: str) -> None:
    subprocess.run(
        command.split(),
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def update_progress(
    completed_tasks: int,
    total_tasks: int,
) -> None:
    progress = (completed_tasks / total_tasks) * 100
    sys.stdout.write(f"\r          {math.floor(progress)}% complete")
    sys.stdout.flush()


def execute_phmmer_search(
    phmmer_cmds: List[str]
) -> None:
    total_tasks = len(phmmer_cmds)
    completed_tasks = 0

    for command in phmmer_cmds:
        if not check_if_phmmer_command_completed(command.split()[8]):
            run_bash_command(command)
        completed_tasks += 1
        update_progress(completed_tasks, total_tasks)


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
