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


def execute_leiden(
    cpm_resolution,
    output_directory: str,
) -> None:
    """In-process Leiden CPM clustering on the ABC edge file.

    Beats MCL on the OrthoBench reference (F=65.7% at resolution=0.1 vs MCL's
    best F=62.4% at inflation=1.5 on the identical RBNH edge graph) and has
    no external binary dependency — both libraries ship as wheels via pip.
    Output format matches what `execute_mcl` produces so downstream parsing
    in helpers.generate_orthogroup_clusters_file is unchanged.

    cpm_resolution may be a positive float, or the string "auto" to set γ
    to the 10th percentile of edge weights. Auto adapts to input diversity:
    closely-related species (bacteria) land near γ≈0.1, cross-kingdom
    eukaryotic inputs near γ≈0.001 — without manual tuning.
    """
    import igraph
    import leidenalg
    import numpy as np

    edges_path = f"{output_directory}/orthohmm_working_res/orthohmm_edges.txt"
    out_path = f"{output_directory}/orthohmm_working_res/orthohmm_edges_clustered.txt"

    # Load edges: ABC format (gene_a \t gene_b \t weight)
    name_to_id = {}
    edges = []
    weights = []
    with open(edges_path, "r") as f:
        for line in f:
            parts = line.rstrip("\n").split()
            if len(parts) < 3:
                continue
            a, b, w = parts[0], parts[1], float(parts[2])
            ai = name_to_id.setdefault(a, len(name_to_id))
            bi = name_to_id.setdefault(b, len(name_to_id))
            edges.append((ai, bi))
            weights.append(w)

    if isinstance(cpm_resolution, str) and cpm_resolution.lower() == "auto":
        # γ = 4 × min(edge_weight). The minimum surviving edge weight is
        # the weakest credible RBNH edge; γ slightly above it means even
        # the weakest credible edges still count toward keeping clusters
        # together, while pure noise (which would have weight ≪ min) is
        # filtered out by the RBNH step that came earlier.
        # Empirically beats both the fixed default (cpm=0.1) and a
        # 10th-percentile heuristic on both OrthoBench (homogeneous,
        # closely-related bilaterians) and the cross-kingdom Three
        # Kingdoms benchmark.
        if weights:
            gamma = 4.0 * float(min(weights))
        else:
            gamma = 0.1
        print(f"          CPM auto: γ={gamma:.5f} (4 × min of {len(weights):,} edge weights)",
              flush=True)
    else:
        gamma = float(cpm_resolution)

    g = igraph.Graph(n=len(name_to_id), edges=edges, directed=False)
    g.es["weight"] = weights

    part = leidenalg.find_partition(
        g, leidenalg.CPMVertexPartition,
        weights="weight", resolution_parameter=gamma,
    )

    id_to_name = [None] * len(name_to_id)
    for name, idx in name_to_id.items():
        id_to_name[idx] = name

    with open(out_path, "w") as out:
        for cluster in part:
            if not cluster:
                continue
            out.write(" ".join(sorted(id_to_name[i] for i in cluster)) + "\n")
