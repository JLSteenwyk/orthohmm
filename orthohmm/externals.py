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
        # The MCL completion sentinel is the citation footer line. We compare
        # against the stripped form because mcl prints the line with leading
        # whitespace that varies between versions.
        if lines and lines[-1].strip() == "( http://link.aip.org/link/?SJMAEL/30/121/1 )":
            return True
    return False


def execute_leiden(
    cpm_resolution,
    output_directory: str,
    edges=None,
    include_isolates: bool = False,
    seed: int = 0,
) -> None:
    """In-process Leiden CPM clustering on the ABC edge file.

    Beats MCL on the OrthoBench reference (F=65.7% at resolution=0.1 vs MCL's
    best F=62.4% at inflation=1.5 on the identical RBNH edge graph) and has
    no external binary dependency — both libraries ship as wheels via pip.
    Output format matches what `execute_mcl` produces so downstream parsing
    in helpers.generate_orthogroup_clusters_file is unchanged.

    cpm_resolution may be a positive float, or the string "auto" to set γ
    to 4 × the smallest positive edge weight. Auto adapts to input diversity:
    closely-related species land near γ≈0.1, cross-kingdom eukaryotic inputs
    near γ≈0.00004 — without manual tuning.
    """
    import igraph
    import leidenalg
    import numpy as np

    edges_path = f"{output_directory}/orthohmm_working_res/orthohmm_edges.txt"
    out_path = f"{output_directory}/orthohmm_working_res/orthohmm_edges_clustered.txt"
    use_auto_resolution = (
        isinstance(cpm_resolution, str)
        and cpm_resolution.lower() == "auto"
    )

    # Load legacy ABC input or consume compact production edge arrays directly.
    name_to_id = {}
    graph_edges = []
    weights = []
    id_to_name = None
    min_positive_weight = None
    positive_weight_count = 0
    if edges is not None and hasattr(edges, "sources"):
        if include_isolates:
            used_ids = np.arange(len(edges.gene_names), dtype=np.int32)
        else:
            used_ids = np.unique(np.concatenate((edges.sources, edges.targets)))
        local_sources = np.searchsorted(used_ids, edges.sources).astype(np.int32)
        local_targets = np.searchsorted(used_ids, edges.targets).astype(np.int32)
        graph_edges = np.column_stack((local_sources, local_targets))
        weights = np.asarray(edges.weights, dtype=np.float64)
        id_to_name = [edges.gene_names[int(gene_id)] for gene_id in used_ids]
        if use_auto_resolution:
            positive = weights[weights > 0]
            positive_weight_count = len(positive)
            if len(positive) > 0:
                min_positive_weight = float(positive.min())
    else:
        with open(edges_path, "r") as f:
            for line in f:
                parts = line.rstrip("\n").split()
                if len(parts) < 3:
                    continue
                a, b, w = parts[0], parts[1], float(parts[2])
                ai = name_to_id.setdefault(a, len(name_to_id))
                bi = name_to_id.setdefault(b, len(name_to_id))
                graph_edges.append((ai, bi))
                weights.append(w)
                if use_auto_resolution and w > 0:
                    positive_weight_count += 1
                    if min_positive_weight is None or w < min_positive_weight:
                        min_positive_weight = w
        id_to_name = [None] * len(name_to_id)
        for name, idx in name_to_id.items():
            id_to_name[idx] = name

    if use_auto_resolution:
        # γ = 4 × min(edge_weight) over edges with strictly positive weight.
        # The weakest credible RBNH edge anchors the resolution; γ slightly
        # above it lets even weak credible edges count toward keeping
        # clusters together. We filter out exact-zero weights — they're
        # numerical artifacts (Viterbi raw score 0) and would collapse
        # γ to 0, which causes Leiden to try to merge everything and
        # crash on large graphs.
        if min_positive_weight is not None:
            gamma = 4.0 * float(min_positive_weight)
        else:
            gamma = 0.1
        print(f"          CPM auto: γ={gamma:.6g} "
              f"(4 × min over {positive_weight_count:,} positive of "
              f"{len(weights):,} edges)",
              flush=True)
    else:
        gamma = float(cpm_resolution)

    g = igraph.Graph(n=len(id_to_name), edges=graph_edges, directed=False)
    g.es["weight"] = weights

    part = leidenalg.find_partition(
        g, leidenalg.CPMVertexPartition,
        weights="weight", resolution_parameter=gamma,
        seed=seed,
    )

    output_clusters = []
    for cluster in part:
        if cluster:
            output_clusters.append(sorted(id_to_name[i] for i in cluster))
    output_clusters.sort()
    with open(out_path, "w") as out:
        for cluster in output_clusters:
            out.write(" ".join(cluster) + "\n")
