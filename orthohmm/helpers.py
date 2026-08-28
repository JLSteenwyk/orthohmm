from collections import defaultdict
from dataclasses import dataclass
from enum import Enum
import itertools
import multiprocessing
import os
import sys
from typing import (
    Tuple,
    List,
    DefaultDict,
    Dict,
)

import numpy as np

from .files import (
    write_fasta_files_for_all_ogs,
    write_fasta_files_for_single_copy_orthologs,
    write_file_of_single_copy_ortholog_names,
    write_clusters_file,
    write_copy_number_file,
)


class StopStep(Enum):
    prepare = "prepare"
    infer = "infer"
    write = "write"


class StartStep(Enum):
    search_res = "search_res"


class SubstitutionMatrix(Enum):
    blosum45 = "BLOSUM45"
    blosum50 = "BLOSUM50"
    blosum62 = "BLOSUM62"
    blosum80 = "BLOSUM80"
    blosum90 = "BLOSUM90"
    pam30 = "PAM30"
    pam70 = "PAM70"
    pam120 = "PAM120"
    pam240 = "PAM240"
    wag = "WAG"
    lg = "LG"


@dataclass
class IndexedEdges:
    """Compact undirected weighted graph keyed by global gene indices."""

    gene_names: List[str]
    sources: np.ndarray
    targets: np.ndarray
    weights: np.ndarray

    def __len__(self) -> int:
        return len(self.sources)


@dataclass
class IndexedEdgeThresholds:
    """Per-gene edge thresholds aligned with the global gene table."""

    values: np.ndarray


def write_indexed_edges(edges: IndexedEdges, output_directory: str) -> None:
    path = f"{output_directory}/orthohmm_working_res/orthohmm_edges.txt"
    with open(path, "w") as handle:
        for source, target, weight in zip(
            edges.sources, edges.targets, edges.weights
        ):
            handle.write(
                f"{edges.gene_names[int(source)]}\t"
                f"{edges.gene_names[int(target)]}\t{float(weight)}\n"
            )


def _write_indexed_edges(edges: IndexedEdges, output_directory: str) -> None:
    """Backward-compatible alias for the former private writer."""
    write_indexed_edges(edges, output_directory)


def _determine_indexed_network_edges(
    search_results,
    gene_lengths: np.ndarray,
    pairwise_rbh_corr: Dict[frozenset, np.float64],
    reciprocal_best_hit_thresholds: Dict[np.str_, np.float64],
    evalue_threshold: float,
    output_directory: str,
) -> IndexedEdges:
    """Build and deduplicate numeric edges without per-edge Python objects."""
    gene_names = [str(row["name"]) for row in gene_lengths]
    gene_to_id = {name: idx for idx, name in enumerate(gene_names)}
    species_global_ids = {
        species: np.fromiter(
            (gene_to_id[gene] for gene in names),
            dtype=np.int32,
            count=len(names),
        )
        for species, names in search_results.species_ids.items()
    }
    if isinstance(reciprocal_best_hit_thresholds, IndexedEdgeThresholds):
        thresholds = reciprocal_best_hit_thresholds.values
        if len(thresholds) != len(gene_names):
            raise ValueError("indexed edge thresholds do not match gene table")
    else:
        thresholds = np.full(len(gene_names), np.nan, dtype=np.float64)
        for gene, threshold in reciprocal_best_hit_thresholds.items():
            gene_id = gene_to_id.get(str(gene))
            if gene_id is not None:
                thresholds[gene_id] = float(threshold)

    source_chunks = []
    target_chunks = []
    weight_chunks = []
    for pair, result in search_results.items():
        mask = result.evalues < evalue_threshold
        if not mask.any():
            continue
        correction = float(pairwise_rbh_corr[frozenset(pair)])
        query_ids = species_global_ids[pair[0]][result.query_indices[mask]]
        target_ids = species_global_ids[pair[1]][result.target_indices[mask]]
        with np.errstate(divide="ignore", invalid="ignore"):
            weights = result.scores[mask] / correction
        keep = (
            ~np.isnan(thresholds[query_ids])
            & (weights >= thresholds[query_ids])
            & (query_ids != target_ids)
        )
        if keep.any():
            source_chunks.append(np.minimum(query_ids[keep], target_ids[keep]))
            target_chunks.append(np.maximum(query_ids[keep], target_ids[keep]))
            weight_chunks.append(weights[keep])

    if not source_chunks:
        edges = IndexedEdges(
            gene_names,
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.float64),
        )
        write_indexed_edges(edges, output_directory)
        return edges

    sources = np.concatenate(source_chunks).astype(np.int32, copy=False)
    targets = np.concatenate(target_chunks).astype(np.int32, copy=False)
    weights = np.concatenate(weight_chunks).astype(np.float64, copy=False)
    packed = (sources.astype(np.int64) << 32) | targets.astype(np.int64)
    order = np.argsort(packed, kind="stable")
    packed_sorted = packed[order]
    weights_sorted = weights[order]
    starts = np.concatenate(
        (np.array([0], dtype=np.int64), np.flatnonzero(np.diff(packed_sorted)) + 1)
    )
    unique_packed = packed_sorted[starts]
    unique_weights = np.maximum.reduceat(weights_sorted, starts)
    edges = IndexedEdges(
        gene_names,
        (unique_packed >> 32).astype(np.int32),
        (unique_packed & 0xFFFFFFFF).astype(np.int32),
        unique_weights,
    )
    write_indexed_edges(edges, output_directory)
    return edges


def generate_phmmer_cmds(
    files: List[str],
    phmmer: str,
    output_directory: str,
    fasta_directory: str,
    cpu: int,
    substitution_matrix: SubstitutionMatrix,
) -> List[str]:
    pairwise_combos = list(itertools.product(files, repeat=2))
    phmmer_cmds = []
    for combo in pairwise_combos:
        phmmer_cmds.append(f"{phmmer} --mx {substitution_matrix.value} --noali --notextw --cpu {cpu} --tblout {output_directory}/orthohmm_working_res/{combo[0]}_2_{combo[1]}.phmmerout.txt {fasta_directory}/{combo[0]} {fasta_directory}/{combo[1]}")

    return phmmer_cmds


def get_sequence_lengths(
    fasta_directory: str,
    files: List[str],
) -> np.ndarray:
    gene_lengths = list()

    for file in files:
        ids = set()
        with open(os.path.join(fasta_directory, file), "r") as fasta_file:
            sequence_id = None
            sequence_length = 0

            for line in fasta_file:
                line = line.strip()
                if line.startswith(">"):
                    if sequence_id:
                        gene_lengths.append(
                            (file, sequence_id, sequence_length)
                        )
                    sequence_id = line[1:].split()[0]
                    if sequence_id in ids:
                        print(f"{sequence_id} appears twice in file {file}. Please ensure all FASTA headers are unique within each file.")
                        sys.exit(1)
                    ids.add(sequence_id)
                    sequence_length = 0
                else:
                    sequence_length += len(line)

            if sequence_id:
                gene_lengths.append((file, sequence_id, sequence_length))

    dtype = [
        ("spp", object),
        ("name", object),
        ("length", int)
    ]
    return np.array(gene_lengths, dtype=dtype)


def merge_with_gene_lengths(
    res: np.ndarray,
    gene_lengths: np.ndarray,
) -> np.ndarray:
    length_dict = {name: length for _, name, length in gene_lengths}

    res_merged = np.empty((len(res), 6), dtype=object)

    res_merged[:, 0] = res["target_name"]
    res_merged[:, 1] = res["query_name"]
    res_merged[:, 2] = res["evalue"]
    res_merged[:, 3] = res["score"]

    res_merged[:, 4] = [
        length_dict.get(name, None) for name in res["target_name"]
    ]
    res_merged[:, 5] = [
        length_dict.get(name, None) for name in res["query_name"]
    ]

    return res_merged


def read_and_filter_phmmer_output(
    taxon_a: str,
    taxon_b: str,
    output_directory: str,
    evalue_threshold: float,
    search_results=None,
) -> np.ndarray:
    # If in-memory search results are available, use them
    if search_results is not None:
        from .search.engine import (
            IndexedSearchResults,
            indexed_results_to_phmmer_format,
            results_to_phmmer_format,
        )
        pair_results = search_results.get((taxon_a, taxon_b))
        if pair_results is None:
            dtype_res = [
                ("target_name", "U50"),
                ("query_name", "U50"),
                ("evalue", float),
                ("score", float)
            ]
            return np.array([], dtype=dtype_res)
        if isinstance(search_results, IndexedSearchResults):
            return indexed_results_to_phmmer_format(
                pair_results,
                search_results.species_ids[taxon_a],
                search_results.species_ids[taxon_b],
                evalue_threshold,
            )
        return results_to_phmmer_format(pair_results, evalue_threshold)

    res_path = f"{output_directory}/orthohmm_working_res/{taxon_a}_2_{taxon_b}.phmmerout.txt"

    dtype_res = [
        ("target_name", "U50"),
        ("query_name", "U50"),
        ("evalue", float),
        ("score", float)
    ]

    res = np.genfromtxt(
        res_path,
        comments="#",
        dtype=dtype_res,
        usecols=[0, 2, 4, 5],
        encoding="utf-8"
    )

    res = res[res["evalue"] < evalue_threshold]

    return res


def normalize_by_gene_length(res_merged: np.ndarray) -> np.ndarray:
    # Sum-of-lengths normalization for phmmer bit scores. Phmmer's bit
    # scores already include implicit length correction via Karlin–Altschul,
    # so dividing by sqrt(L_q*L_t) over-normalizes and degrades clustering.
    # The built-in Viterbi engine pre-normalizes its raw integer scores
    # inside results_to_phmmer_format (geometric-mean) and skips this step.
    res_merged[:, 3] = res_merged[:, 3] / (res_merged[:, 4] + res_merged[:, 5])

    return res_merged


def correct_by_phylogenetic_distance(
    best_hits_A_to_B: Dict[np.str_, Dict[np.str_, np.float64]],
    best_hits_B_to_A: Dict[np.str_, Dict[np.str_, np.float64]],
    pair: Tuple[str, str],
    pairwise_rbh_corr: Dict[frozenset, np.float64],
) -> Tuple[
    Dict[np.str_, np.float64],
    Dict[np.str_, np.float64],
    Dict[frozenset, np.float64],
]:
    pair_set = frozenset(pair)
    # get rbh scores
    rbh_scores = []
    rbh_pairs_identified = 0
    for query, best_hit in best_hits_A_to_B.items():
        target = best_hit["target"]
        # Check if the reciprocal hit is also the best hit
        if target in best_hits_B_to_A and best_hits_B_to_A[target]["target"] == query:
            score_between_rbh_pair = (
                (best_hit["score"] + best_hits_B_to_A[target]["score"]) / 2
            )
            rbh_scores.append(score_between_rbh_pair)
            rbh_pairs_identified += 1

    mean_rbh_scores = np.mean(rbh_scores) if rbh_scores else 0

    if pair_set not in pairwise_rbh_corr:
        pairwise_rbh_corr[pair_set] = mean_rbh_scores
    else:
        pairwise_rbh_corr[pair_set] = (pairwise_rbh_corr[pair_set] + mean_rbh_scores) / 2

    # Phylogenetic correction
    correction_factor = pairwise_rbh_corr[pair_set]
    best_hit_scores_A_to_B = {key: value["score"] / correction_factor for key, value in best_hits_A_to_B.items()}
    best_hit_scores_B_to_A = {key: value["score"] / correction_factor for key, value in best_hits_B_to_A.items()}

    return best_hit_scores_A_to_B, best_hit_scores_B_to_A, pairwise_rbh_corr


def get_best_hits_and_scores(
    res_merged: np.ndarray
) -> DefaultDict[
    np.str_,
    Dict[np.str_, np.float64],
]:
    """
    get dictionaries of scores for best hit and best hit
    """
    best_hits = defaultdict(dict)

    for record in res_merged:
        query = record[1]
        target = record[0]
        score = record[3]
        if query not in best_hits or best_hits[query]["score"] < score:
            best_hits[query] = {"target": target, "score": score}

    return best_hits


def get_threshold_per_gene(
    best_hits_A_to_B: Dict[np.str_, Dict[np.str_, np.float64]],
    best_hits_B_to_A: Dict[np.str_, Dict[np.str_, np.float64]],
    best_hit_scores_A_to_B: Dict[np.str_, np.float64],
    best_hit_scores_B_to_A: Dict[np.str_, np.float64],
    reciprocal_best_hit_thresholds: Dict[np.str_, np.float64],
) -> Dict[str, np.float64]:
    for gene_A, data_A in best_hits_A_to_B.items():
        gene_B = data_A["target"]
        data_B = best_hits_B_to_A.get(gene_B)
        if data_B is not None and data_B["target"] == gene_A:
            score = (best_hit_scores_A_to_B[gene_A] + best_hit_scores_B_to_A[gene_B]) / 2

            current_threshold = reciprocal_best_hit_thresholds.get(gene_A)

            if current_threshold is None or score < current_threshold:
                reciprocal_best_hit_thresholds[gene_A] = score

    return reciprocal_best_hit_thresholds


def process_pair_edge_thresholds(
    pair: Tuple[str, str],
    output_directory: str,
    gene_lengths: np.ndarray,
    evalue_threshold: float,
    search_results=None,
) -> Tuple[
    Dict[str, np.float64],
    Dict[frozenset, np.float64],
]:
    fwd_res = read_and_filter_phmmer_output(
        pair[0], pair[1],
        output_directory,
        evalue_threshold,
        search_results=search_results,
    )
    rev_res = read_and_filter_phmmer_output(
        pair[1], pair[0],
        output_directory,
        evalue_threshold,
        search_results=search_results,
    )

    fwd_res_merged = merge_with_gene_lengths(fwd_res, gene_lengths)
    rev_res_merged = merge_with_gene_lengths(rev_res, gene_lengths)

    # Built-in search results are pre-normalized by sqrt(L_q*L_t) inside the
    # engine; phmmer bit scores still need the legacy sum-of-lengths step.
    if search_results is None:
        fwd_res_merged = normalize_by_gene_length(fwd_res_merged)
        rev_res_merged = normalize_by_gene_length(rev_res_merged)

    best_hits_A_to_B = get_best_hits_and_scores(fwd_res_merged)
    best_hits_B_to_A = get_best_hits_and_scores(rev_res_merged)

    best_hit_scores_A_to_B, best_hit_scores_B_to_A, pairwise_rbh_corr = \
        correct_by_phylogenetic_distance(
            best_hits_A_to_B,
            best_hits_B_to_A,
            pair,
            {}
        )

    reciprocal_best_hit_thresholds = \
        get_threshold_per_gene(
            best_hits_A_to_B,
            best_hits_B_to_A,
            best_hit_scores_A_to_B,
            best_hit_scores_B_to_A,
            {}
        )

    return reciprocal_best_hit_thresholds, pairwise_rbh_corr


def _best_indexed_hits(result, n_queries, evalue_threshold):
    """Return the first maximum-scoring target for each local query index."""
    best_targets = np.full(n_queries, -1, dtype=np.int32)
    best_scores = np.full(n_queries, -np.inf, dtype=np.float64)
    if result is None or len(result.scores) == 0:
        best_scores.fill(np.nan)
        return best_targets, best_scores

    mask = result.evalues < evalue_threshold
    if not mask.any():
        best_scores.fill(np.nan)
        return best_targets, best_scores

    query_indices = result.query_indices[mask]
    target_indices = result.target_indices[mask]
    scores = result.scores[mask]
    np.maximum.at(best_scores, query_indices, scores)

    # np.unique returns the first matching row for each query, preserving the
    # legacy strict-maximum behavior when multiple targets have equal scores.
    winning_rows = np.flatnonzero(scores == best_scores[query_indices])
    winning_queries, first = np.unique(
        query_indices[winning_rows], return_index=True
    )
    selected_rows = winning_rows[first]
    best_targets[winning_queries] = target_indices[selected_rows]
    best_scores[best_targets == -1] = np.nan
    return best_targets, best_scores


def _determine_indexed_edge_thresholds(
    files,
    gene_lengths,
    evalue_threshold,
    search_results,
):
    """Calculate reciprocal-best-hit thresholds using compact integer arrays."""
    gene_names = [str(row["name"]) for row in gene_lengths]
    gene_to_id = {name: idx for idx, name in enumerate(gene_names)}
    species_global_ids = {
        species: np.fromiter(
            (gene_to_id[gene] for gene in names),
            dtype=np.int32,
            count=len(names),
        )
        for species, names in search_results.species_ids.items()
    }
    best_hits = {
        pair: _best_indexed_hits(
            search_results.get(pair),
            len(search_results.species_ids[pair[0]]),
            evalue_threshold,
        )
        for pair in itertools.product(files, repeat=2)
    }

    threshold_values = np.full(len(gene_names), np.nan, dtype=np.float64)
    threshold_initialized = np.zeros(len(gene_names), dtype=bool)
    pairwise_corr = {}
    file_pairs = list(itertools.product(files, repeat=2))

    for idx, pair in enumerate(file_pairs):
        fwd_targets, fwd_scores = best_hits[pair]
        rev_targets, rev_scores = best_hits[(pair[1], pair[0])]
        query_indices = np.flatnonzero(fwd_targets >= 0)
        target_indices = fwd_targets[query_indices]
        reciprocal = (
            (target_indices < len(rev_targets))
            & (rev_targets[target_indices] == query_indices)
        )
        reciprocal_queries = query_indices[reciprocal]
        reciprocal_targets = target_indices[reciprocal]

        if len(reciprocal_queries):
            pair_scores = (
                fwd_scores[reciprocal_queries]
                + rev_scores[reciprocal_targets]
            ) / 2.0
            correction = np.mean(pair_scores)
        else:
            correction = np.float64(0.0)

        pair_key = frozenset(pair)
        if pair_key in pairwise_corr:
            pairwise_corr[pair_key] = (
                pairwise_corr[pair_key] + correction
            ) / 2.0
        else:
            pairwise_corr[pair_key] = correction

        if len(reciprocal_queries):
            with np.errstate(divide="ignore", invalid="ignore"):
                new_values = (
                    fwd_scores[reciprocal_queries] / correction
                    + rev_scores[reciprocal_targets] / correction
                ) / 2.0
            global_ids = species_global_ids[pair[0]][reciprocal_queries]
            initialized = threshold_initialized[global_ids]
            current = threshold_values[global_ids]
            update = ~initialized | (new_values < current)
            threshold_values[global_ids[update]] = new_values[update]
            threshold_initialized[global_ids] = True

        progress = ((idx + 1) / len(file_pairs)) * 100
        sys.stdout.write(f"\r          {progress:.2f}% complete")
        sys.stdout.flush()

    return IndexedEdgeThresholds(threshold_values), pairwise_corr


def update_progress(
    lock,
    completed_tasks,
    total_tasks: int,
) -> None:
    with lock:
        completed_tasks.value += 1
        progress = (completed_tasks.value / total_tasks) * 100
        sys.stdout.write(f"\r          {progress:.2f}% complete")
        sys.stdout.flush()


def determine_edge_thresholds(
    files: List[str],
    fasta_directory: str,
    output_directory: str,
    cpu: int,
    evalue_threshold: float,
    search_results=None,
) -> Tuple[
    np.ndarray,
    Dict[str, float],
    Dict[frozenset, float],
]:
    gene_lengths = get_sequence_lengths(fasta_directory, files)
    file_pairs = [(file1, file2) for file1 in files for file2 in files]

    if search_results is not None:
        from .search.engine import IndexedSearchResults
        if isinstance(search_results, IndexedSearchResults):
            thresholds, pairwise_corr = _determine_indexed_edge_thresholds(
                files,
                gene_lengths,
                evalue_threshold,
                search_results,
            )
            return gene_lengths, thresholds, pairwise_corr

        # In-memory path: process sequentially (results already computed)
        results_list = []
        total_tasks = len(file_pairs)
        for idx, pair in enumerate(file_pairs):
            res = process_pair_edge_thresholds(
                pair, output_directory, gene_lengths,
                evalue_threshold, search_results=search_results,
            )
            results_list.append(res)
            progress = ((idx + 1) / total_tasks) * 100
            sys.stdout.write(f"\r          {progress:.2f}% complete")
            sys.stdout.flush()

        final_reciprocal_thresholds = {}
        final_pairwise_corr = {}

        for thresholds, pairwise_corr in results_list:
            for key, value in thresholds.items():
                if key in final_reciprocal_thresholds:
                    final_reciprocal_thresholds[key] = min(
                        final_reciprocal_thresholds[key], value
                    )
                else:
                    final_reciprocal_thresholds[key] = value
            for key, value in pairwise_corr.items():
                if key in final_pairwise_corr:
                    final_pairwise_corr[key] = (
                        final_pairwise_corr[key] + value
                    ) / 2
                else:
                    final_pairwise_corr[key] = value

        return gene_lengths, final_reciprocal_thresholds, final_pairwise_corr

    pool = multiprocessing.Pool(processes=cpu)
    completed_tasks = multiprocessing.Value("i", 0)
    total_tasks = len(file_pairs)
    lock = multiprocessing.Lock()

    results = [pool.apply_async(
        process_pair_edge_thresholds,
        args=(
            pair,
            output_directory,
            gene_lengths,
            evalue_threshold
        ),
        callback=lambda _: update_progress(
                lock, completed_tasks, total_tasks
            )
        )
        for pair in file_pairs
    ]

    pool.close()
    pool.join()

    final_reciprocal_thresholds = {}
    final_pairwise_corr = {}

    for res in results:
        thresholds, pairwise_corr = res.get()
        for key, value in thresholds.items():
            if key in final_reciprocal_thresholds:
                final_reciprocal_thresholds[key] = min(
                    final_reciprocal_thresholds[key], value
                )
            else:
                final_reciprocal_thresholds[key] = value
        for key, value in pairwise_corr.items():
            if key in final_pairwise_corr:
                final_pairwise_corr[key] = (
                    final_pairwise_corr[key] + value
                ) / 2
            else:
                final_pairwise_corr[key] = value

    return gene_lengths, final_reciprocal_thresholds, final_pairwise_corr


def process_pair_determine_network_edges(
    pair,
    output_directory,
    gene_lengths,
    pairwise_rbh_corr,
    reciprocal_best_hit_thresholds,
    evalue_threshold,
    search_results=None,
):
    edges = {}
    res = read_and_filter_phmmer_output(
        pair[0], pair[1], output_directory, evalue_threshold,
        search_results=search_results,
    )

    # Phmmer bit scores: divide by (L_q + L_t) — Karlin–Altschul-style sum
    # normalization. Built-in Viterbi raw scores are pre-normalized by
    # sqrt(L_q*L_t) inside results_to_phmmer_format, so skip the divisor here.
    use_builtin = search_results is not None
    for hit in res:
        query_length = gene_lengths[hit["query_name"]]
        target_length = gene_lengths[hit["target_name"]]
        if use_builtin:
            length_norm = 1.0  # already applied in the search engine
        else:
            length_norm = float(query_length + target_length)
        norm_score = (hit["score"] / length_norm) / pairwise_rbh_corr[frozenset(pair)]

        try:
            if norm_score >= reciprocal_best_hit_thresholds[hit["query_name"]]:
                genes = frozenset([hit["query_name"], hit["target_name"]])
                if len(genes) == 2:
                    if genes in edges and edges[genes] < norm_score:
                        edges[genes] = norm_score
                    elif genes not in edges:
                        edges[genes] = norm_score
        except KeyError:
            continue
    return edges


def determine_network_edges(
    files: List[str],
    output_directory: str,
    gene_lengths: np.ndarray,
    pairwise_rbh_corr: Dict[frozenset, np.float64],
    reciprocal_best_hit_thresholds: Dict[np.str_, np.float64],
    evalue_threshold: float,
    cpu: int,
    search_results=None,
) -> Dict[frozenset, np.float64]:

    from .search.engine import IndexedSearchResults
    if isinstance(search_results, IndexedSearchResults):
        return _determine_indexed_network_edges(
            search_results,
            gene_lengths,
            pairwise_rbh_corr,
            reciprocal_best_hit_thresholds,
            evalue_threshold,
            output_directory,
        )

    gene_lengths = {str(row["name"]): int(row["length"]) for row in gene_lengths}
    file_pairs = [(file1, file2) for file1 in files for file2 in files]

    if search_results is not None:
        # In-memory path: process sequentially
        edges = {}
        total_tasks = len(file_pairs)
        for idx, pair in enumerate(file_pairs):
            edge_result = process_pair_determine_network_edges(
                pair, output_directory, gene_lengths,
                pairwise_rbh_corr, reciprocal_best_hit_thresholds,
                evalue_threshold, search_results=search_results,
            )
            for key, value in edge_result.items():
                if key in edges:
                    edges[key] = max(edges[key], value)
                else:
                    edges[key] = value
            progress = ((idx + 1) / total_tasks) * 100
            sys.stdout.write(f"\r          {progress:.2f}% complete")
            sys.stdout.flush()

        with open(f"{output_directory}/orthohmm_working_res/orthohmm_edges.txt", "w") as file:
            for key, value in edges.items():
                key_str = "\t".join(map(str, key))
                file.write(f"{key_str}\t{value}\n")

        return edges

    total_tasks = len(file_pairs)
    completed_tasks = multiprocessing.Value("i", 0)
    lock = multiprocessing.Lock()

    pool = multiprocessing.Pool(processes=cpu)

    results = []
    for pair in file_pairs:
        result = pool.apply_async(
            process_pair_determine_network_edges,
            args=(
                pair,
                output_directory,
                gene_lengths,
                pairwise_rbh_corr,
                reciprocal_best_hit_thresholds,
                evalue_threshold,
            ),
            callback=lambda _: update_progress(
                lock, completed_tasks, total_tasks
            )
        )
        results.append(result)

    pool.close()
    pool.join()

    edges = {}
    for result in results:
        edge_result = result.get()
        for key, value in edge_result.items():
            if key in edges:
                edges[key] = max(edges[key], value)
            else:
                edges[key] = value

    with open(f"{output_directory}/orthohmm_working_res/orthohmm_edges.txt", "w") as file:
        for key, value in edges.items():
            key_str = "\t".join(map(str, key))
            file.write(f"{key_str}\t{value}\n")

    return edges


def get_singletons(
    gene_lengths: np.ndarray,
    clustering_res: List[List[str]],
) -> Tuple[
    List[List[str]],
    List[List[str]],
]:
    singletons = sorted(
        set(gene_lengths["name"]) - set([j for i in clustering_res for j in i])
    )
    singletons = [[str(i)] for i in singletons]
    clustering_res.extend(singletons)

    return clustering_res, singletons


def get_all_fasta_entries(
    fasta_directory: str,
    files: List[str],
) -> Dict[str, str]:
    entries = {}
    for fasta_file in files:
        fasta_file_entries = {}

        header = None
        sequence = []

        with open(f"{fasta_directory}/{fasta_file}", "r") as file:
            for line in file:
                line = line.strip()
                if line.startswith(">"):
                    if header:
                        fasta_file_entries[header] = "".join(sequence)
                    header = line[1:].split()[0]  # remove the ">" character
                    sequence = []
                else:
                    sequence.append(line)

            # don't forget to add the last entry to the dictionary
            if header:
                fasta_file_entries[header] = "".join(sequence)

        entries[fasta_file] = fasta_file_entries

    return entries


def get_orthogroup_information(
    files: List[str],
    gene_lengths: np.ndarray,
    clustering_res: List[List[str]],
    single_copy_threshold: float,
    entries: Dict[str, str],
) -> Tuple[
    List[List[str]],
    Dict[str, List[str]],
    Dict[str, List[str]],
    List[str],
]:
    ogs_dat = dict()
    og_cn = dict()
    single_copy_ogs = list()

    total_ogs = len(clustering_res)
    width = len(str(total_ogs))
    total_number_of_taxa = len(np.unique(gene_lengths["spp"]))
    gene_to_row = {
        str(row["name"]): row_idx
        for row_idx, row in enumerate(gene_lengths)
    }

    og_cn["files:"] = files

    for i in range(total_ogs):
        # Format the OG number with leading zeros
        og_id = f"OG{i:0{width}}:"
        row_indices = sorted(
            gene_to_row[gene]
            for gene in clustering_res[i]
            if gene in gene_to_row
        )
        og_rows = gene_lengths[np.asarray(row_indices, dtype=np.intp)]
        # test if single-copy
        unique_species = np.unique(og_rows["spp"])
        if len(unique_species) == len(og_rows["spp"]):
            # test if sufficient occupancy
            if len(unique_species) / total_number_of_taxa > single_copy_threshold:
                single_copy_ogs.append(f"OG{i}")
        og_dat = list()
        for row in og_rows:
            lines = [entries[row["spp"]][row["name"]][i:i+70] for i in range(0, len(entries[row["spp"]][row["name"]]), 70)]
            og_dat.append(f">{row['name']}")
            og_dat.extend(lines)
        ogs_dat[f"OG{i}"] = og_dat

        spp_values, counts = np.unique(og_rows["spp"], return_counts=True)
        spp_counts = dict(zip(spp_values, counts))

        cnts = []
        for file in files:
            try:
                cnts.append(str(spp_counts.get(file, 0)))
            except IndexError:
                cnts.append("0")
        og_cn[og_id] = cnts

        clustering_res[i].insert(0, og_id)

    return clustering_res, og_cn, ogs_dat, single_copy_ogs


def generate_orthogroup_clusters_file(
    output_directory: str,
    gene_lengths: np.ndarray,
    files: List[str],
    single_copy_threshold: float,
    fasta_directory: str,
) -> Tuple[
    List[List[str]],
    Dict[str, List[str]],
    Dict[str, List[str]],
    List[str],
]:
    clustering_res = list()

    entries = get_all_fasta_entries(fasta_directory, files)

    with open(
        f"{output_directory}/orthohmm_working_res/orthohmm_edges_clustered.txt",
        "r",
    ) as file:
        for line in file:
            line = line.strip()
            if line:
                clustering_res.append(line.split())

    # get singletons - i.e., genes that aren't in groups with other genes
    clustering_res, singletons = get_singletons(
        gene_lengths,
        clustering_res,
    )

    clustering_res, og_cn, ogs_dat, single_copy_ogs = \
        get_orthogroup_information(
            files,
            gene_lengths,
            clustering_res,
            single_copy_threshold,
            entries,
        )

    write_clusters_file(output_directory, clustering_res)

    return singletons, og_cn, ogs_dat, single_copy_ogs


def generate_orthogroup_files(
    output_directory: str,
    gene_lengths: np.ndarray,
    og_cn: Dict[str, List[str]],
    ogs_dat: Dict[str, List[str]],
    single_copy_ogs: List[str],
) -> None:
    extensions = (".fa", ".faa", ".fas", ".fasta", ".pep", ".prot")
    write_copy_number_file(output_directory, og_cn)
    write_file_of_single_copy_ortholog_names(output_directory, og_cn)
    write_fasta_files_for_all_ogs(output_directory, ogs_dat)
    write_fasta_files_for_single_copy_orthologs(
        output_directory, ogs_dat, gene_lengths, single_copy_ogs, extensions
    )
