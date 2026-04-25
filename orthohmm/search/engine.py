"""Search engine orchestrator.

Ties together pre-filtering, profile construction, Viterbi scoring,
and E-value estimation into a complete phmmer replacement.
"""

import itertools
import multiprocessing
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

import numpy as np

from .matrices import get_matrix, get_background_freqs, get_ka_params
from .sequences import SequenceStore, SpeciesSequences
from .profile import build_profiles_batch
from .prefilter import prefilter_candidates
from .viterbi import (
    batch_viterbi_score, batch_viterbi_c, batch_viterbi_multipair_c,
    is_c_available,
)
try:
    from .viterbi_cuda_ctypes import is_cuda_available, batch_viterbi_cuda
except Exception:  # CUDA build optional
    def is_cuda_available():
        return False
    batch_viterbi_cuda = None
from .evalue import batch_estimate_evalues


@dataclass
class SpeciesPairResults:
    """Search results for one species pair (queries from A, targets from B)."""
    query_species: str
    target_species: str
    target_names: np.ndarray   # U50 string array
    query_names: np.ndarray    # U50 string array
    evalues: np.ndarray        # float64
    scores: np.ndarray         # float64


def _adaptive_max_candidates(query_species, target_species, kmer_k, min_total_hits,
                             use_reduced_alphabet, base_mc=20):
    """Determine max_candidates based on k-mer coverage between species.

    Runs a lightweight Phase A probe (no extension scoring) to measure
    what fraction of queries share enough k-mers with the target database.
    Close species → high coverage → low mc. Distant → low coverage → high mc.
    """
    from .prefilter import (build_kmer_index, compute_freq_threshold,
                            REDUCED_ALPHA, REDUCED_ALPHA_SIZE, _remap_flat)
    from numba import int32

    if use_reduced_alphabet:
        alpha_size = REDUCED_ALPHA_SIZE
        target_flat = _remap_flat(target_species.flat_sequences, REDUCED_ALPHA,
                                  len(target_species.flat_sequences))
        query_flat = _remap_flat(query_species.flat_sequences, REDUCED_ALPHA,
                                 len(query_species.flat_sequences))
    else:
        alpha_size = 20
        target_flat = target_species.flat_sequences
        query_flat = query_species.flat_sequences

    kmer_offsets, kmer_entries, kmer_freqs = build_kmer_index(
        target_flat, target_species.offsets, target_species.lengths,
        k=kmer_k, alpha_size=alpha_size,
    )
    freq_thresh = compute_freq_threshold(kmer_freqs, target_species.num_sequences)

    # Quick probe: count how many queries have >= min_total_hits k-mer matches
    # Use the Numba Phase A on a sample of queries for speed
    from .prefilter import _score_query_two_stage
    N = query_species.num_sequences
    sample_size = min(N, 500)
    step = max(1, N // sample_size)

    queries_with_cands = 0
    for qi in range(0, N, step):
        q_off = int(query_species.offsets[qi])
        q_len = int(query_species.lengths[qi])
        q_seq = query_flat[q_off:q_off + q_len]
        cids, _ = _score_query_two_stage(
            q_seq, int32(q_len), int32(kmer_k), int32(alpha_size),
            kmer_offsets, kmer_entries, kmer_freqs, int32(freq_thresh),
            int32(target_species.num_sequences), int32(min_total_hits),
            int32(1), int32(10),
        )
        if len(cids) > 0:
            queries_with_cands += 1

    coverage = queries_with_cands / sample_size

    if coverage > 0.8:
        mc = base_mc            # close species (e.g., human-chimp)
    elif coverage > 0.5:
        mc = base_mc * 3        # moderate (e.g., human-mouse)
    elif coverage > 0.2:
        mc = base_mc * 5        # distant (e.g., human-zebrafish)
    else:
        mc = base_mc * 10       # very distant (e.g., human-worm)

    return mc, coverage


def search_species_pair(
    query_species: SpeciesSequences,
    target_species: SpeciesSequences,
    matrix_name: str,
    band_width: int = 64,
    kmer_k: int = 4,
    use_reduced_alphabet: bool = False,
    min_total_hits: int = 4,
    min_diag_hits: int = 1,
    diag_bin_width: int = 10,
    max_candidates_per_query: int = 0,
) -> SpeciesPairResults:
    """Execute profile HMM search for one species pair.

    Steps:
      1. Pre-filter candidates (adaptive mc based on species distance)
      2. Build profiles for all queries
      3. Batch Viterbi scoring on candidate pairs
      4. Compute E-values
      5. Return results

    If max_candidates_per_query=0 (default), mc is determined adaptively
    based on k-mer coverage between the species pair.
    """
    sub_matrix = get_matrix(matrix_name)
    bg_freqs = get_background_freqs(matrix_name)
    lam, K = get_ka_params(matrix_name)

    N_query = query_species.num_sequences
    N_target = target_species.num_sequences

    if N_query == 0 or N_target == 0:
        return SpeciesPairResults(
            query_species=query_species.species_file,
            target_species=target_species.species_file,
            target_names=np.empty(0, dtype="U50"),
            query_names=np.empty(0, dtype="U50"),
            evalues=np.empty(0, dtype=np.float64),
            scores=np.empty(0, dtype=np.float64),
        )

    # Step 1: Determine adaptive mc if not specified
    if max_candidates_per_query <= 0:
        adaptive_mc, coverage = _adaptive_max_candidates(
            query_species, target_species, kmer_k, min_total_hits,
            use_reduced_alphabet,
        )
        max_candidates_per_query = max(adaptive_mc, 20)

    # Step 2: Pre-filter candidates using k-mer index + extension re-ranking
    cand_ids, cand_offsets = prefilter_candidates(
        query_species, target_species,
        k=kmer_k,
        use_reduced_alphabet=use_reduced_alphabet,
        min_total_hits=min_total_hits,
        min_diag_hits=min_diag_hits,
        diag_bin_width=diag_bin_width,
        max_candidates_per_query=max_candidates_per_query,
        sub_matrix=sub_matrix,
    )

    # Build pairs array from prefilter results
    pairs_list = []
    for qi in range(N_query):
        start = int(cand_offsets[qi])
        end = int(cand_offsets[qi + 1])
        for j in range(start, end):
            pairs_list.append((qi, int(cand_ids[j])))

    if not pairs_list:
        return SpeciesPairResults(
            query_species=query_species.species_file,
            target_species=target_species.species_file,
            target_names=np.empty(0, dtype="U50"),
            query_names=np.empty(0, dtype="U50"),
            evalues=np.empty(0, dtype=np.float64),
            scores=np.empty(0, dtype=np.float64),
        )

    pairs = np.array(pairs_list, dtype=np.int32)

    # Step 2: Build profiles for all queries
    flat_match_emit, insert_emit, transitions, profile_offsets, profile_lengths = \
        build_profiles_batch(query_species, sub_matrix, bg_freqs)

    # Step 3: Batch Viterbi scoring. Routing strategy:
    #   - Pairs whose target length T <= CUDA_T_MAX go to the GPU (if
    #     available) — ~2.4x faster than CPU multipair.
    #   - Pairs with T > CUDA_T_MAX can't fit in shared memory, so they
    #     go to the CPU multipair kernel.
    #   - If CUDA isn't available, everything goes to multipair AVX2.
    #   - If AVX2 isn't available either, fall back to scalar C then Numba.
    # The partition means one very long protein can't force the whole
    # batch to the slower path.
    CUDA_T_MAX = 1998  # kernel shared-memory cap
    scores = np.empty(len(pairs), dtype=np.int32)

    # Determine per-pair target lengths and partition into GPU-eligible / CPU-only
    use_gpu = is_cuda_available()
    if use_gpu:
        pair_t_lens = target_species.lengths[pairs[:, 1]]
        gpu_mask = pair_t_lens <= CUDA_T_MAX
    else:
        gpu_mask = np.zeros(len(pairs), dtype=bool)

    n_gpu = int(gpu_mask.sum())
    n_cpu = len(pairs) - n_gpu

    def _run_cpu(sub_pairs):
        """CPU fallback: multipair AVX2, then scalar C, then Numba."""
        if not is_c_available():
            return batch_viterbi_score(
                flat_match_emit, insert_emit, transitions,
                profile_offsets, profile_lengths,
                target_species.flat_sequences, target_species.offsets,
                target_species.lengths, sub_pairs, band_width,
            )
        try:
            return batch_viterbi_multipair_c(
                flat_match_emit, insert_emit, transitions,
                profile_offsets, profile_lengths,
                target_species.flat_sequences, target_species.offsets,
                target_species.lengths, sub_pairs, band_width,
            )
        except RuntimeError:
            return batch_viterbi_c(
                flat_match_emit, insert_emit, transitions,
                profile_offsets, profile_lengths,
                target_species.flat_sequences, target_species.offsets,
                target_species.lengths, sub_pairs, band_width,
            )

    if n_gpu > 0:
        gpu_pairs = pairs[gpu_mask]
        try:
            gpu_scores = batch_viterbi_cuda(
                flat_match_emit, insert_emit, transitions,
                profile_offsets, profile_lengths,
                target_species.flat_sequences, target_species.offsets,
                target_species.lengths, gpu_pairs, band_width,
            ).astype(np.int32)
            scores[gpu_mask] = gpu_scores
        except RuntimeError:
            # GPU failed — route all pairs to CPU
            scores = _run_cpu(pairs)
            n_gpu = 0
            n_cpu = len(pairs)

    if n_cpu > 0 and n_gpu > 0:
        # Some pairs still need CPU processing (long targets)
        cpu_pairs = pairs[~gpu_mask]
        scores[~gpu_mask] = _run_cpu(cpu_pairs)
    elif n_cpu > 0 and n_gpu == 0 and not use_gpu:
        # No GPU available at all
        scores = _run_cpu(pairs)

    # Step 4: Compute E-values
    db_total_residues = int(target_species.lengths.sum())
    query_lens_for_pairs = np.array(
        [query_species.lengths[pairs[i, 0]] for i in range(len(pairs))],
        dtype=np.int32,
    )
    evalues = batch_estimate_evalues(
        scores, query_lens_for_pairs,
        db_total_residues, lam, K,
    )

    # Step 5: Build result arrays with names
    target_names = np.array(
        [target_species.ids[pairs[i, 1]] for i in range(len(pairs))],
        dtype="U50",
    )
    query_names = np.array(
        [query_species.ids[pairs[i, 0]] for i in range(len(pairs))],
        dtype="U50",
    )

    return SpeciesPairResults(
        query_species=query_species.species_file,
        target_species=target_species.species_file,
        target_names=target_names,
        query_names=query_names,
        evalues=evalues,
        scores=scores.astype(np.float64),
    )


def _filter_significant(r: SpeciesPairResults,
                        evalue_threshold: float) -> SpeciesPairResults:
    """Keep only hits with evalue < threshold. Reduces IPC payload ~20x.

    The prefilter emits up to max_candidates_per_query targets per query;
    most of those are not homologs and end up with evalue >> 1e-4. Pickling
    them back to the parent process is wasted bandwidth and RAM — the
    downstream code drops them anyway when it applies the same threshold.
    """
    mask = r.evalues < evalue_threshold
    if mask.all():
        return r
    return SpeciesPairResults(
        query_species=r.query_species,
        target_species=r.target_species,
        target_names=r.target_names[mask],
        query_names=r.query_names[mask],
        evalues=r.evalues[mask],
        scores=r.scores[mask],
    )


def _search_pair_worker(args):
    """Worker function for multiprocessing."""
    (query_file, target_file, fasta_directory,
     matrix_name, band_width, kmer_k, evalue_threshold) = args

    query_sp = SpeciesSequences.from_fasta(
        f"{fasta_directory}/{query_file}", query_file
    )
    target_sp = SpeciesSequences.from_fasta(
        f"{fasta_directory}/{target_file}", target_file
    )

    r = search_species_pair(
        query_sp, target_sp, matrix_name,
        band_width=band_width, kmer_k=kmer_k,
    )
    return _filter_significant(r, evalue_threshold)


def execute_builtin_search(
    files: List[str],
    fasta_directory: str,
    output_directory: str,
    cpu: int,
    substitution_matrix,
    band_width: int = 64,
    kmer_k: int = 5,
    evalue_threshold: float = 1e-4,
) -> Dict[Tuple[str, str], SpeciesPairResults]:
    """Replace execute_phmmer_search with built-in search engine.

    Parameters
    ----------
    files : list of FASTA filenames
    fasta_directory : path to directory containing FASTA files
    output_directory : path for working results
    cpu : number of parallel workers
    substitution_matrix : SubstitutionMatrix enum value
    band_width : Viterbi band width
    kmer_k : k-mer size for pre-filtering

    Returns
    -------
    dict mapping (query_file, target_file) -> SpeciesPairResults
    """
    matrix_name = substitution_matrix.value
    file_pairs = list(itertools.product(files, repeat=2))

    total_tasks = len(file_pairs)
    results = {}

    # Build worker args. Worker filters hits by evalue_threshold before
    # returning, which cuts IPC payload ~20x at typical max_candidates settings.
    worker_args = [
        (qf, tf, fasta_directory, matrix_name, band_width, kmer_k,
         evalue_threshold)
        for qf, tf in file_pairs
    ]

    # Drain as-completed instead of holding all AsyncResults — lets the OS
    # reclaim each worker's pickled payload as soon as the parent has it.
    pool = multiprocessing.Pool(processes=cpu)
    completed = 0
    try:
        for r in pool.imap_unordered(
            _search_pair_worker, worker_args, chunksize=1,
        ):
            results[(r.query_species, r.target_species)] = r
            completed += 1
            progress = (completed / total_tasks) * 100
            sys.stdout.write(f"\r          {progress:.2f}% complete")
            sys.stdout.flush()
    finally:
        pool.close()
        pool.join()
    sys.stdout.write("\n")

    return results


def results_to_phmmer_format(
    results: SpeciesPairResults,
    evalue_threshold: float,
) -> np.ndarray:
    """Convert search results to the same structured array format
    that read_and_filter_phmmer_output returns.

    Returns numpy array with dtype:
      [("target_name", "U50"), ("query_name", "U50"),
       ("evalue", float), ("score", float)]
    """
    if len(results.scores) == 0:
        dtype = [
            ("target_name", "U50"),
            ("query_name", "U50"),
            ("evalue", float),
            ("score", float),
        ]
        return np.array([], dtype=dtype)

    # Filter by evalue
    mask = results.evalues < evalue_threshold

    n = mask.sum()
    dtype = [
        ("target_name", "U50"),
        ("query_name", "U50"),
        ("evalue", float),
        ("score", float),
    ]
    out = np.empty(n, dtype=dtype)
    out["target_name"] = results.target_names[mask]
    out["query_name"] = results.query_names[mask]
    out["evalue"] = results.evalues[mask]
    out["score"] = results.scores[mask]

    return out
