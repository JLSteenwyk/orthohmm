"""E-value estimation from Viterbi scores.

Uses Karlin-Altschul statistics: E = K * m * n * exp(-lambda * S)

The E-value is only used as a coarse pre-filter (threshold 0.0001)
in OrthoHMM's downstream pipeline. The actual edge weights use
length-normalized, phylogenetically-corrected bit scores.
"""

import numpy as np
from numba import njit, prange, float64, int32

from .matrices import KA_PARAMS


@njit(cache=True)
def estimate_evalue(score, query_len, db_total_residues, lam, K):
    """Estimate E-value for a single hit.

    Parameters
    ----------
    score : int, Viterbi score
    query_len : int, query sequence length
    db_total_residues : int, total residues in target database
    lam : float, Karlin-Altschul lambda
    K : float, Karlin-Altschul K

    Returns
    -------
    float, estimated E-value
    """
    if score <= 0:
        return 1e10
    return K * query_len * db_total_residues * np.exp(-lam * score)


@njit(parallel=True, cache=True)
def batch_estimate_evalues(
    scores,              # int32 array
    query_lengths,       # int32 array, one per pair
    db_total_residues,   # int64
    lam,                 # float64
    K,                   # float64
):
    """Estimate E-values for a batch of hits.

    Parameters
    ----------
    scores : (M,) int32
    query_lengths : (M,) int32, query length for each pair
    db_total_residues : total residues in target species
    lam, K : Karlin-Altschul parameters

    Returns
    -------
    evalues : (M,) float64
    """
    M = len(scores)
    evalues = np.empty(M, dtype=np.float64)
    for i in prange(M):
        if scores[i] <= 0:
            evalues[i] = 1e10
        else:
            evalues[i] = (
                K * float64(query_lengths[i]) * float64(db_total_residues)
                * np.exp(-lam * float64(scores[i]))
            )
    return evalues
