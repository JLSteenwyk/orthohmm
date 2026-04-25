"""Built-in center-star MSA construction, replacing the MAFFT subprocess.

Algorithm: pick a center sequence (the median-length member of the
orthogroup), align every other member to it with banded Needleman-Wunsch
(C/OpenMP), then project all pairwise alignments onto the center's
coordinate system to produce a consistent multiple alignment.

This is a standard approach for closely related sequences (which is what
orthogroups are by definition) and avoids the MAFFT dependency entirely.

For highly divergent protein families (multi-domain, large indels),
progressive alignment with a guide tree is more accurate — but that's
rarely the case for orthogroups within a single clade.
"""

from __future__ import annotations

import ctypes
from pathlib import Path as _Path
from typing import List, Optional, Tuple

import numpy as np

from .matrices import FULL_CHAR_TO_IDX


_pa_lib = None


def _load_pair_align():
    """Lazy-load the C pair_align.so shared library."""
    global _pa_lib
    if _pa_lib is not None:
        return _pa_lib
    so_path = _Path(__file__).resolve().parent / "csrc" / "pair_align.so"
    lib = ctypes.cdll.LoadLibrary(str(so_path))
    lib.batch_pair_align_c.restype = None
    lib.batch_pair_align_c.argtypes = [
        ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,  # seqs_flat, offsets, lengths
        ctypes.c_void_p,                                     # sub_matrix
        ctypes.c_int32, ctypes.c_int32,                     # gap_open, gap_extend
        ctypes.c_void_p, ctypes.c_int32,                    # pairs, M
        ctypes.c_int32,                                      # band_width
        ctypes.c_void_p, ctypes.c_void_p,                   # out_a, out_b
        ctypes.c_void_p, ctypes.c_int32,                    # out_lens, max_align_len
        ctypes.c_void_p,                                     # out_scores
    ]
    lib.pair_align_set_num_threads.restype = None
    lib.pair_align_set_num_threads.argtypes = [ctypes.c_int32]
    _pa_lib = lib
    return _pa_lib


def _encode_sequence(seq: str) -> np.ndarray:
    """Encode a protein string to uint8 with 20 for unknown."""
    out = np.empty(len(seq), dtype=np.uint8)
    for i, c in enumerate(seq):
        out[i] = FULL_CHAR_TO_IDX.get(ord(c), 20)
    return out


def center_star_msa(
    sequences: List[str],
    sub_matrix: np.ndarray,
    gap_open: int = -11,
    gap_extend: int = -1,
    band_width: int = 0,
    n_threads: int = 1,
) -> Optional[List[str]]:
    """Build a center-star MSA from a list of protein strings.

    Returns aligned sequences (with '-' as gap) in the input order, or
    None if the alignment could not be built.

    Steps
    -----
    1. Pick center = the member with median length.
    2. Align every other member to center using banded NW.
    3. Scan all pairwise alignments: for each center position, find the
       max number of inserted residues across all members. Insert that
       many columns into the output MSA.
    4. Project each member onto the merged coordinate system.
    """
    n = len(sequences)
    if n < 2:
        return sequences

    # Encode everything to uint8
    encoded = [_encode_sequence(s) for s in sequences]
    lens = np.array([len(e) for e in encoded], dtype=np.int32)

    # Center = median length (argmedian)
    sorted_idx = np.argsort(lens)
    center_idx = int(sorted_idx[n // 2])

    # Build flat array for C
    offsets = np.zeros(n, dtype=np.int64)
    offsets[1:] = np.cumsum(lens.astype(np.int64))[:-1]
    flat = np.concatenate(encoded).astype(np.uint8)

    # Pairs: (other, center) for each other != center
    other_idx = [i for i in range(n) if i != center_idx]
    if not other_idx:
        return [sequences[center_idx]]
    pairs = np.array([(i, center_idx) for i in other_idx], dtype=np.int32)
    M = len(pairs)

    # Allocate output buffers. Worst case alignment length = la + lb.
    max_align_len = int(lens.max() * 2 + 4)

    out_a = np.zeros(M * max_align_len, dtype=np.uint8)
    out_b = np.zeros(M * max_align_len, dtype=np.uint8)
    out_lens = np.zeros(M, dtype=np.int32)
    out_scores = np.zeros(M, dtype=np.int32)

    lib = _load_pair_align()
    lib.pair_align_set_num_threads(n_threads)

    _sm = np.ascontiguousarray(sub_matrix, dtype=np.int8)
    _pr = np.ascontiguousarray(pairs.flatten(), dtype=np.int32)

    lib.batch_pair_align_c(
        flat.ctypes.data, offsets.ctypes.data, lens.ctypes.data,
        _sm.ctypes.data,
        int(gap_open), int(gap_extend),
        _pr.ctypes.data, M,
        int(band_width),
        out_a.ctypes.data, out_b.ctypes.data,
        out_lens.ctypes.data, int(max_align_len),
        out_scores.ctypes.data,
    )

    # Project: for each center position (0..L_center), determine how many
    # "inserted" (non-center) residues appear there across all pairs.
    L_center = int(lens[center_idx])

    # For each pair, walk the alignment and compute insert counts per center pos.
    # center_gaps[k] = number of insertions BEFORE center residue k (or after last, for k=L_center)
    max_insert = np.zeros(L_center + 1, dtype=np.int32)
    pair_walks: List[List[Tuple[int, int]]] = []  # per pair: list of (center_pos_before, kind)
    for p in range(M):
        L = int(out_lens[p])
        if L == 0:
            return None  # alignment overflow — give up gracefully
        a = out_a[p * max_align_len : p * max_align_len + L]
        b = out_b[p * max_align_len : p * max_align_len + L]
        center_col = 0  # position in center seq
        inserts_here = 0
        walk = []
        for col in range(L):
            if b[col] == 20:
                # center has gap here → 'other' residue is an insertion
                inserts_here += 1
                walk.append((center_col, 1))  # 1 = insertion column (other has AA)
            else:
                # center has a residue
                if inserts_here > max_insert[center_col]:
                    max_insert[center_col] = inserts_here
                inserts_here = 0
                walk.append((center_col, 0))  # 0 = match column
                center_col += 1
        # Handle trailing insertions (past end of center)
        if inserts_here > max_insert[center_col]:
            max_insert[center_col] = inserts_here
        pair_walks.append(walk)

    # Merged MSA length = L_center + sum(max_insert)
    msa_len = L_center + int(max_insert.sum())

    # Cumulative insert positions: insertion slot before center residue k
    # starts at column = k + sum(max_insert[:k])
    insert_start = np.zeros(L_center + 1, dtype=np.int32)
    for k in range(1, L_center + 1):
        insert_start[k] = insert_start[k - 1] + 1 + max_insert[k - 1]
    # center residue k lands at column: insert_start[k]

    # Build aligned sequences
    aligned = [None] * n
    GAP_CHAR = ord('-')

    # Center sequence: residues at fixed columns, gaps everywhere else
    center_enc = encoded[center_idx]
    center_aligned = np.full(msa_len, GAP_CHAR, dtype=np.uint8)
    for k in range(L_center):
        # Convert encoded byte back to ASCII
        aa = center_enc[k]
        center_aligned[insert_start[k]] = _aa_to_char(aa)
    aligned[center_idx] = bytes(center_aligned).decode("ascii")

    # Each non-center member
    for w_idx, pair_idx in enumerate(other_idx):
        L = int(out_lens[w_idx])
        a = out_a[w_idx * max_align_len : w_idx * max_align_len + L]
        b = out_b[w_idx * max_align_len : w_idx * max_align_len + L]
        walk = pair_walks[w_idx]

        member_aligned = np.full(msa_len, GAP_CHAR, dtype=np.uint8)

        # Track: for each center residue we've seen, how many insertions did we write?
        center_col = 0
        insert_count = 0
        for col in range(L):
            a_aa = a[col]
            if walk[col][1] == 0:
                # Match column: place a_aa at insert_start[center_col]
                if a_aa < 20:
                    member_aligned[insert_start[center_col]] = _aa_to_char(a_aa)
                # else: a had a gap here — leave as '-'
                center_col += 1
                insert_count = 0
            else:
                # Insertion column: place a_aa at insert_start[center_col] + 1 + insert_count
                slot = insert_start[center_col] - max_insert[center_col] + insert_count
                # insert_start[center_col] points to center residue col; insertions
                # sit in columns [insert_start[center_col] - max_insert[center_col], insert_start[center_col])
                if 0 <= slot < msa_len and a_aa < 20:
                    member_aligned[slot] = _aa_to_char(a_aa)
                insert_count += 1
        aligned[pair_idx] = bytes(member_aligned).decode("ascii")

    # Every position should be filled. Anything still None → fallback.
    for k in range(n):
        if aligned[k] is None:
            return None

    return aligned


# ASCII codes for AA alphabet ACDEFGHIKLMNPQRSTVWY
_AA_CHARS = np.array(
    [ord(c) for c in "ACDEFGHIKLMNPQRSTVWY"] + [ord("X")] * 12,
    dtype=np.uint8,
)


def _aa_to_char(aa_idx: int) -> int:
    if 0 <= aa_idx < 32:
        return int(_AA_CHARS[aa_idx])
    return ord("X")
