"""Opt-in position-based profile weights and experimental effective counts."""

from numbers import Real

import numpy as np


PROFILE_SEQUENCE_WEIGHTING = (
    "none", "henikoff", "henikoff_cluster", "species_henikoff_cluster",
)


def validate_profile_sequence_weighting(mode):
    if not isinstance(mode, str) or mode not in PROFILE_SEQUENCE_WEIGHTING:
        raise ValueError(
            "Unknown profile sequence weighting; expected one of: "
            + ", ".join(PROFILE_SEQUENCE_WEIGHTING)
        )
    return mode


def _alignment(encoded_msa):
    msa = np.asarray(encoded_msa)
    if (msa.ndim != 2 or not np.issubdtype(msa.dtype, np.integer)
            or np.any(msa < 0) or np.any(msa > 20)):
        raise ValueError("Expected a 2D integer alignment with symbols 0 through 20")
    return msa


def position_based_weights(encoded_msa):
    """Henikoff contributions, normalized for each row's nonmissing length.

    Gaps/unknowns (20) contribute nothing. The result is relative, not yet
    normalized to an effective sequence count. Empty rows have zero weight.
    """
    msa = _alignment(encoded_msa)
    weights = np.zeros(len(msa), dtype=np.float64)
    lengths = (msa < 20).sum(axis=1)
    for column in msa.T:
        valid = column < 20
        counts = np.bincount(column[valid], minlength=20)
        distinct = np.count_nonzero(counts)
        if distinct:
            weights[valid] += 1.0 / (distinct * counts[column[valid]])
    np.divide(weights, lengths, out=weights, where=lengths > 0)
    return weights


def effective_cluster_count(encoded_msa, identity=0.62):
    """Count single-linkage components using aligned matches / shorter length.

    Only valid residues enter matches or lengths. This experimental count
    is not entropy weighting or a reproduction of all HMMER MSA semantics.
    """
    msa = _alignment(encoded_msa)
    if (isinstance(identity, (bool, np.bool_)) or not isinstance(identity, Real)
            or not np.isfinite(identity) or not 0 < identity <= 1):
        raise ValueError("Identity must be finite and in (0, 1]")
    lengths = (msa < 20).sum(axis=1)
    usable = np.flatnonzero(lengths > 0)
    parents = np.arange(len(msa))

    def root(index):
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    for offset, left in enumerate(usable):
        right = usable[offset + 1:]
        matches = ((msa[right] == msa[left]) & (msa[left] < 20)).sum(axis=1)
        linked = matches >= identity * np.minimum(lengths[left], lengths[right])
        for neighbor in right[linked]:
            parents[root(int(neighbor))] = root(int(left))
    return len({root(int(index)) for index in usable})


def profile_sequence_weights(encoded_msa, mode, species=None):
    """Return row weights summing to actual N or the specified effective N."""
    validate_profile_sequence_weighting(mode)
    msa = _alignment(encoded_msa)
    if mode == "none":
        return np.ones(len(msa), dtype=np.float64)
    weights = position_based_weights(msa)
    if mode == "species_henikoff_cluster":
        labels = np.asarray(species) if species is not None else None
        if (labels is None or labels.ndim != 1 or len(labels) != len(msa)
                or any(label is None for label in labels)
                or (labels.dtype.kind in "fc" and not np.isfinite(labels).all())):
            raise ValueError("Species weighting requires one valid species label per row")
        for label in np.unique(labels):
            rows = labels == label
            total = weights[rows].sum()
            if total > 0:
                weights[rows] /= total
    effective_n = len(msa) if mode == "henikoff" else effective_cluster_count(msa)
    total = weights.sum()
    if total > 0:
        weights *= effective_n / total
    return weights
