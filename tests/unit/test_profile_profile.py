import numpy as np

from orthohmm.search.matrices import get_background_freqs, get_matrix
from orthohmm.search.msa_profile import build_msa_profile
from orthohmm.search.profile_profile import build_profile_profile_edges


def _profile(sequences, prefix):
    return build_msa_profile(
        sequences,
        [f"{prefix}_{index}" for index in range(len(sequences))],
        get_matrix("BLOSUM62"),
        get_background_freqs("BLOSUM62"),
    )


def test_profile_profile_edges_join_mutually_closest_split_family():
    profiles = {
        0: _profile([
            "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY",
            "ACDEFGHIKLMNPQRSTVWYACDEYGHIKLMNPQRSTVWY",
            "ACDEYGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY",
        ], "left"),
        1: _profile([
            "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQKSTVWY",
            "ACDEFGHIKLMNPQRSTVWYACDEYGHIKLMNPQKSTVWY",
            "ACDEYGHIKLMNPQRSTVWYACDEFGHIKLMNPQKSTVWY",
        ], "right"),
        2: _profile([
            "WYVTSRQPNMLKIHGFEDCAWYVTSRQPNMLKIHGFEDCA",
            "WYVTSRQPNMLKIHGFEDCAWYVTSRQPNMLKIHGYEDCA",
            "WYVTSRQPNMLKIHGFEDCAWYVTSRQPNMLKIHAFEDCA",
        ], "other"),
    }
    clusters = [[0, 1, 2], [3, 4, 5], [6, 7, 8]]

    result = build_profile_profile_edges(
        profiles,
        clusters,
        [f"gene_{index}" for index in range(9)],
        get_background_freqs("BLOSUM62"),
        cpu=1,
        similarity_threshold=0.6,
    )

    assert result.candidate_pairs >= 1
    assert result.reciprocal_pairs == 1
    observed = {
        (int(source), int(target))
        for source, target in zip(result.edges.sources, result.edges.targets)
    }
    assert observed == {
        (source, target) for source in range(3) for target in range(3, 6)
    }


def test_profile_profile_edges_validate_bounds():
    with np.testing.assert_raises_regex(ValueError, "similarity_threshold"):
        build_profile_profile_edges({}, [], [], np.full(20, 0.05), 1, 0.0)
