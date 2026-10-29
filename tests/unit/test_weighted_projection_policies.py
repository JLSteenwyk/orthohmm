from dataclasses import replace
import json

import numpy as np
import pytest

from orthohmm.args_processing import process_args
from orthohmm.parser import create_parser
from orthohmm.search import msa_profile
from orthohmm.search import profile_expansion as expansion
from orthohmm.search.matrices import get_background_freqs, get_matrix
from orthohmm.search.msa_center_star import (
    WEIGHTED_PROJECTION_MODES, center_star_msa, validate_profile_projection_combination,
)


def constituent_digest(builder):
    import hashlib
    import numpy as np
    from orthohmm.search.matrices import get_background_freqs, get_matrix

    rng = np.random.default_rng(20261019)
    digest = hashlib.sha256()
    matrix, background = get_matrix("BLOSUM62"), get_background_freqs("BLOSUM62")
    alphabet = np.array(list("ACDEFGHIKLMNPQRSTVWY"))
    cases = (("none", "legacy"), ("henikoff", "legacy"),
             ("henikoff_cluster", "legacy"), ("none", "residue_preserving"))
    for index in range(96):
        base = "".join(rng.choice(alphabet, size=int(rng.integers(12, 48))))
        sequences = [base, base, "GGX" + base, base[:len(base) // 2] + "CC" + base[len(base) // 2:]]
        if index % 2: sequences[1] = base[2:]
        if index % 3 == 0: sequences[2] = sequences[2][:-3]
        for mode, policy in cases:
            profile = builder(sequences=sequences, ids=["a", "b", "c", "d"], sub_matrix=matrix, bg_freqs=background,
                              sequence_weighting=mode, projection_policy=policy)
            assert profile is not None and profile.member_ids == ["a", "b", "c", "d"] and profile.length >= 6
            digest.update(mode.encode() + b":" + policy.encode())
            digest.update(int(profile.length).to_bytes(8, "little"))
            for field, dtype in (("match_emissions", "i1"), ("insert_emissions", "i1"),
                                 ("transitions", "<i4"), ("consensus", "u1")):
                digest.update(getattr(profile, field).astype(dtype).tobytes())
            digest.update(np.array([profile.effective_sequence_count], dtype="<f8").tobytes())
    return digest.hexdigest()


def test_384_constituent_profiles_match_frozen_precombination_source():
    assert constituent_digest(msa_profile.build_msa_profile) == "aa66538d3a2561d3f8899a3b67346b3f502c77a8626a641771872199df146743"




@pytest.mark.parametrize("mode", ("none", *WEIGHTED_PROJECTION_MODES))
def test_whitelist_supports_only_registered_joint_constructors(mode):
    assert validate_profile_projection_combination("residue_preserving", mode) == "residue_preserving"
    with pytest.raises(ValueError, match="separate experiment"):
        validate_profile_projection_combination("residue_preserving", mode, "strict_span_guard")
    with pytest.raises(ValueError, match="separate experiment"):
        validate_profile_projection_combination("residue_preserving", "species_henikoff_cluster")


@pytest.mark.parametrize("mode,weights,count", [
    ("henikoff", [9 / 14, 9 / 14, 12 / 7], 3),
    ("henikoff_cluster", [3 / 7, 3 / 7, 8 / 7], 2),
])
def test_repaired_column_weights_and_unknowns_match_hand_calculation(monkeypatch, mode, weights, count):
    aligned = ["AA-X", "AA--", "CCGX"]
    calls = []
    def align(*args, **kwargs):
        assert args[0] == ["AAX", "AA", "CCGX"]
        calls.append(kwargs["projection_policy"])
        return aligned
    monkeypatch.setattr(msa_profile, "center_star_msa", align)
    matrix, background = get_matrix("BLOSUM62"), get_background_freqs("BLOSUM62")
    profile = msa_profile.build_msa_profile(sequences=["AAX", "AA", "CCGX"], ids=["a", "b", "c"],
        sub_matrix=matrix, bg_freqs=background,
        sequence_weighting=mode, projection_policy="residue_preserving")
    encoded = msa_profile.encode_alignment(aligned)
    actual = msa_profile.profile_sequence_weights(encoded, mode)
    np.testing.assert_allclose(actual, weights)
    expected, mask, _ = msa_profile.compute_pssm(encoded, matrix, background, sequence_weights=weights)
    np.testing.assert_array_equal(profile.match_emissions, expected)
    assert mask.tolist() == [True, True, False, False]
    assert profile.effective_sequence_count == pytest.approx(count)
    assert profile.member_ids == ["a", "b", "c"]
    assert profile.projection_policy == calls[0] == "residue_preserving" and profile.sequence_weighting == mode


@pytest.mark.parametrize("mode", WEIGHTED_PROJECTION_MODES)
def test_real_joint_expansion_records_actual_nonempty_profiles_and_recruitment(tmp_path, mode):
    names = ["a0", "a1", "a2", "full", "fragment", "partial"]
    seqs = ["A" * 100, "A" * 100, "A" * 20, "A" * 100, "A" * 20, "C" * 50 + "A" * 20 + "C" * 50]
    (tmp_path / "input.fa").write_text("".join(f">{name}\n{seq}\n" for name, seq in zip(names, seqs)))
    result = expansion.expand_profiles([[0, 1, 2], [3], [4], [5]], names, str(tmp_path), ["input.fa"],
        "BLOSUM62", 1, 1e-4, np.array([3, 4, 5]), np.array([0, 0, 0]), np.array([.9, .8, .7]),
        sequence_weighting=mode, projection_policy="residue_preserving")
    assert result.profiles_built == 1 and result.profile_candidates > 0
    assert (0, 3, .9) in list(zip(result.edges.sources.tolist(), result.edges.targets.tolist(), result.edges.weights.tolist()))
    records = result.profile_construction_policies
    assert len(records) == 1 and records[0]["cluster_id"] == 0
    assert records[0]["sequence_weighting"] == mode and records[0]["projection_policy"] == "residue_preserving"
    assert records[0]["effective_sequence_count"] == pytest.approx(3 if mode == "henikoff" else 1)
    assert records[0]["effective_sequence_count"] == result.profile_effective_counts[0][1]
    json.dumps(records, allow_nan=False)


@pytest.mark.parametrize("change", [{"sequence_weighting": "none"}, {"projection_policy": "legacy"},
                                  {"effective_sequence_count": float("nan")}, {"effective_sequence_count": 0.}])
def test_silently_ignored_policy_or_invalid_effective_count_fails_before_scoring(tmp_path, monkeypatch, change):
    names = ["a", "b", "c"]
    base = "ACDEFGHIKLMNPQRSTVWY" * 3
    (tmp_path / "input.fa").write_text("".join(f">{name}\n{base}\n" for name in names))
    model = msa_profile.build_msa_profile(sequences=[base] * 3, ids=names,
        sub_matrix=get_matrix("BLOSUM62"), bg_freqs=get_background_freqs("BLOSUM62"),
        sequence_weighting="henikoff", projection_policy="residue_preserving")
    monkeypatch.setattr(expansion, "build_cluster_profiles", lambda *a, **kw: {0: replace(model, **change)})
    def forbidden(*args, **kwargs): raise AssertionError("Malformed joint profiles reached scoring")
    monkeypatch.setattr(expansion, "search_genes_against_profiles", forbidden)
    with pytest.raises(ValueError, match="both requested policies"):
        expansion.expand_profiles([[0, 1, 2]], names, str(tmp_path), ["input.fa"], "BLOSUM62", 1, 1e-4,
            [], [], [], sequence_weighting="henikoff", projection_policy="residue_preserving")


def test_unknown_residues_survive_real_alignment_before_existing_missing_value_weight_encoding():
    seqs = ["ACDEFGHIKLMNPQRSTVWY", "ACDEXFGHIKLMNPQRSTVWY", "GGGACDEFGHIKLMNPQRSTVWY"]
    aligned = center_star_msa(seqs, get_matrix("BLOSUM62"), projection_policy="residue_preserving")
    assert [row.replace("-", "") for row in aligned] == seqs
    encoded = msa_profile.encode_alignment(aligned)
    assert encoded.min() >= 0 and encoded.max() <= 20
    for mode in WEIGHTED_PROJECTION_MODES:
        assert np.isfinite(msa_profile.profile_sequence_weights(encoded, mode)).all()




@pytest.mark.parametrize("weighting", WEIGHTED_PROJECTION_MODES)
def test_joint_worker_construction_errors_propagate_instead_of_becoming_missing_profiles(monkeypatch, weighting):
    def fail(*args, **kwargs): raise ValueError("fixture construction failure")
    monkeypatch.setattr(expansion, "build_msa_profile", fail)
    task = (7, ["a", "b"], [np.array([0, 1], dtype=np.uint8)] * 2,
            get_matrix("BLOSUM62"), get_background_freqs("BLOSUM62"),
            {"sequence_weighting": weighting, "projection_policy": "residue_preserving"})
    with pytest.raises(RuntimeError, match="cluster 7: fixture construction failure"):
        expansion._build_profile_worker(task)
