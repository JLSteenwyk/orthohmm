from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest

from orthohmm.search.matrices import get_background_freqs, get_matrix
from orthohmm.search.msa_center_star import center_star_msa
from orthohmm.search.msa_profile import build_msa_profile


BASE = "ACDEFGHIKLMNPQRSTVWY" * 3


def align(sequences, **options):
    return center_star_msa(sequences, get_matrix("BLOSUM62"),
                           gap_open=-12, gap_extend=-3, n_threads=1, **options)


def assert_preserved(sequences, aligned):
    assert aligned is not None and len(aligned) == len(sequences)
    assert len({len(row) for row in aligned}) == 1
    assert [row.replace("-", "") for row in aligned] == sequences


@pytest.mark.parametrize("position", [0, 1, 29, 30, 59, 60])
@pytest.mark.parametrize("size", [1, 3, 8])
def test_every_source_residue_survives_insertions(position, size):
    sequences = [BASE, BASE, BASE[:position] + "G" * size + BASE[position:]]
    assert_preserved(sequences, align(sequences, projection_policy="residue_preserving"))


@pytest.mark.parametrize("position", [0, 30, 60])
def test_default_legacy_projection_is_unchanged_even_when_lossy(position):
    sequences = [BASE, BASE, BASE[:position] + "GGG" + BASE[position:]]
    default = align(sequences)
    assert default == align(sequences, projection_policy="legacy")
    assert len(default[-1].replace("-", "")) == 60
    assert len(sequences[-1]) == 63
    assert default != align(sequences, projection_policy="residue_preserving")


@pytest.mark.parametrize("sequences", [
    [BASE, BASE, "XXX" + BASE + "XX"],
    [BASE[:20] + "XXX" + BASE[20:], BASE[:20] + "X" + BASE[20:],
     BASE[:20] + "XXXXXX" + BASE[20:]],
    ["X" * 20, "X" * 25, "X" * 30],
    [BASE, BASE[:15] + "GGG" + BASE[15:45] + "PPP" + BASE[45:],
     BASE[:10] + "XX" + BASE[10:40] + BASE[45:]],
])
def test_unknowns_and_multiple_indels_do_not_alias_traceback_gaps(sequences):
    assert_preserved(sequences, align(sequences, projection_policy="residue_preserving"))


@pytest.mark.parametrize("seed", range(8))
def test_randomized_projection_preserves_rows_and_threaded_calls(seed):
    rng = np.random.default_rng(seed)
    alphabet = np.array(list("ACDEFGHIKLMNPQRSTVWYX"))
    ancestor = list(rng.choice(alphabet, 35))
    families = []
    for _ in range(3):
        sequences = []
        for _ in range(5):
            sequence = ancestor.copy()
            for _ in range(4):
                position = int(rng.integers(len(sequence) + 1))
                sequence[position:position] = list(rng.choice(alphabet, int(rng.integers(1, 5))))
            del sequence[5:8]
            sequences.append("".join(sequence))
        families.append(sequences)
    def run(sequences):
        return align(sequences, projection_policy="residue_preserving")
    serial = [run(sequences) for sequences in families]
    with ThreadPoolExecutor(max_workers=3) as executor:
        threaded = list(executor.map(run, families))
    assert threaded == serial
    for sequences, aligned in zip(families, serial):
        assert_preserved(sequences, aligned)


def test_profile_records_opt_in_policy_without_changing_legacy_default():
    sequences = [BASE, BASE, "GGG" + BASE]
    args = (sequences, ["a", "b", "c"], get_matrix("BLOSUM62"), get_background_freqs("BLOSUM62"))
    default = build_msa_profile(*args)
    legacy = build_msa_profile(*args, projection_policy="legacy")
    fixed = build_msa_profile(*args, projection_policy="residue_preserving")
    assert default.projection_policy == legacy.projection_policy == "legacy"
    assert fixed.projection_policy == "residue_preserving"
    assert default.member_ids == fixed.member_ids == ["a", "b", "c"]
    for name in ("match_emissions", "insert_emissions", "transitions", "consensus"):
        assert np.array_equal(getattr(default, name), getattr(legacy, name))


def test_unknown_projection_policy_is_rejected_before_construction():
    with pytest.raises(ValueError, match="projection policy"):
        align([BASE], projection_policy="unregistered")
    with pytest.raises(ValueError, match="projection policy"):
        build_msa_profile([BASE], ["a"], get_matrix("BLOSUM62"),
                          get_background_freqs("BLOSUM62"), projection_policy="unregistered")
