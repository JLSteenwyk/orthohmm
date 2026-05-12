"""Unit tests for orthohmm/search/profile.py (profile HMM construction)."""
import numpy as np
import pytest

from orthohmm.search.matrices import (
    ALPHABET, ALPHABET_SIZE, get_background_freqs, get_matrix,
)
from orthohmm.search.profile import (
    ProfileHMM,
    build_profile,
    build_profiles_batch,
)
from orthohmm.search.sequences import SpeciesSequences, encode_sequence


# ─── ProfileHMM dataclass ───────────────────────────────────────────────


class TestProfileHMMDataclass:
    def test_holds_arrays_with_expected_shapes(self):
        p = ProfileHMM(
            length=10,
            match_emissions=np.zeros((10, ALPHABET_SIZE), dtype=np.int8),
            insert_emissions=np.zeros(ALPHABET_SIZE, dtype=np.int8),
            transitions=np.zeros(7, dtype=np.int32),
        )
        assert p.length == 10
        assert p.match_emissions.shape == (10, 20)
        assert p.insert_emissions.shape == (20,)
        assert p.transitions.shape == (7,)


# ─── build_profile ──────────────────────────────────────────────────────


class TestBuildProfile:
    def _setup(self):
        sub = get_matrix("BLOSUM62")
        bg = get_background_freqs("BLOSUM62")
        return sub, bg

    def test_match_emissions_copy_substitution_rows(self):
        sub, bg = self._setup()
        seq_str = "ACDE"
        encoded = encode_sequence(seq_str)
        prof = build_profile(encoded, len(seq_str), sub, bg)
        assert prof.length == 4
        # Each position's match emission row is the substitution matrix row
        # for the residue at that position
        for i, c in enumerate(seq_str):
            aa = ALPHABET.index(c)
            assert np.array_equal(prof.match_emissions[i], sub[aa])

    def test_unknown_residue_emits_zeros(self):
        sub, bg = self._setup()
        # Mix in an 'X' (unknown) at position 1
        encoded = encode_sequence("AXC")
        prof = build_profile(encoded, 3, sub, bg)
        # Position 0 = A, position 2 = C → substitution rows
        assert np.array_equal(prof.match_emissions[0], sub[ALPHABET.index("A")])
        assert np.array_equal(prof.match_emissions[2], sub[ALPHABET.index("C")])
        # Position 1 = unknown → all zeros
        assert np.all(prof.match_emissions[1] == 0)

    def test_transitions_have_expected_layout(self):
        sub, bg = self._setup()
        prof = build_profile(encode_sequence("AAAA"), 4, sub, bg,
                             gap_open=-12, gap_extend=-3)
        # transitions = [T_MM, T_MI, T_MD, T_IM, T_II, T_DM, T_DD]
        t = prof.transitions
        assert t.shape == (7,)
        assert t[0] == 0          # T_MM
        assert t[1] == -12        # T_MI = gap_open
        assert t[2] == -12        # T_MD = gap_open
        assert t[4] == -3         # T_II = gap_extend
        assert t[6] == -3         # T_DD = gap_extend

    def test_insert_emissions_have_uniform_penalty(self):
        sub, bg = self._setup()
        prof = build_profile(encode_sequence("ACDE"), 4, sub, bg)
        # All 20 amino acids get the same insert emission score
        assert prof.insert_emissions.shape == (20,)
        assert (prof.insert_emissions == prof.insert_emissions[0]).all()


# ─── build_profiles_batch ───────────────────────────────────────────────


class TestBuildProfilesBatch:
    def _setup(self, tmp_path, records):
        p = tmp_path / "x.fa"
        p.write_text("".join(f">{h}\n{s}\n" for h, s in records))
        return SpeciesSequences.from_fasta(str(p), "x.fa")

    def test_flat_layout_matches_individual_profiles(self, tmp_path):
        sp = self._setup(tmp_path, [("g1", "ACDE"), ("g2", "FGHI")])
        sub = get_matrix("BLOSUM62")
        bg = get_background_freqs("BLOSUM62")

        flat_match, ins_emit, trans, offsets, lens = build_profiles_batch(
            sp, sub, bg,
        )
        # Total residues = 4 + 4 = 8
        assert flat_match.shape == (8, 20)
        assert list(lens) == [4, 4]
        assert list(offsets) == [0, 4]
        # Verify g1's first position = substitution row for 'A'
        a = ALPHABET.index("A")
        assert np.array_equal(flat_match[0], sub[a])
        # g2 starts at offset 4; first position = 'F'
        f = ALPHABET.index("F")
        assert np.array_equal(flat_match[4], sub[f])

    def test_transitions_and_insert_shared_across_profiles(self, tmp_path):
        sp = self._setup(tmp_path, [("g1", "AAA")])
        sub = get_matrix("BLOSUM62")
        bg = get_background_freqs("BLOSUM62")
        _, ins, trans, _, _ = build_profiles_batch(
            sp, sub, bg, gap_open=-15, gap_extend=-2,
        )
        assert trans.shape == (7,)
        assert ins.shape == (20,)
        assert trans[1] == -15   # T_MI
        assert trans[4] == -2    # T_II
