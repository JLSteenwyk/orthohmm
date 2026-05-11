"""Unit tests for orthohmm/search/sequences.py."""
import numpy as np
import pytest

from orthohmm.search.matrices import ALPHABET, ALPHABET_SIZE
from orthohmm.search.sequences import (
    SequenceStore,
    SpeciesSequences,
    encode_sequence,
)


class TestEncodeSequence:
    def test_alphabet_in_order(self):
        enc = encode_sequence(ALPHABET)  # "ACDEFGHIKLMNPQRSTVWY"
        assert enc.dtype == np.uint8
        assert list(enc) == list(range(ALPHABET_SIZE))

    def test_unknown_residues_get_sentinel(self):
        # X is not in the canonical 20; should map to ALPHABET_SIZE (20)
        enc = encode_sequence("AXC")
        assert list(enc) == [0, ALPHABET_SIZE, 1]

    def test_empty_string(self):
        enc = encode_sequence("")
        assert enc.shape == (0,)
        assert enc.dtype == np.uint8


class TestSpeciesSequencesFromFasta:
    def test_two_record_file(self, tmp_path):
        p = tmp_path / "x.fa"
        p.write_text(
            ">geneA\nACDE\nFGHI\n"     # 8 residues split across two lines
            ">geneB extra info\nKLMN\n"  # description after first whitespace is dropped
        )
        sp = SpeciesSequences.from_fasta(str(p), "x.fa")
        assert sp.species_file == "x.fa"
        assert sp.ids == ["geneA", "geneB"]
        assert sp.num_sequences == 2
        assert list(sp.lengths) == [8, 4]
        # offsets: 0, then cumulative
        assert list(sp.offsets) == [0, 8]
        # geneA is ACDEFGHI; encoded indices: 0,1,2,3,4,5,6,7
        assert list(sp.get_sequence(0)) == [0, 1, 2, 3, 4, 5, 6, 7]
        # geneB is KLMN -> K=8, L=9, M=10, N=11
        assert list(sp.get_sequence(1)) == [8, 9, 10, 11]

    def test_empty_fasta_returns_zero_length_arrays(self, tmp_path):
        p = tmp_path / "empty.fa"
        p.write_text("")
        sp = SpeciesSequences.from_fasta(str(p), "empty.fa")
        assert sp.ids == []
        assert sp.num_sequences == 0
        assert sp.flat_sequences.shape == (0,)
        assert sp.offsets.shape == (0,)
        assert sp.lengths.shape == (0,)

    def test_duplicate_header_within_species_exits(self, tmp_path):
        p = tmp_path / "dup.fa"
        p.write_text(">a\nACDE\n>a\nACDE\n")
        with pytest.raises(SystemExit):
            SpeciesSequences.from_fasta(str(p), "dup.fa")


class TestSequenceStoreFromDirectory:
    def _write(self, tmp_path, name, records):
        p = tmp_path / name
        p.write_text("".join(f">{h}\n{s}\n" for h, s in records))
        return p

    def test_loads_multiple_species(self, tmp_path):
        self._write(tmp_path, "h.fa", [("h1", "ACDE")])
        self._write(tmp_path, "m.fa", [("m1", "FGHI"), ("m2", "KLMN")])
        store = SequenceStore.from_fasta_directory(
            str(tmp_path), ["h.fa", "m.fa"],
        )
        assert sorted(store.species.keys()) == ["h.fa", "m.fa"]
        assert store.species["h.fa"].num_sequences == 1
        assert store.species["m.fa"].num_sequences == 2

    def test_get_gene_lengths_structured_array_shape(self, tmp_path):
        self._write(tmp_path, "h.fa", [("h1", "ACDE")])
        self._write(tmp_path, "m.fa", [("m1", "FGHI"), ("m2", "KLMN")])
        store = SequenceStore.from_fasta_directory(
            str(tmp_path), ["h.fa", "m.fa"],
        )
        gl = store.get_gene_lengths()
        # 3 total entries (1 + 2)
        assert gl.shape == (3,)
        # dtype matches OrthoHMM's get_sequence_lengths convention
        assert gl.dtype.names == ("spp", "name", "length")
        names = sorted(gl["name"].tolist())
        assert names == ["h1", "m1", "m2"]
        # Each row reports its source species + length
        lens = {row["name"]: row["length"] for row in gl}
        assert lens == {"h1": 4, "m1": 4, "m2": 4}
