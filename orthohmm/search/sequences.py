"""Compact sequence storage for the built-in search engine.

Stores protein sequences as integer-encoded flat arrays with offset
indexing, adapted from ClustKIT's SequenceDataset pattern.
"""

import os
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

import numpy as np

from .matrices import ALPHABET_SIZE, FULL_CHAR_TO_IDX


FASTA_EXTENSIONS = (".fa", ".faa", ".fas", ".fasta", ".pep", ".prot")


def encode_sequence(seq_str: str) -> np.ndarray:
    """Encode a protein sequence string to uint8 array.

    Uses ACDEFGHIKLMNPQRSTVWY alphabet (0-19).
    Unknown residues map to 20.
    """
    encoded = np.empty(len(seq_str), dtype=np.uint8)
    for i, c in enumerate(seq_str):
        encoded[i] = FULL_CHAR_TO_IDX.get(ord(c), ALPHABET_SIZE)
    return encoded


@dataclass
class SpeciesSequences:
    """All sequences for one species FASTA file, stored compactly."""
    species_file: str
    ids: List[str]
    flat_sequences: np.ndarray  # uint8, concatenated integer-encoded
    offsets: np.ndarray         # int64, start offset per sequence
    lengths: np.ndarray         # int32, length per sequence

    @property
    def num_sequences(self) -> int:
        return len(self.ids)

    def get_sequence(self, idx: int) -> np.ndarray:
        """Return the encoded sequence at index idx."""
        start = self.offsets[idx]
        length = self.lengths[idx]
        return self.flat_sequences[start:start + length]

    @classmethod
    def from_fasta(cls, fasta_path: str, species_file: str) -> "SpeciesSequences":
        """Parse a FASTA file into compact storage."""
        ids = []
        sequences = []
        seen_ids = set()

        header = None
        seq_parts = []

        with open(fasta_path, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if header is not None:
                        seq_str = "".join(seq_parts)
                        sequences.append(encode_sequence(seq_str))
                    header = line[1:].split()[0]
                    if header in seen_ids:
                        print(
                            f"{header} appears twice in file {species_file}. "
                            "Please ensure all FASTA headers are unique within each file."
                        )
                        sys.exit(1)
                    seen_ids.add(header)
                    ids.append(header)
                    seq_parts = []
                else:
                    seq_parts.append(line)

            if header is not None:
                seq_str = "".join(seq_parts)
                sequences.append(encode_sequence(seq_str))

        if not sequences:
            return cls(
                species_file=species_file,
                ids=[],
                flat_sequences=np.empty(0, dtype=np.uint8),
                offsets=np.empty(0, dtype=np.int64),
                lengths=np.empty(0, dtype=np.int32),
            )

        lengths = np.array([len(s) for s in sequences], dtype=np.int32)
        offsets = np.zeros(len(sequences), dtype=np.int64)
        offsets[1:] = np.cumsum(lengths[:-1])
        flat_sequences = np.concatenate(sequences)

        return cls(
            species_file=species_file,
            ids=ids,
            flat_sequences=flat_sequences,
            offsets=offsets,
            lengths=lengths,
        )


@dataclass
class SequenceStore:
    """All species loaded into memory for search."""
    species: Dict[str, SpeciesSequences] = field(default_factory=dict)

    @classmethod
    def from_fasta_directory(
        cls, fasta_directory: str, files: List[str]
    ) -> "SequenceStore":
        """Load all FASTA files from a directory."""
        store = cls()
        for filename in files:
            path = os.path.join(fasta_directory, filename)
            store.species[filename] = SpeciesSequences.from_fasta(path, filename)
        return store

    def get_gene_lengths(self) -> np.ndarray:
        """Return structured array matching OrthoHMM's get_sequence_lengths output.

        dtype = [("spp", object), ("name", object), ("length", int)]
        """
        entries = []
        for filename, sp in self.species.items():
            for i, seq_id in enumerate(sp.ids):
                entries.append((filename, seq_id, int(sp.lengths[i])))

        dtype = [("spp", object), ("name", object), ("length", int)]
        return np.array(entries, dtype=dtype)
