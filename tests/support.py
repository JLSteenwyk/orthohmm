"""Small deterministic fixtures for application tests, without research drivers."""

import random

import numpy as np

from orthohmm.search.sequences import SpeciesSequences, encode_sequence


def fixture():
    rng = random.Random(20260926)
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    ancestors = ["".join(rng.choice(alphabet) for _ in range(200)) for _ in range(3)]
    result = {}
    for taxon in ("S1", "S2", "S3", "S4"):
        genes = {}
        for family, ancestor in enumerate(ancestors):
            sequence = list(ancestor)
            for position in rng.sample(range(200), 10):
                sequence[position] = rng.choice(alphabet.replace(sequence[position], ""))
            sequence = "".join(sequence)
            genes[f"{taxon}_family{family}"] = sequence
            if family == 2:
                genes[f"{taxon}_family{family}_duplicate"] = sequence
        result[taxon] = genes
    return result


def species_sequences(ids, sequences, name):
    encoded = [encode_sequence(sequence) for sequence in sequences]
    lengths = np.asarray([len(sequence) for sequence in encoded], dtype=np.int32)
    offsets = np.zeros(len(encoded), dtype=np.int64)
    if len(encoded) > 1:
        offsets[1:] = np.cumsum(lengths[:-1], dtype=np.int64)
    flat = (np.concatenate(encoded).astype(np.uint8, copy=False)
            if encoded else np.empty(0, dtype=np.uint8))
    return SpeciesSequences(name, list(ids), flat, offsets, lengths)
