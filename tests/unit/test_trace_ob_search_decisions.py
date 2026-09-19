import numpy as np

from benchmark_tools.trace_ob_search_decisions import subset_queries, SETTINGS
from orthohmm.search import engine
from orthohmm.search.sequences import SpeciesSequences


def test_query_subset_preserves_full_target_search(tmp_path, monkeypatch):
    monkeypatch.setattr(engine, "is_cuda_available", lambda: False)
    path = tmp_path / "input.fa"
    sequence = "ACDEFGHIKLMNPQRSTVWY" * 4
    path.write_text(f">first\n{'A' * 80}\n>second\n{sequence}\n>third\n{sequence[::-1]}\n")
    sp = SpeciesSequences.from_fasta(str(path), path.name)
    subset = subset_queries(sp, {"second"})
    assert subset.ids == ["second"]
    assert len(subset.flat_sequences) == sum(subset.lengths)
    np.testing.assert_array_equal(subset.get_sequence(0), sp.get_sequence(1))
    full = engine.search_species_pair_indexed(sp, sp, n_threads=1, **SETTINGS)
    observed = engine.search_species_pair_indexed(subset, sp, n_threads=1, **SETTINGS)
    keep = full.query_indices == 1
    assert np.count_nonzero(keep) > 0
    np.testing.assert_array_equal(observed.target_indices, full.target_indices[keep])
    np.testing.assert_array_equal(observed.scores, full.scores[keep])
    np.testing.assert_array_equal(observed.evalues, full.evalues[keep])


def test_subset_retains_original_order_and_repacks_offsets():
    sp = SpeciesSequences("x", ["c", "a", "b"], np.arange(12, dtype=np.uint8),
                          np.array([0, 3, 7]), np.array([3, 4, 5]))
    subset = subset_queries(sp, {"b", "c"})
    assert subset.ids == ["c", "b"]
    np.testing.assert_array_equal(subset.offsets, [0, 3])
    np.testing.assert_array_equal(subset.lengths, [3, 5])
    np.testing.assert_array_equal(subset.get_sequence(1), sp.get_sequence(2))
