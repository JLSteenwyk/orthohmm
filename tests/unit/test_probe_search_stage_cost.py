import numpy as np
import pytest

from benchmark_tools import probe_search_stage_cost as probe
from orthohmm.search.sequences import SpeciesSequences


def species():
    return SpeciesSequences("A.faa", ["z", "a", "h", "b"], np.arange(14, dtype=np.uint8),
                            np.array([0, 3, 8, 10], dtype=np.int64), np.array([3, 5, 2, 4], dtype=np.int32))


def test_hash_sample_is_order_independent_and_subset_restores_original_gene_order():
    original = species()
    selected = probe.selected_names(original.ids, 3)
    assert selected == probe.selected_names(list(reversed(original.ids)), 3)
    subset, indices = probe.subset_species(original, 3)
    assert indices == sorted(indices) and subset.ids == [original.ids[i] for i in indices]
    assert set(subset.ids) == set(selected) and subset.species_file == original.species_file
    assert subset.offsets.dtype == np.int64 and subset.lengths.dtype == np.int32
    for new, old in enumerate(indices):
        assert np.array_equal(subset.get_sequence(new), original.get_sequence(old))
    assert np.array_equal(original.flat_sequences, np.arange(14, dtype=np.uint8))


@pytest.mark.parametrize("names,count", [(["a", "a"], 1), (["a"], 2), (["a"], 0), (["a"], True)])
def test_invalid_sampling_is_not_silently_changed(names, count):
    with pytest.raises(ValueError):
        probe.selected_names(names, count)


def test_timing_wrapper_preserves_argument_identity_return_identity_and_exception():
    value, result = object(), object()
    measurements = {}
    def operation(argument, *, setting):
        assert argument is value and setting == 3
        return result
    wrapped = probe.timed_wrapper(operation, measurements, "stage")
    assert wrapped(value, setting=3) is result
    assert measurements["stage"]["calls"] == 1 and measurements["stage"]["wall_s"] >= 0
    def fails():
        raise RuntimeError("native failure")
    with pytest.raises(RuntimeError, match="native failure"):
        probe.timed_wrapper(fails, measurements, "failure")()
    assert measurements["failure"]["calls"] == 1


def result():
    return probe.engine.IndexedSpeciesPairResults("A.faa", "B.faa", np.array([1], dtype=np.int32),
        np.array([0], dtype=np.int32), np.array([1e-6]), np.array([2.]), 1)


def test_search_instrumentation_preserves_native_settings_and_restores_functions(monkeypatch):
    query, target = species(), species()
    expected = result()
    before = probe.engine.build_profiles_batch
    def search(q, t, matrix, **kwargs):
        assert q is query and t is target and matrix == "BLOSUM62"
        assert kwargs == dict(band_width=64, kmer_k=3, max_candidates_per_query=100,
                              n_threads=4, prefilter_phase_a_topk=5000)
        return expected
    monkeypatch.setattr(probe.engine, "search_species_pair_indexed", search)
    a = probe.search_pair(query, target, 100, 4, False)
    b = probe.search_pair(query, target, 100, 4, True)
    assert a["fingerprint"] == b["fingerprint"] and b["significant_hits"] == 1
    assert probe.engine.build_profiles_batch is before


def test_native_array_changes_cannot_pass_fingerprint_equivalence():
    a, b = result(), result()
    assert probe.fingerprint_result(a) == probe.fingerprint_result(b)
    b.scores[0] += 1e-6
    assert probe.fingerprint_result(a) != probe.fingerprint_result(b)


@pytest.mark.parametrize("mutation", ["none", "python", "binary", "inventory"])
def test_probe_requires_exact_live_search_source_and_cpu_extensions(tmp_path, mutation):
    old, new = tmp_path / "old", tmp_path / "new"
    for root in (old, new):
        folder = root / "orthohmm/search/csrc"
        folder.mkdir(parents=True)
        (folder.parent / "engine.py").write_text("search source\n")
        for name in ("hmm_viterbi.so", "kmer_prefilter.so", "pair_align.so"):
            (folder / name).write_text("identity fixture\n")
    if mutation == "python":
        (new / "orthohmm/search/engine.py").write_text("changed\n")
    elif mutation == "binary":
        (new / "orthohmm/search/csrc/hmm_viterbi.so").write_text("changed\n")
    elif mutation == "inventory":
        (new / "orthohmm/search/extra.py").write_text("new\n")
    if mutation == "none":
        assert len(probe.search_source_records(new, old)) == 8
    else:
        with pytest.raises(ValueError):
            probe.search_source_records(new, old)
