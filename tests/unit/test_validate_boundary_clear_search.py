from types import SimpleNamespace
import json
from pathlib import Path

import numpy as np
import pytest

from benchmark_tools import validate_boundary_clear_search as worker
from orthohmm.search.sequences import SpeciesSequences


def species(name="query"):
    return SpeciesSequences(name,["a","b"],np.array([1,2,3],dtype=np.uint8),
        np.array([0,1],dtype=np.int64),np.array([1,2],dtype=np.int32))


def result(q=(0,1), t=(1,0), scores=(2.,3.), evalues=(.01,.02)):
    return SimpleNamespace(query_species="query",target_species="target", candidate_count=2,
        query_indices=np.array(q,dtype=np.int32),target_indices=np.array(t,dtype=np.int32),
        scores=np.array(scores),evalues=np.array(evalues))


def test_query_reversal_preserves_sequences_and_canonical_outputs():
    original = species()
    reversed_query = worker.reordered_species(original,[1,0])
    assert reversed_query.ids == ["b","a"]
    assert reversed_query.get_sequence(0).tolist() == [2,3]
    assert reversed_query.get_sequence(1).tolist() == [1]
    assert worker.canonical_result(result(),original,species("target")) == worker.canonical_result(
        result(q=(1,0)),reversed_query,species("target"))


@pytest.mark.parametrize("order", [[0,0],[1],[0,2]])
def test_query_order_rejects_nonpermutations(order):
    with pytest.raises(ValueError,match="permutation"):
        worker.reordered_species(species(),order)


@pytest.mark.parametrize("changed", [result(q=(-1,1)),result(t=(1,2)),result(scores=(float("nan"),3.)),
    result(evalues=(-.01,.02)),result(q=(0,0),t=(1,1))])
def test_canonical_result_rejects_invalid_or_duplicate_pairs(changed):
    with pytest.raises(ValueError):
        worker.canonical_result(changed,species(),species("target"))


def test_complete_search_binding_tracks_integer_kernel_and_restores_library_state(monkeypatch):
    prior = worker.viterbi._c_hmm_lib,worker.viterbi._c_hmm_available
    library = object()
    arrays = [np.ones(1,dtype=np.int32) for _ in range(8)] + [np.array([[0,1],[1,0]],dtype=np.int32)]
    def kernel(*args,**kwargs):
        assert worker.viterbi._c_hmm_lib is library
        return np.array([2,3],dtype=np.int32)
    def search(*args,**kwargs):
        assert kwargs["kmer_k"]==4 and kwargs["max_candidates_per_query"]==400
        worker.engine.batch_viterbi_multipair_c(*arrays,64,n_threads=4)
        return result()
    monkeypatch.setattr(worker.engine,"batch_viterbi_multipair_c",kernel)
    monkeypatch.setattr(worker.engine,"search_species_pair_indexed",search)
    actual = worker.search_with_library(species(),species("target"),400,library,capture_inputs=True)
    assert actual["kernel_calls"]==1 and len(actual["kernel_inputs"])==9
    assert (worker.viterbi._c_hmm_lib,worker.viterbi._c_hmm_available)==prior


def test_complete_search_rejects_silent_fallback_and_restores_library_state(monkeypatch):
    prior = worker.viterbi._c_hmm_lib,worker.viterbi._c_hmm_available
    def search(*args,**kwargs):
        worker.engine.batch_viterbi_c()
    monkeypatch.setattr(worker.engine,"search_species_pair_indexed",search)
    with pytest.raises(AssertionError,match="fell back"):
        worker.search_with_library(species(),species("target"),100,object())
    assert (worker.viterbi._c_hmm_lib,worker.viterbi._c_hmm_available)==prior


def test_complete_search_requires_real_kernel_invocation(monkeypatch):
    monkeypatch.setattr(worker.engine,"search_species_pair_indexed",lambda *args,**kwargs:result())
    with pytest.raises(ValueError,match="actual multipair"):
        worker.search_with_library(species(),species("target"),100,object())


def test_warm_capture_rejects_kernel_input_mutation(monkeypatch):
    arrays = [np.ones(1,dtype=np.int32) for _ in range(8)] + [np.array([[0,1],[1,0]],dtype=np.int32)]
    def kernel(*args,**kwargs):
        args[0][0] += 1
        return np.array([2,3],dtype=np.int32)
    def search(*args,**kwargs):
        worker.engine.batch_viterbi_multipair_c(*arrays,64,n_threads=4)
        return result()
    monkeypatch.setattr(worker.engine,"batch_viterbi_multipair_c",kernel)
    monkeypatch.setattr(worker.engine,"search_species_pair_indexed",search)
    with pytest.raises(ValueError,match="modified search input"):
        worker.search_with_library(species(),species("target"),100,object(),capture_inputs=True)


def test_all_search_settings_match_prespecified_protocol():
    protocol = json.loads(Path(worker.__file__).with_name(worker.PROTOCOL).read_text())
    assert protocol["settings"] == worker.SETTINGS
