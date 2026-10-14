import numpy as np
import pytest

from benchmark_tools.audit_hit_arrival_order import block_layout,block_permutation,toy_reproducer


def test_toy_confirms_unchanged_exact_tie_depends_on_input_order():
    result = toy_reproducer()
    assert result["original_assignment"]==1 and result["reversed_assignment"]==2


def test_block_permutations_preserve_all_rows_and_internal_order():
    species = np.array([0,0,1,1],dtype=np.int32)
    q,t = np.array([2,3,0,1]),np.array([0,1,2,3])
    starts,ends,keys = block_layout(species,q,t)
    assert keys.tolist()==[2,1]
    forward = block_permutation(starts,ends,keys)
    backward = block_permutation(starts,ends,keys,True)
    assert forward.tolist()==[2,3,0,1]
    assert backward.tolist()==[0,1,2,3]
    assert sorted(forward)==list(range(4))


def test_block_layout_rejects_interleaved_species_pairs():
    with pytest.raises(ValueError,match="contiguous"):
        block_layout(np.array([0,1]),np.array([0,1,0]),np.array([1,0,1]))


def test_block_layout_rejects_empty_hits():
    with pytest.raises(ValueError,match="nonempty"):
        block_layout(np.array([0,1]),np.array([],dtype=int),np.array([],dtype=int))
