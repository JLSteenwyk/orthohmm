import numpy as np
import pytest

from benchmark_tools.audit_ob_canonical_candidates import independently_sorted,identities
from benchmark_tools.candidate_hit_order_policy import canonical_hit_order_v1
from benchmark_tools.probe_ob_canonical_candidates import array_identity


def test_independent_sort_matches_directed_keys_without_changing_values():
    q,t,s=np.array([2,0,1,0]),np.array([0,1,0,0]),np.array([.1,.2,.3,4.])
    expected=canonical_hit_order_v1(list("abc"),q,t,s)
    observed=independently_sorted(3,q,t,s)
    assert all(np.array_equal(a,b) for a,b in zip(expected,observed))
    assert identities(observed)==[array_identity(a) for a in expected]


@pytest.mark.parametrize("n",[0,3037000500])
def test_integer_encoding_bounds(n):
    with pytest.raises(ValueError):
        independently_sorted(n,np.array([0]),np.array([0]),np.array([1.]))


@pytest.mark.parametrize("problem",["duplicate","negative","outside","float_index","nan","length","shape"])
def test_corrupt_numeric_inputs(problem):
    q,t,s=np.array([0,1]),np.array([1,0]),np.array([1.,2.])
    if problem=="duplicate": q[:]=0; t[:]=1
    elif problem=="negative": q[0]=-1
    elif problem=="outside": t[1]=2
    elif problem=="float_index": q=q.astype(float)
    elif problem=="nan": s[0]=np.nan
    elif problem=="length": s=s[:1]
    elif problem=="shape": q=q.reshape(1,2)
    with pytest.raises(ValueError):
        independently_sorted(2,q,t,s)


def test_empty_hits_and_order_invariance():
    empty=independently_sorted(2,np.array([],dtype=int),np.array([],dtype=int),np.array([]))
    assert all(len(a)==0 for a in empty)
    q,t,s=np.array([0,1]),np.array([1,0]),np.array([.1,.2])
    assert identities(independently_sorted(2,q,t,s))==identities(independently_sorted(2,q[::-1],t[::-1],s[::-1]))


def test_large_unsigned_keys_do_not_promote_to_float():
    q=np.array([2000000000,2000000000],dtype=np.uint64)
    t=np.array([2,1],dtype=np.uint64)
    s=np.array([1.,2.])
    assert independently_sorted(3000000000,q,t,s)[1].tolist()==[1,2]
