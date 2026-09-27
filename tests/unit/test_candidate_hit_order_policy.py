import numpy as np
import pytest

from benchmark_tools.candidate_hit_order_policy import (
    canonical_hit_order_v1, fixture_inputs, fixture_orders,
)
from orthohmm.refinement import merge_supported_satellite_candidate_clusters


def test_canonicalization_preserves_values_self_hits_dtypes_and_inputs():
    q, t, s = np.array([1,0,0],dtype=np.int32),np.array([0,1,0],dtype=np.int32),np.array([.2,.1,3.])
    before = [a.copy() for a in (q,t,s)]
    actual = canonical_hit_order_v1(["a","b"],q,t,s)
    assert [a.tolist() for a in actual] == [[0,0,1],[0,1,0],[3.,.1,.2]]
    for a,b,c in zip((q,t,s),before,actual):
        assert np.array_equal(a,b) and a.dtype == c.dtype
    assert all(np.array_equal(a,b) for a,b in zip(actual,canonical_hit_order_v1(["a","b"],*actual)))


def test_empty_hits_allowed_with_defined_universe():
    result = canonical_hit_order_v1(["a"],np.array([],dtype=int),np.array([],dtype=int),np.array([]))
    assert all(len(a)==0 for a in result)


@pytest.mark.parametrize("names", [[],["b","a"],["a","a"],[""],[1]])
def test_invalid_indexing_rejected(names):
    with pytest.raises(ValueError):
        canonical_hit_order_v1(names,np.array([0]),np.array([0]),np.array([1.]))


@pytest.mark.parametrize("problem", ["shape","length","float_index","negative","outside","nan","inf","zero","duplicate","complex"])
def test_invalid_hits_rejected(problem):
    q,t,s = np.array([0,1]),np.array([1,0]),np.array([1.,2.])
    if problem == "shape": q=q.reshape(1,2)
    elif problem == "length": t=t[:1]
    elif problem == "float_index": q=q.astype(float)
    elif problem == "negative": q[0]=-1
    elif problem == "outside": t[0]=2
    elif problem in ("nan","inf","zero"): s[0]={"nan":np.nan,"inf":np.inf,"zero":0.}[problem]
    elif problem == "duplicate": q[:]=0; t[:]=1
    elif problem == "complex": s=s.astype(complex)
    with pytest.raises(ValueError):
        canonical_hit_order_v1(["a","b"],q,t,s)


def test_36_orders_reproduce_sensitivity_and_experimental_invariance():
    inputs = fixture_inputs()
    baseline, canonical = set(),set()
    for order in fixture_orders():
        arrays = tuple(a[order] for a in inputs)
        for seen, values in ((baseline,arrays),(canonical,canonical_hit_order_v1(list("abcde"),*arrays))):
            groups,merges,_,_ = merge_supported_satellite_candidate_clusters(
                [[0,1,2],[3],[4]],*values,np.arange(5),max_satellites_per_anchor=1,max_iterations=1)
            assert merges == 1
            seen.add(tuple(sorted(tuple(sorted(g)) for g in groups)))
    assert len(baseline)==2
    assert len(canonical)==1
    assert canonical <= baseline


def test_order_can_change_float64_cluster_total():
    a = np.add.reduceat(np.array([.1,.2,.3]),[0])[0]
    b = np.add.reduceat(np.array([.3,.2,.1]),[0])[0]
    assert a != b and abs(a-b) < 1e-15
