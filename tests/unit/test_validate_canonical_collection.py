import numpy as np
import pytest

from benchmark_tools.validate_canonical_collection import remap_queries,independent_rows
from orthohmm.search.engine import IndexedSpeciesPairResults


def result():
    return IndexedSpeciesPairResults('q','t',np.array([1,0],dtype=np.int32),np.array([0,1],dtype=np.int32),
        np.array([0.,1e-4]),np.array([5.,4.]))


def test_query_remapping_preserves_raw_scores_evalues_targets_and_original_result():
    original = result()
    remapped = remap_queries(original,[4,2])
    assert remapped.query_indices.tolist()==[4,2] and original.query_indices.tolist()==[0,1]
    assert remapped.scores is original.scores and remapped.evalues is original.evalues
    assert remapped.target_indices is original.target_indices


def test_invalid_sample_index_is_rejected():
    with pytest.raises(ValueError,match='outside'):
        remap_queries(result(),[2])


def test_independent_row_reader_preserves_boundary_filtering_and_index_tables():
    pair = result()
    rows = independent_rows({('q','t'):pair},{'q':['qa','qb'],'t':['ta','tb']},
        {'qa':0,'qb':1,'ta':2,'tb':3},1e-4)
    assert rows==[(0,3,5.)]
