import numpy as np
import pytest

from orthohmm.refinement import merge_supported_satellite_candidate_clusters












def test_order_can_change_float64_cluster_total():
    a = np.add.reduceat(np.array([.1,.2,.3]),[0])[0]
    b = np.add.reduceat(np.array([.3,.2,.1]),[0])[0]
    assert a != b and abs(a-b) < 1e-15
