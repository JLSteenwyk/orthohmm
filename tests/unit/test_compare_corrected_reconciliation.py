import hashlib
import json

import pytest

from benchmark_tools.compare_corrected_reconciliation import compare_trees, load, compare


@pytest.mark.parametrize('left,right,distance', [
    ('((A,B),(C,D));', '((D,C),(B,A));', 0),
    ('((A:1,B:2):3,(C:4,D:5):6);', '((A:9,B:8),(C:7,D:6));', 0),
    ('((A,B),(C,D));', '((A,C),(B,D));', 4),
    ('((A,B),C,D);', '((A,B),(C,D));', 1),
])
def test_rooted_topology_comparison(tmp_path, left, right, distance):
    a,b=tmp_path/'a.nwk',tmp_path/'b.nwk'
    a.write_text(left); b.write_text(right)
    result=compare_trees(a,b)
    assert result['species']==4
    assert result['rooted_symmetric_difference']==distance
    assert result['left_only_clades']+result['right_only_clades']==distance
    assert not result['either_tree_is_ground_truth']


def test_mismatched_taxa_rejected(tmp_path):
    a,b=tmp_path/'a.nwk',tmp_path/'b.nwk'
    a.write_text('((A,B),(C,D));'); b.write_text('((A,B),(C,E));')
    with pytest.raises(ValueError,match='universes'):
        compare_trees(a,b)


def test_checksum_rejection(tmp_path):
    p=tmp_path/'data.json'; p.write_text('{}')
    with pytest.raises(ValueError,match='checksum'):
        load(p,'0'*64)


def test_nonadmitted_input_rejected(tmp_path):
    p=tmp_path/'data.json'; p.write_text(json.dumps(dict(status='not_admitted')))
    pin=(p,hashlib.sha256(p.read_bytes()).hexdigest())
    with pytest.raises(ValueError,match='admitted'):
        compare(pin,pin)


def admission(tmp_path, name, provenance='shared', pair_count=2):
    tree=tmp_path/(name+'.nwk'); tree.write_text('((A,B),(C,D));')
    metrics=tmp_path/(name+'_metrics.json')
    metrics.write_text(json.dumps(dict(counts=dict(ortholog_pairs=2,candidate_families=3,root_hogs=4))))
    def record(p):
        return dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
    report=dict(status='corrected_qfo_native_pair_output_verified',cell=name,native_pair_count=pair_count,
                native_group_integrity=dict(native_metrics=record(metrics),species_tree=record(tree),
                    partition=dict(candidate_families=3,root_hogs=4),membership=None))
    for key in ('candidate_admission','environment','prepared'):
        report[key]=dict(sha256=provenance)
    path=tmp_path/(name+'.json'); path.write_text(json.dumps(report))
    return path,hashlib.sha256(path.read_bytes()).hexdigest()


def test_comparison_is_diagnostic_not_accuracy(tmp_path):
    result=compare(admission(tmp_path,'left'),admission(tmp_path,'right'))
    assert result['tree_comparison']['rooted_symmetric_difference']==0
    assert [r['native_pairs'] for r in result['rows']]==[2,2]
    assert not result['accuracy_evaluated']
    assert not result['publication_ready']


@pytest.mark.parametrize('fault',['provenance','counts','tree'])
def test_inconsistent_comparison_rejected(tmp_path,fault):
    left=admission(tmp_path,'left')
    right=admission(tmp_path,'right',provenance='other' if fault=='provenance' else 'shared',
                    pair_count=3 if fault=='counts' else 2)
    if fault=='tree':
        (tmp_path/'right.nwk').write_text('(A,B,C,D);')
    with pytest.raises(ValueError):
        compare(left,right)
