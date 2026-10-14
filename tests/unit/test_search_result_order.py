from collections import Counter
import json
from types import SimpleNamespace

import numpy as np
import pytest

from orthohmm.accuracy import build_singleton_assignment_edges,validate_search_result_order
from orthohmm.args_processing import process_args
from orthohmm.helpers import StartStep,SubstitutionMatrix
from orthohmm.orthohmm import _collect_search_hit_arrays,execute
from orthohmm.parser import create_parser
from orthohmm.search.engine import IndexedSearchResults,IndexedSpeciesPairResults


def results(reverse=False):
    rows = [(('Q','Beta'),IndexedSpeciesPairResults('Q','Beta',np.array([0,1],dtype=np.int32),
             np.array([0,0],dtype=np.int32),np.array([1e-6,1e-6]),np.array([5.,4.]))),
            (('Q','Alpha'),IndexedSpeciesPairResults('Q','Alpha',np.array([1,0],dtype=np.int32),
             np.array([0,0],dtype=np.int32),np.array([1e-6,1e-6]),np.array([3.,5.])))]
    return IndexedSearchResults(dict(rows[::-1] if reverse else rows),{'Q':['s'],'Alpha':['a1','a2'],'Beta':['b1','b2']})


def collect(search,order='arrival'):
    return _collect_search_hit_arrays(search,1e-4,{gene:i for i,gene in enumerate(['s','a1','a2','b1','b2'])},result_order=order)


def records(arrays):
    return Counter(zip(*(a.tolist() for a in arrays)))


def test_canonical_collection_is_insertion_order_independent_and_preserves_internal_rows_and_hits():
    first,second = collect(results(),'canonical'),collect(results(True),'canonical')
    assert first[1].tolist()==[2,1,3,4]
    for a,b in zip(first,second):
        np.testing.assert_array_equal(a,b)
    assert records(first)==records(collect(results()))==records(collect(results(True)))
    assert first[0].dtype==first[1].dtype==np.int32 and first[2].dtype==np.float64


def test_default_arrival_is_preserved_and_canonical_exact_tie_assignment_is_stable():
    names = ['s','a1','a2','b1','b2']
    clusters = [[0],[1,2],[3,4]]
    def edges(data):
        result = build_singleton_assignment_edges(names,clusters,*data)
        return list(zip(result.sources.tolist(),result.targets.tolist(),result.weights.tolist()))
    assert edges(collect(results()))!=edges(collect(results(True)))
    assert edges(collect(results(),'canonical'))==edges(collect(results(True),'canonical'))


def test_canonical_collection_preserves_mask_boundary_and_hit_multiplicity():
    search = results()
    row = search.pair_results[('Q','Alpha')]
    row.query_indices = np.zeros(5,dtype=np.int32)
    row.target_indices = np.zeros(5,dtype=np.int32)
    row.evalues = np.array([0.,0.,1e-4,float('nan'),float('inf')])
    row.scores = np.array([5.,5.,9.,10.,11.])
    collected = collect(search,'canonical')
    assert collected[1].tolist()==[1,1,3,4]
    assert records(collected)[(0,1,5.)]==2


def test_empty_indexed_results_remain_typed_empty_arrays():
    arrays = collect(IndexedSearchResults({},{}),'canonical')
    assert [len(a) for a in arrays]==[0,0,0]
    assert [a.dtype for a in arrays]==[np.dtype('int32'),np.dtype('int32'),np.dtype('float64')]


def test_legacy_named_results_keep_default_behavior_and_reject_canonical():
    search = {'pair':SimpleNamespace(query_names=['s'],target_names=['a1'],evalues=[0.],scores=[5.])}
    assert _collect_search_hit_arrays(search,1e-4,{'s':0,'a1':1})==([0],[1],[5.])
    with pytest.raises(ValueError,match='indexed'):
        _collect_search_hit_arrays(search,1e-4,{'s':0,'a1':1},result_order='canonical')


@pytest.mark.parametrize('order',[None,True,[],{},'unknown'])
def test_invalid_collection_modes_are_explicit_errors(order):
    with pytest.raises(ValueError,match='collection order'):
        validate_search_result_order(order)
    with pytest.raises(ValueError):
        collect(results(),order)


def test_cli_is_opt_in_and_handoff_keeps_default_call_shape(tmp_path):
    parser = create_parser()
    assert parser.parse_args([str(tmp_path)]).search_result_order=='arrival'
    assert 'search_result_order' not in process_args(parser.parse_args([str(tmp_path)]))
    configured = process_args(parser.parse_args([str(tmp_path),'--accuracy_profile','high_sensitivity','--search_result_order','canonical']))
    assert configured['search_result_order']=='canonical'


@pytest.mark.parametrize('extra',[
    [],['--accuracy_profile','high_sensitivity','--start','search_res'],
    ['--accuracy_profile','high_sensitivity','--search_mode','phmmer'],
])
def test_cli_rejects_inactive_canonical_configurations(tmp_path,monkeypatch,extra):
    import orthohmm.args_processing as processing
    monkeypatch.setattr(processing,'find_executable',lambda name:name)
    with pytest.raises(SystemExit):
        process_args(create_parser().parse_args([str(tmp_path),'--search_result_order','canonical',*extra]))


@pytest.mark.parametrize('profile,search,start',[
    ('standard','builtin',None),('high_sensitivity','phmmer',None),('high_sensitivity','builtin',StartStep.search_res),
])
def test_api_rejects_inactive_canonical_configurations_before_search(tmp_path,profile,search,start):
    with pytest.raises(ValueError,match='fresh built-in high_sensitivity'):
        execute(str(tmp_path),str(tmp_path/'out'),None,2,.5,'mcl',1.5,start,None,
            SubstitutionMatrix.blosum62,1e-4,accuracy_profile=profile,search_mode=search,search_result_order='canonical')


def test_pipeline_records_actual_canonical_policy_without_changing_search_settings(tmp_path,monkeypatch):
    import orthohmm.search.engine as engine
    class SearchReached(Exception): pass
    captured = {}
    def search(*args,**kwargs):
        captured.update(kwargs)
        raise SearchReached()
    monkeypatch.setattr(engine,'execute_builtin_search',search)
    (tmp_path/'species.faa').write_text('>gene\nACDEFGHIK\n')
    metrics = tmp_path/'metrics.json'
    with pytest.raises(SearchReached):
        execute(str(tmp_path),str(tmp_path/'out'),None,2,.5,'mcl',1.5,None,None,
            SubstitutionMatrix.blosum62,1e-4,accuracy_profile='high_sensitivity',search_result_order='canonical',metrics_json=str(metrics))
    assert 'search_result_order' not in captured
    assert captured['kmer_k']==4 and captured['max_candidates_per_query']==100
    assert json.loads(metrics.read_text())['metadata']['search_result_order']=='canonical'


def test_pipeline_uses_canonical_collector_before_writing_checkpoint(tmp_path,monkeypatch):
    import orthohmm.orthohmm as pipeline
    import orthohmm.search.engine as engine
    (tmp_path/'Q.fa').write_text('>s\nACDEFGHIK\n')
    (tmp_path/'Alpha.fa').write_text('>a1\nACDEFGHIK\n>a2\nACDEFGHIK\n')
    (tmp_path/'Beta.fa').write_text('>b1\nACDEFGHIK\n>b2\nACDEFGHIK\n')
    monkeypatch.setattr(engine,'execute_builtin_search',lambda *args,**kwargs:results())
    captured = {}
    class CheckpointReached(Exception): pass
    def checkpoint(output,names,species,q,t,scores):
        captured['targets'] = [names[index] for index in t]
        captured['queries'] = [names[index] for index in q]
        raise CheckpointReached()
    monkeypatch.setattr(pipeline,'write_accuracy_checkpoint',checkpoint)
    with pytest.raises(CheckpointReached):
        execute(str(tmp_path),str(tmp_path/'out'),None,2,.5,'mcl',1.5,None,None,
            SubstitutionMatrix.blosum62,1e-4,accuracy_profile='high_sensitivity',search_result_order='canonical')
    assert captured['targets']==['a2','a1','b1','b2']
    assert captured['queries']==['s']*4


def test_benchmark_wrapper_rejects_inactive_canonical_mode():
    from benchmark_tools.benchmark_production import parse_args
    with pytest.raises(SystemExit):
        parse_args(['inputs','output','metrics.json','--cpu','1','--search-result-order','canonical'])


def test_refinement_recollection_uses_the_same_policy_without_forwarding_it_to_graph_rules(tmp_path,monkeypatch):
    import orthohmm.orthohmm as pipeline
    names = ['s','a1','a2','b1','b2']
    working = tmp_path/'orthohmm_working_res'
    working.mkdir()
    (working/'orthohmm_edges_clustered.txt').write_text('s\na1 a2\nb1 b2\n')
    lengths = np.array([('Q' if name=='s' else 'Alpha' if name.startswith('a') else 'Beta',name,10)
                       for name in names],dtype=[('spp',object),('name',object),('length',int)])
    captured = {}
    def refine(clusters,q,t,scores,*args,**kwargs):
        captured['targets'] = [names[index] for index in t]
        assert 'search_result_order' not in kwargs
        return clusters
    monkeypatch.setattr(pipeline,'refine_cluster_indices',refine)
    pipeline._refine_cluster_file(str(tmp_path),lengths,results(),{},1e-4,search_result_order='canonical')
    assert captured['targets']==['a2','a1','b1','b2']
