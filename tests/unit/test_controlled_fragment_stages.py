"""Check observed stage rules without inferring an unsupported biological cause."""

import json

import numpy as np
import pytest

from benchmark_tools import trace_controlled_fragment_stages as current


def row(**kwargs):
    value=dict(predicted=False, same_seed_group=False, same_candidate=True,
        same_root_hog=True,pair_event="speciation",membership_filter_active=False,detachment_events=[])
    value.update(kwargs)
    return value


@pytest.mark.parametrize("event,active,root,detached,predicted,location", [
    ("speciation",False,False,False,True,"event_rule_retention"),
    ("uncertain",False,True,False,True,"event_rule_retention"),
    ("duplication",False,True,False,False,"observed_duplication_exclusion"),
    ("unambiguous_bypass",False,True,False,True,"unambiguous_bypass_retention"),
    ("speciation",True,True,False,True,"event_rule_retention"),
    ("speciation",True,False,False,False,"root_partition_exclusion"),
    ("speciation",True,False,True,False,"unsupported_satellite_separation"),
])
def test_exact_observed_phylogeny_rules(event,active,root,detached,predicted,location):
    value=row()
    value.update(pair_event=event,membership_filter_active=active,same_root_hog=root,
                 detachment_events=[dict(event_index=0)] if detached else [],predicted=predicted)
    assert current.pair_location(current.METHODS[1],value)==location
    value["predicted"]=not predicted
    with pytest.raises(ValueError): current.pair_location(current.METHODS[1],value)


def test_across_candidate_never_native_pair():
    value=row()
    value.update(same_candidate=False,pair_event=None)
    assert current.pair_location(current.METHODS[1],value)=="candidate_separation"
    value["predicted"]=True
    with pytest.raises(ValueError): current.pair_location(current.METHODS[1],value)


@pytest.mark.parametrize("method",[current.METHODS[0],current.METHODS[3]])
def test_group_derived_identity_and_no_orthology_claim(method):
    value=row()
    assert current.pair_location(method,value) in ("group_separation","mcl_group_separation")
    value["predicted"]=True
    with pytest.raises(ValueError): current.pair_location(method,value)


def test_full_of_native_exclusion_not_a_duplication_call():
    value=row()
    value["same_seed_group"]=True
    assert current.pair_location(current.METHODS[2],value)=="within_mcl_native_pair_exclusion"


def test_group_index_complete_and_disjoint():
    assert current.index(dict(G=["a","b"]),{"a","b"})==dict(a="G",b="G")
    with pytest.raises(ValueError): current.index(dict(G=["a"]),{"a","b"})
    with pytest.raises(ValueError): current.index(dict(G=["a"],H=["a","b"]),{"a","b"})


def test_root_groups_bound_to_source_and_gene_universe(tmp_path):
    path=tmp_path/"roots.tsv"
    path.write_text("root_hog\tsource_family\tgenes\nR0\tF0\ta,b\nR1\tF1\tc\n")
    groups,sources,index=current.root_groups(path,dict(F0={"a","b"},F1={"c"}),{"a","b","c"})
    assert groups==dict(R0=["a","b"],R1=["c"]) and sources==dict(R0="F0",R1="F1")
    assert index==dict(a="R0",b="R0",c="R1")
    path.write_text("root_hog\tsource_family\tgenes\nR0\tF0\ta,c\nR1\tF1\tb\n")
    with pytest.raises(ValueError): current.root_groups(path,dict(F0={"a","b"},F1={"c"}),{"a","b","c"})


def test_original_metrics_have_a_distinct_admission_binding():
    admission=dict(status="admitted",native_validation=dict(metrics=dict(sha256="original")))
    dataset=dict(baseline_admission=[dict(method=current.METHODS[0],admission=admission)])
    assert current.metrics_evidence({},dataset,"baseline",current.METHODS[0])==dict(sha256="original")
    newer=dict(admission=dict(status="admitted",native_validation=dict(metrics=dict(sha256="fragment"))))
    assert current.metrics_evidence(newer,dataset,"fragment",current.METHODS[0])==dict(sha256="fragment")
    with pytest.raises(ValueError): current.metrics_evidence({},dataset,"baseline",current.METHODS[1])


def test_of_unavailable_search_not_false_hit():
    context=dict(seeds=dict(a="G",b="G"))
    value=current.observe(current.METHODS[2],context,("a","b"),False)
    assert value["search_status"]=="unavailable_matched_adapter"
    assert value["hit_forward"] is value["hit_reverse"] is value["graph_direct"] is None
    assert value["pair_event"] is None


def test_source_artifact_binding_and_changed_bytes(tmp_path):
    output=tmp_path/"method"
    output.mkdir()
    artifact=output/"groups.txt"
    artifact.write_text("a b\n")
    ref=current.record(artifact)
    ref=dict(absolute_path=ref["path"],bytes=ref["bytes"],sha256=ref["sha256"])
    status=tmp_path/"status.json"
    status.write_text(json.dumps(dict(methods={current.METHODS[0]:dict(status="process_succeeded",exit_code=0,outputs=[ref])})))
    binding=dict(output=str(output),execution=current.record(status))
    _,reader,_=current.artifact_reader(binding,current.METHODS[0],{})
    assert reader(artifact)==artifact
    artifact.write_text("changed")
    with pytest.raises(ValueError): reader(artifact)


def test_existing_output_guard(tmp_path):
    with pytest.raises(FileExistsError): current.run(tmp_path,tmp_path/"missing.json","0"*64,tmp_path)


def native_fixture(tmp_path, method):
    output=tmp_path/method
    working=output/"orthohmm_working_res"
    checkpoint=working/"high_sensitivity_checkpoint"
    checkpoint.mkdir(parents=True)
    (checkpoint/"gene_names.txt").write_text("a\nb\n")
    for name,values in (("gene_to_species",[0,1]),("hit_queries",[0,1]),("hit_targets",[1,0]),("hit_scores",[15.,16.])):
        np.save(checkpoint/(name+".npy"),np.array(values,dtype="float64" if name=="hit_scores" else "int32"))
    files={p.name:{k:current.record(p)[k] for k in ("bytes","sha256")} for p in checkpoint.iterdir()}
    (checkpoint/"manifest.json").write_text(json.dumps(dict(schema_version=1,complete=True,genes=2,hits=2,files=files)))
    (working/"orthohmm_edges.txt").write_text("a\tb\t1\n")
    (output/"orthohmm_orthogroups.txt").write_text("OG0: a b\n")
    metadata=dict(output_directory=str(output))
    if method==current.METHODS[1]:
        candidate=working/"phylogeny_candidate_superfamilies.txt"
        candidate.write_text("a b\n")
        (working/"phylogeny_candidate_merges.json").write_text("[]")
        (working/"phylogeny_candidate_seeds.tsv").write_text("candidate_family\tseed_families\nFamily0000000\tSeed0000000\n")
        metadata["phylogeny_candidate_profile"]=dict(seed_families=1,candidate_families=1,merges=0)
        phylogeny=output/"orthohmm_phylogeny"
        phylogeny.mkdir()
        (phylogeny/"provenance_manifest.json").write_text(json.dumps(dict(input_cluster_sha256=current.record(candidate)["sha256"],
            pair_orthology_rule="positive_paralogy",root_duplication_rule="species_overlap",membership_reconciliation=None)))
        (phylogeny/"orthohmm_root_hogs.tsv").write_text("root_hog\tsource_family\tgenes\nRootHOG0000000\tFamily0000000\ta,b\n")
        (phylogeny/"orthohmm_reconciliation_nodes.tsv").write_text("\t".join(current.reconciliation.NODE_HEADER)+"\n")
    metric=tmp_path/(method+".json")
    metric.write_text(json.dumps(dict(status="complete",metadata=metadata,counts=dict(significant_hits=2,network_edges=1))))
    refs=[]
    for p in output.rglob("*"):
        if p.is_file():
            r=current.record(p)
            refs.append(dict(absolute_path=r["path"],bytes=r["bytes"],sha256=r["sha256"]))
    status=tmp_path/"status.json"
    status.write_text(json.dumps(dict(methods={method:dict(status="process_succeeded",exit_code=0,outputs=refs)})))
    return dict(output=str(output),execution=current.record(status),configured=dict(metrics=str(metric))),current.record(metric)


@pytest.mark.parametrize("method",current.METHODS[:2])
def test_real_context_checkpoint_graph_and_bypass_loading(tmp_path,method):
    binding,metric=native_fixture(tmp_path,method)
    evidence={}
    context=current.hmm_context(binding,method,{("a","b")},dict(a="A",b="B"),metric,evidence)
    value=current.observe(method,context,("a","b"),True)
    assert value["hit_forward"] is value["hit_reverse"] is value["graph_direct"] is value["graph_connected"] is True
    assert value["directed_hits"]==dict(forward=[dict(row=0,score=15.)],reverse=[dict(row=1,score=16.)])
    if method==current.METHODS[1]:
        assert value["pair_event"]=="unambiguous_bypass" and value["membership_filter_active"] is False
    assert all(current.record(r["path"])==r for r in evidence.values())


def test_changed_graph_bytes_not_silently_reconstructed(tmp_path):
    method=current.METHODS[0]
    binding,metric=native_fixture(tmp_path,method)
    (tmp_path/method/"orthohmm_working_res/orthohmm_edges.txt").write_text("a\tb\t2\n")
    with pytest.raises(ValueError): current.hmm_context(binding,method,{("a","b")},dict(a="A",b="B"),metric,{})
