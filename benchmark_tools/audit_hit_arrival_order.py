"""Audit order sensitivity using unchanged helpers and complete accepted-hit data."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
import numpy as np
from benchmark_tools.audit_accuracy_checkpoint import audit
from benchmark_tools.prepare_ob_candidate_neighborhood import check,record
from benchmark_tools.run_membership_hillclimb import write_json
from benchmark_tools.run_search_recovery import runtime_records
from orthohmm.accuracy import build_rbnh_edges,build_singleton_assignment_edges,load_accuracy_checkpoint,read_index_clusters
from orthohmm.externals import execute_leiden


def block_layout(species,queries,targets):
    if not len(queries) or queries.shape != targets.shape:
        raise ValueError("Require nonempty aligned hit endpoints")
    base = int(np.max(species))+1
    keys = species[queries].astype(np.int64)*base+species[targets]
    starts = np.r_[0,np.flatnonzero(keys[1:]!=keys[:-1])+1]
    ends = np.r_[starts[1:],len(keys)]
    block_keys = keys[starts]
    if len(np.unique(block_keys)) != len(block_keys):
        raise ValueError("Species-pair hit blocks must be contiguous and unique")
    return starts,ends,block_keys


def block_permutation(starts,ends,keys,reverse=False):
    order = np.argsort(keys,kind="stable")
    if reverse:
        order = order[::-1]
    return np.concatenate([np.arange(starts[i],ends[i],dtype=np.int64) for i in order])


def assignments(edges,clusters):
    n = len(edges.gene_names)
    membership = np.full(n,-1,dtype=np.int32)
    singleton = np.zeros(n,dtype=bool)
    for index,genes in enumerate(clusters):
        membership[genes] = index
        singleton[genes] = len(genes)==1
    if np.any(membership<0):
        raise ValueError("Initial partition lost genes")
    output = {}
    for a,b in zip(edges.sources,edges.targets):
        if singleton[a] == singleton[b]:
            raise ValueError("Assignment edge does not connect singleton to nonsingleton")
        query,target = (a,b) if singleton[a] else (b,a)
        query,cluster = int(query),int(membership[target])
        if query in output and output[query]!=cluster:
            raise ValueError("Singleton assigned to multiple clusters")
        output[query] = cluster
    return output


def edge_identity(edges):
    return {"edges":len(edges),"arrays":[hashlib.sha256(a.tobytes()).hexdigest()
        for a in (edges.sources,edges.targets,edges.weights)]}


def toy_reproducer():
    names = ["s","a1","a2","b1","b2"]
    clusters = [[0],[1,2],[3,4]]
    q,t,s = np.array([0,0],dtype=np.int32),np.array([1,3],dtype=np.int32),np.array([5.,5.])
    first = assignments(build_singleton_assignment_edges(names,clusters,q,t,s),clusters)
    reverse = assignments(build_singleton_assignment_edges(names,clusters,q[::-1],t[::-1],s[::-1]),clusters)
    if first==reverse:
        raise ValueError("Unchanged helper no longer demonstrates the prescribed exact tie")
    return {"original_assignment":first[0],"reversed_assignment":reverse[0],"same_input_hits":True}


def initial_partition(output,edges):
    working = output/"orthohmm_working_res"
    working.mkdir(exist_ok=True)
    execute_leiden(.1,str(output),edges=edges,include_isolates=True,seed=4)
    return read_index_clusters(str(working/"orthohmm_edges_clustered.txt"),edges.gene_names)


def run(repo,output):
    if output.exists() or os.environ.get("CUDA_VISIBLE_DEVICES")!="":
        raise ValueError("Require new CPU-only order audit namespace")
    started = time.monotonic()
    baseline_path = repo/"benchmarks/results/full_search_hillclimb_r1_v3_20261012/OrthoBench/results.json"
    baseline = json.loads(baseline_path.read_text())
    control = next(row for row in baseline["arms"] if row["label"]=="control100")
    if baseline["status"]!="finished_all_arms_scored" or control["status"]!="scored_development_result":
        raise ValueError("Completed native baseline required")
    check(baseline["execution_plan"])
    plan = json.loads(Path(baseline["execution_plan"]["path"]).read_text())
    case = next(row for row in plan["cases"] if row["dataset"]=="OrthoBench" and row["label"]=="control100")
    directory = Path(case["output"])/"output/orthohmm_working_res/high_sensitivity_checkpoint"
    checkpoint_pin = control["native"]["numeric_checkpoint"]["manifest"]
    check(checkpoint_pin)
    checkpoint_audit = audit(directory,checkpoint_pin["sha256"])
    names,species,q,t,scores = load_accuracy_checkpoint(directory,verify=True)
    if (len(names),len(np.unique(species)),len(scores))!=(251378,12,18486508):
        raise ValueError("Full native checkpoint universe changed")
    manifest = json.loads((directory/"manifest.json").read_text())
    checked = [record(baseline_path),baseline["execution_plan"],checkpoint_pin]
    for name,identity in manifest["files"].items():
        pin = record(directory/name)
        if pin["sha256"]!=identity["sha256"] or pin["bytes"]!=identity["bytes"]:
            raise ValueError("Checkpoint source identity changed")
        checked.append(pin)
    source = Path(__file__).resolve().parent.parent
    sources = [record(source/path) for path in ("orthohmm/accuracy.py","orthohmm/externals.py","orthohmm/leiden_worker.py",
        "benchmark_tools/audit_hit_arrival_order.py","benchmark_tools/hillclimb_hit_arrival_order_audit_20261014.json")]
    runtime = runtime_records()
    output.mkdir(parents=True)
    result = {"status":"hit_arrival_order_audit_running","job_id":os.environ.get("SLURM_JOB_ID"),
        "checked":checked,"sources":sources,"runtime":runtime,"checkpoint_audit":checkpoint_audit,"toy":toy_reproducer(),
        "rows":[],"defaults_promoted":False,"goal_completed":False}
    write_json(output/"report.json",result)
    try:
        starts,ends,keys = block_layout(species,q,t)
        reference = build_rbnh_edges(names,species,q,t,scores)
        clusters = initial_partition(output,reference)
        original_assignments = None
        reference_identity = edge_identity(reference)
        for mode in ("arrival","canonical_blocks","reversed_blocks"):
            if mode=="arrival":
                arrays = (q,t,scores)
            else:
                permutation = block_permutation(starts,ends,keys,mode=="reversed_blocks")
                if len(permutation)!=len(q) or len(np.unique(permutation))!=len(q):
                    raise ValueError("Block permutation lost/duplicated hits")
                arrays = tuple(values[permutation] for values in (q,t,scores))
            graph = reference if mode=="arrival" else build_rbnh_edges(names,species,*arrays)
            if edge_identity(graph)!=reference_identity:
                raise ValueError("Block-only permutation changed RBNH graph")
            edges = build_singleton_assignment_edges(names,clusters,*arrays)
            selected = assignments(edges,clusters)
            if original_assignments is None:
                original_assignments = selected
            changed = [{"gene":names[gene],"original_cluster":original_assignments.get(gene),
                        "permuted_cluster":selected.get(gene)} for gene in sorted(original_assignments.keys()|selected.keys())
                       if original_assignments.get(gene)!=selected.get(gene)]
            result["rows"].append({"mode":mode,"rbnh_graph_unchanged":True,"singleton_edges":edge_identity(edges),
                "assigned_singletons":len(selected),"changed_assignments":changed})
            write_json(output/"report.json",result)
        for pin in (*checked,*sources,*runtime["files"]):
            check(pin)
        if runtime_records()!=runtime:
            raise ValueError("Scientific runtime changed during order audit")
        result.update(status="complete_native_hit_block_order_audit",species_pair_blocks=len(keys),
            initial_rbnh_graph=reference_identity,initial_clusters=len(clusters),wall_s=time.monotonic()-started,
            scope="Unchanged initial graph/singleton helpers on a fixed accepted-hit checkpoint. No full inference, new scores, historical-difference causal attribution or default promotion.")
    except Exception as error:
        result.update(status="failed_native_hit_block_order_audit",error_type=type(error).__name__,error=str(error))
    write_json(output/"report.json",result)
    return result


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args = parser.parse_args()
    result = run(args.repo.resolve(),args.output.resolve())
    print(result["status"])
    raise SystemExit(0 if result["status"]=="complete_native_hit_block_order_audit" else 1)
