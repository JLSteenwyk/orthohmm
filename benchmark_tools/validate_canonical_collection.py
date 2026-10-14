"""Verify canonical collection and checkpoints from real sampled native searches."""

import argparse
from collections import Counter
from dataclasses import replace
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
import numpy as np
from benchmark_tools.audit_accuracy_checkpoint import audit
from benchmark_tools.prepare_ob_candidate_neighborhood import check,record
from benchmark_tools.probe_search_stage_cost import selected_names,subset_species,fingerprint_result
from benchmark_tools.run_membership_hillclimb import write_json
from benchmark_tools.run_search_recovery import runtime_records
from orthohmm.accuracy import write_accuracy_checkpoint
from orthohmm.orthohmm import _collect_search_hit_arrays
from orthohmm.search import engine
from orthohmm.search.engine import IndexedSearchResults
from orthohmm.search.sequences import SpeciesSequences


def remap_queries(result,original_indices):
    indices = np.asarray(original_indices,dtype=np.int32)
    if np.any(result.query_indices<0) or np.any(result.query_indices>=len(indices)):
        raise ValueError("Sampled query index outside original mapping")
    return replace(result,query_indices=indices[result.query_indices])


def independent_rows(results,species_ids,gene_to_id,threshold):
    rows = []
    for key in sorted(results):
        value = results[key]
        for q,t,evalue,score in zip(value.query_indices,value.target_indices,value.evalues,value.scores):
            if evalue<threshold:
                rows.append((gene_to_id[species_ids[key[0]][q]],gene_to_id[species_ids[key[1]][t]],float(score)))
    return rows


def run(repo,output,source_revision):
    if output.exists() or os.environ.get("CUDA_VISIBLE_DEVICES")!="" or engine.is_cuda_available():
        raise ValueError("Require new CPU-only canonical collection probe namespace")
    source = Path(__file__).resolve().parent.parent
    revision = subprocess.check_output(["git","rev-parse","--verify",source_revision+"^{commit}"],cwd=repo,text=True).strip()
    paths = [source/"orthohmm/orthohmm.py",source/"orthohmm/accuracy.py",Path(__file__),
        source/"benchmark_tools/hillclimb_canonical_collection_native_probe_20261014.json",
        *sorted((source/"orthohmm/search").rglob("*.py"))]
    sources = [record(path) for path in paths]
    for path,pin in zip(paths,sources):
        blob = subprocess.check_output(["git","show",revision+":"+str(path.relative_to(source))],cwd=repo)
        if hashlib.sha256(blob).hexdigest()!=pin["sha256"]:
            raise ValueError("Wrong frozen collector/search source")
    live_path = repo/"benchmarks/results/full_search_hillclimb_r1_v3_20261012/execution_plan.json"
    live = json.loads(live_path.read_text())
    case = next(case for case in live["cases"] if case["dataset"]=="QfO" and case["label"]=="control100")
    pins = {Path(pin["path"]).name:pin for pin in case["inputs"]}
    names = selected_names(list(pins),3)
    checked = [record(live_path),*[pins[name] for name in names]]
    for name in ("hmm_viterbi.so","kmer_prefilter.so","pair_align.so"):
        relative = Path("orthohmm/search/csrc")/name
        old,new = record(Path(live["source_root"])/relative),record(source/relative)
        if old["sha256"]!=new["sha256"]:
            raise ValueError("Probe does not use the unchanged original CPU libraries")
        checked.extend((old,new))
    for pin in (*checked,*sources): check(pin)
    runtime = runtime_records()
    output.mkdir(parents=True)
    loaded = {name:SpeciesSequences.from_fasta(pins[name]["path"],name) for name in names}
    species_ids = {name:species.ids for name,species in loaded.items()}
    gene_names = sorted(gene for ids in species_ids.values() for gene in ids)
    if len(gene_names)!=len(set(gene_names)):
        raise ValueError("Nonunique selected proteome gene universe")
    gene_to_id = {gene:index for index,gene in enumerate(gene_names)}
    owners = {gene:index for index,name in enumerate(names) for gene in species_ids[name]}
    pairs,searches = {},[]
    for name in names:
        query,indices = subset_species(loaded[name],32)
        for target in names:
            if target==name: continue
            result = engine.search_species_pair_indexed(query,loaded[target],"BLOSUM62",kmer_k=4,
                band_width=64,max_candidates_per_query=100,n_threads=4,prefilter_phase_a_topk=5000)
            pairs[name,target] = remap_queries(result,indices)
            searches.append({"query":name,"target":target,"query_indices":indices,"query_ids":query.ids,
                "result":fingerprint_result(pairs[name,target])})
    reverse = dict(reversed(list(pairs.items())))
    first = _collect_search_hit_arrays(IndexedSearchResults(pairs,species_ids),1e-4,gene_to_id,result_order="canonical")
    second = _collect_search_hit_arrays(IndexedSearchResults(reverse,species_ids),1e-4,gene_to_id,result_order="canonical")
    expected = independent_rows(pairs,species_ids,gene_to_id,1e-4)
    if list(zip(*(values.tolist() for values in first)))!=expected or any(not np.array_equal(a,b) for a,b in zip(first,second)):
        raise ValueError("Canonical collector does not reproduce independent rows/completion-order invariance")
    arrival = _collect_search_hit_arrays(IndexedSearchResults(reverse,species_ids),1e-4,gene_to_id)
    if Counter(zip(*(values.tolist() for values in arrival)))!=Counter(expected):
        raise ValueError("Canonical collection changed the accepted hit multiset")
    checkpoints,audits = [],[]
    for label,arrays in (("forward",first),("reverse",second)):
        directory = output/label
        directory.mkdir()
        checkpoint = write_accuracy_checkpoint(str(directory),gene_names,[owners[gene] for gene in gene_names],*arrays)
        manifest = record(checkpoint/"manifest.json")
        checkpoints.append({"directory":str(checkpoint),"manifest":manifest})
        audits.append(audit(checkpoint,manifest["sha256"]))
    if checkpoints[0]["manifest"]["sha256"]!=checkpoints[1]["manifest"]["sha256"]:
        raise ValueError("Canonical accepted-hit checkpoint hashes differ")
    for pin in (*checked,*sources,*runtime["files"]): check(pin)
    if runtime_records()!=runtime: raise ValueError("Scientific runtime changed")
    report = {"status":"real_native_canonical_collection_and_checkpoints_verified","source_revision":revision,
        "job_id":os.environ.get("SLURM_JOB_ID"),"checked":checked,"sources":sources,"runtime":runtime,
        "searches":searches,"genes":len(gene_names),"accepted_hits":len(expected),"checkpoints":checkpoints,"audits":audits,
        "scope":"Six sampled-query searches against full real targets; real HMM scores/E-values, no fabricated significance values. Query indices are remapped into full species gene tables. Collector/order/checkpoint correctness only, not complete proteome inference or new biological accuracy scores.",
        "defaults_promoted":False,"new_accuracy_scores":False,"goal_completed":False}
    write_json(output/"report.json",report)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo",required=True,type=Path)
    parser.add_argument("--output",required=True,type=Path)
    parser.add_argument("--source-revision",required=True)
    args = parser.parse_args()
    run(args.repo.resolve(),args.output.resolve(),args.source_revision)
