"""Reconstruct canonical inputs and independently check candidate outputs."""

import argparse
import hashlib
import json
from pathlib import Path
import pickle
import shlex
import subprocess

import numpy as np

from benchmark_tools.prepare_ob_candidate_neighborhood import record,check
from benchmark_tools.probe_installed_ob_clustering import write_json
from benchmark_tools.probe_ob_candidate_order_scores import LABELS
from benchmark_tools.audit_ob_candidate_order_scores import compare_all,validate_summary
from benchmark_tools.audit_ob_dependency_replay import candidate_readback
from benchmark_tools.compare_installed_ob_search import partition
from benchmark_tools.audit_installed_orthobench import compare_partitions

PLAN_SHA = "6925ccc4684e1687a1f03ee8589d45d0f72b53605be5406f2a975c49acc8ecb7"


def independently_sorted(n, q, t, s):
    if n <= 0 or n > 3037000499:
        raise ValueError("Gene universe cannot be encoded safely")
    if (any(a.ndim != 1 for a in (q,t,s)) or not len(q)==len(t)==len(s)
            or q.dtype.kind not in "iu" or t.dtype.kind not in "iu" or s.dtype.kind not in "fiu"
            or not np.isfinite(s).all() or np.any(s<=0)):
        raise ValueError("Invalid numeric hits")
    if len(q) and (min(q.min(),t.min())<0 or max(q.max(),t.max())>=n):
        raise ValueError("Out-of-universe hit")
    # Independent integer-key sort, not the adapter's lexsort implementation.
    keys = q.astype(np.int64)*n+t.astype(np.int64)
    order = np.argsort(keys)
    if len(order)>1 and np.any(np.diff(keys[order])==0):
        raise ValueError("Duplicate directed key")
    return q[order],t[order],s[order]


def identities(arrays):
    return [dict(shape=list(a.shape),dtype=a.dtype.str,
                 sha256=hashlib.sha256(a.tobytes(order="C")).hexdigest()) for a in arrays]


def input_readback(plan):
    check(plan["cache"])
    # Only the previously admitted, hash-verified local cache is deserialized.
    with Path(plan["cache"]["path"]).open("rb") as stream:
        payload = pickle.load(stream)
    names = sorted(payload["all_gene_ids"])
    if len(names)!=len(set(names)):
        raise ValueError("Duplicate historical genes")
    ids = {g:i for i,g in enumerate(names)}
    hits = payload["all_hits"]
    old = (np.fromiter((ids[a] for a,b in hits),dtype=np.int32,count=len(hits)),
           np.fromiter((ids[b] for a,b in hits),dtype=np.int32,count=len(hits)),
           np.fromiter(hits.values(),dtype=np.float64,count=len(hits)))
    del hits,payload
    if np.any(old[0]==old[1]):
        raise ValueError("Unexpected historical self hits")
    old_records = identities(independently_sorted(len(names),*old))
    del old
    checkpoint = Path(plan["checkpoint"])
    if names != (checkpoint / "gene_names.txt").read_text().splitlines():
        raise ValueError("Gene indexing differs")
    arrays = tuple(np.load(checkpoint / f"hit_{name}.npy",allow_pickle=False,mmap_mode="r")
                   for name in ("queries","targets","scores"))
    full = independently_sorted(len(names),*arrays)
    full_records = identities(full)
    nonself = full[0]!=full[1]
    new_records = identities(tuple(a[nonself] for a in full))
    if old_records[:2] != new_records[:2]:
        raise ValueError("Historical/fresh nonself keys differ")
    return dict(policy="canonical_directed_pair_v1", checks=dict(historical_scores_order_invariant=True,
                fresh_scores_order_invariant=True,self_control_nonself_equal=True),
        canonical_arrays=dict(zip(LABELS,(old_records,new_records,old_records,new_records,full_records))),
        self_hits=int(np.count_nonzero(~nonself)),scores_rounded=False,scoring=False)


def audit(repo,directory):
    records=[]
    def read(path,sha=None):
        item=record(path)
        if sha is not None and item["sha256"]!=sha:
            raise ValueError("Changed pinned report")
        records.append(item)
        return json.loads(Path(path).read_text())
    path=directory / "plan.json"
    plan=read(path,PLAN_SHA)
    records.extend(plan["checked_records"])
    for item in records: check(item)
    if plan["fixture"] is not False or plan["output"]!=str(directory) or plan["repo"]!=str(repo):
        raise ValueError("Wrong full-data plan scope")
    scheduler_command=["sacct","-j","22323","--format=JobID,State,ExitCode,Elapsed","-n","-P"]
    scheduler=subprocess.check_output(scheduler_command,text=True,timeout=60)
    status=[line.split("|") for line in scheduler.splitlines() if line.startswith("22323|")]
    if len(status)!=1 or status[0][1:3]!=["COMPLETED","0:0"] or (directory/"failure.json").exists():
        raise ValueError("Native execution not successfully complete")
    submission,started=read(directory/"submission.json"),read(directory/"started.json")
    report,wrapper=read(directory/"report.json"),read(directory/"wrapper_complete.json")
    expected=["env","-i","PATH=/usr/bin:/bin","HOME=/home/bizon","LANG=C.UTF-8",
              *[k+"="+v for k,v in plan["env"].items()],plan["runtime"]["executable"],"-I",
              plan["wrapper"]["path"],"--run",str(path),"--plan-sha256",PLAN_SHA]
    if (submission["job_id"]!="22323" or submission["plan"]!=record(path)
            or submission["command"][-2]!="--wrap" or shlex.split(submission["command"][-1])!=expected
            or started["command"]!=expected[-5:] or started["plan"]!=record(path)
            or started["runtime"]!={k:plan["runtime"][k] for k in ("versions","sources","native")}
            or started["job_id"] is not None):
        raise ValueError("Submitted execution or runtime differs")
    validate_summary(report,record(path),plan["source"])
    if (wrapper["status"]!="canonical_candidate_native_complete_pending_readback"
            or wrapper["plan"]!=record(path) or wrapper["wrapper"]!=plan["wrapper"]
            or wrapper["native_report"]!=record(directory/"report.json")
            or wrapper["canonical_inputs"]!=record(directory/"canonical_inputs.json")
            or wrapper["production_changed"] is not False):
        raise ValueError("Adapter completion identity differs")
    recorded_inputs=read(directory/"canonical_inputs.json")
    reconstructed=input_readback(plan)
    if reconstructed!=recorded_inputs:
        raise ValueError("Independent canonical inputs differ")
    parent=read(repo/"benchmark_tools/results/ob_candidate_order_score_readback_22322.json",
                "1fcfa5ebd0d99c8a859345c6f269392a3b4738396f81920ed838825f878304dd")
    parent_records={r["path"]:r for r in parent["checked_records"]}
    stage_path=repo/"benchmarks/work/ob_dependency_replay_v2_20260926/leiden011/stage_report.json"
    prior=read(stage_path,parent_records[str(stage_path)]["sha256"])
    universe=set((Path(plan["checkpoint"])/"gene_names.txt").read_text().splitlines())
    groups,rows={},[]
    for row in report["rows"]:
        label=row["label"]
        root=directory/label
        work=root/"replay/orthohmm_working_res"
        if (read(root/"result.json")!=row or row["output_records"]!=[record(p) for p in sorted(work.iterdir())]
                or row["prediction"]!=record(work/"phylogeny_candidate_superfamilies.txt")
                or row["candidate_summary"]["parameters"]!=prior["candidates"]["parameters"]):
            raise ValueError("Output or parameter identity differs")
        evidence=candidate_readback(root,dict(candidates=row["candidate_summary"]),plan["seed"],universe)
        groups[label]=partition(Path(row["prediction"]["path"]),universe)
        records.extend(row["output_records"])
        rows.append(dict(label=label,prediction=row["prediction"],readback=evidence))
    comparisons=compare_all(groups)
    if comparisons!=report["comparisons"]:
        raise ValueError("Native contrasts do not reproduce")
    baselines={}
    for old in parent["rows"]:
        item=old["prediction"]
        check(item)
        records.append(item)
        baseline=partition(Path(item["path"]),universe)
        baselines[old["label"]]={k:compare_partitions(baseline,v) for k,v in groups.items()}
    records.append(record(directory/"slurm-22323.log"))
    for item in records: check(item)
    return dict(status="canonical_candidate_inputs_and_outputs_read_back",source=record(__file__),
        checked_records=records,canonical_input_reconstruction=reconstructed,rows=rows,comparisons=comparisons,
        original_factorial_comparisons=baselines,scheduler=dict(command=scheduler_command,output=scheduler),
        accuracy_evaluated=False,production_changed=False,
        limitations=["Fixed runtime, score vectors, gene indexing and seeds; not universal reproducibility.",
            "Independent integer-key input sort and candidate seed/merge reconstruction, not recalculated search support.",
            "No phylogenetic inference, biological accuracy or controlled timing comparison."])


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo",type=Path,required=True)
    p.add_argument("--directory",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    if a.output.exists(): raise FileExistsError(a.output)
    write_json(a.output,audit(a.repo.resolve(),a.directory.resolve()))
