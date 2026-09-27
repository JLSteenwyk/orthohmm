"""Compare the directed gene-set constraints consumed by membership filtering."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from benchmark_tools.prepare_ob_candidate_neighborhood import record,check
from benchmark_tools.probe_installed_ob_clustering import write_json


def semantic_constraints(events):
    if not isinstance(events,list):
        raise ValueError("Constraints must be a list")
    result=[]
    for event in events:
        sides=[]
        for name in ("source_genes","target_genes"):
            genes=event[name]
            if (not isinstance(genes,list) or not genes
                    or any(not isinstance(g,str) or not g for g in genes)
                    or len(set(genes))!=len(genes)):
                raise ValueError("Require nonempty unique named genes")
            sides.append(tuple(sorted(genes)))
        if set(sides[0]) & set(sides[1]):
            raise ValueError("Source and target overlap")
        result.append(tuple(sides))
    return result


def compare_constraints(old,new):
    a,b=semantic_constraints(old),semantic_constraints(new)
    ca,cb=Counter(a),Counter(b)
    digest=lambda values: hashlib.sha256(json.dumps(sorted(values),separators=(",",":"),ensure_ascii=True).encode()).hexdigest()
    return dict(historical_constraints=len(a),current_constraints=len(b),semantic_multiset_equal=ca==cb,
        historical_only=sum((ca-cb).values()),current_only=sum((cb-ca).values()),
        semantic_sequence_equal=a==b,positions_changed=sum(x!=y for x,y in zip(a,b))+abs(len(a)-len(b)),
        historical_semantic_sha256=digest(a),current_semantic_sha256=digest(b))


def audit(repo):
    records=[]
    def read(path,sha):
        item=record(path)
        if item["sha256"]!=sha:
            raise ValueError("Changed retained evidence")
        records.append(item)
        return json.loads(Path(path).read_text())
    parent=read(repo/"benchmark_tools/results/ob_canonical_candidate_readback_22323.json",
                "cccd9bb7eab5a14ee23ddc9530406e8f9b879a72442a4d19dbde690595569a53")
    prepared=read(repo/"benchmark_tools/results/orthobench_factorial_prepared_20260916.json",
                  "5c325f4d77865e0c7571fe4bb4df0be0959977a49e1d22169f192518700f9382")
    source=record(repo/"orthohmm/phylogeny_pipeline.py")
    if source["sha256"]!="44e00316546b5df78354badd2b1a6bb595b685e98b86f686f1a32543d7f15e4f":
        raise ValueError("Membership consumer differs from reviewed frozen source")
    records.append(source)
    arm=prepared["candidate_arms"]["p1_c1"]
    old=arm["membership_constraints"]
    check(old)
    records.extend([old,arm["candidate_partition"]])
    historical=json.loads(Path(old["path"]).read_text())
    rows=[]
    for row in parent["rows"]:
        current=next(r for r in row["readback"]["original_records"] if r["path"].endswith("/phylogeny_candidate_merges.json"))
        for item in (current,row["prediction"]):
            check(item)
            records.append(item)
        if row["prediction"]["sha256"]!=arm["candidate_partition"]["sha256"]:
            raise ValueError("Candidate file differs before constraint comparison")
        events=json.loads(Path(current["path"]).read_text())
        rows.append(dict(label=row["label"],trace=current,
            raw_trace_byte_equal=current["sha256"]==old["sha256"],**compare_constraints(historical,events)))
    for item in records: check(item)
    return dict(status="canonical_membership_constraints_compared",source=record(__file__),checked_records=records,
        consumer="apply_satellite_membership_constraints",consumed_constraint_fields=["source_genes","target_genes"],
        rows=rows,all_semantic_multisets_equal=all(r["semantic_multiset_equal"] for r in rows),
        outcome_ordering_preserved=all(r["semantic_sequence_equal"] for r in rows),
        accuracy_evaluated=False,phylogeny_run=False,
        limitations=["Reviewed consumer source ignores trace scores, margins, iteration and cluster labels.",
            "With identical family outcomes, directed constraint multiset equality preserves membership and pair decisions.",
            "Constraint order can change root-group ordering/labels; byte equality of downstream outputs is not claimed.",
            "Does not establish identical newly inferred trees, family outcomes, final scores or cross-runtime behavior."])


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    if a.output.exists(): raise FileExistsError(a.output)
    write_json(a.output,audit(a.repo.resolve()))
