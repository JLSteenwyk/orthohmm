"""Fresh installed phylogeny stage from validated canonical candidates."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def replay_arguments(plan):
    return ["--fasta-directory",plan["inputs"],"--candidate-clusters",plan["candidate"]["path"],
        "--membership-constraints",plan["constraints"]["path"],"--output-directory",plan["inference"],
        "--json",str(Path(plan["output"])/"replay.json"),"--cpu",str(plan["cpu"]),
        "--species-tree-mode","infer","--species-tree-rooting","min_variance",
        "--root-rule","species_overlap","--pair-rule","positive_paralogy",
        "--aligner",plan["aligner"],"--tree-builder",plan["tree_builder"]]


def runtime_inventory():
    import importlib.metadata as metadata
    import orthohmm
    import dendropy
    from benchmark_tools.prepare_ob_candidate_neighborhood import record
    return dict(executable=sys.executable,python=sys.version,
        versions={n:metadata.version(n) for n in ("numpy","igraph","leidenalg","DendroPy")},
        scientific_sources=[record(p) for p in sorted(Path(orthohmm.__file__).parent.rglob("*.py"))],
        dendropy_sources=[record(p) for p in sorted(Path(dendropy.__file__).parent.rglob("*.py"))])


def prepare(repo,directory,fixture=False):
    sys.path.insert(0,str(repo))
    from benchmark_tools.prepare_ob_candidate_neighborhood import record,check
    from benchmark_tools.probe_installed_ob_clustering import write_json
    if directory.exists(): raise FileExistsError(directory)
    records=[]
    def read(path,sha):
        item=record(path)
        if item["sha256"]!=sha: raise ValueError("Changed parent evidence")
        records.append(item)
        return json.loads(Path(path).read_text())
    constraints_audit=read(repo/"benchmark_tools/results/canonical_membership_constraint_audit_20260926.json",
        "6751450a54a646db1573dcbdee1570e62fdc4c16d8a532935dc5bbb24bb5d29e")
    if not constraints_audit["all_semantic_multisets_equal"] or not constraints_audit["outcome_ordering_preserved"]:
        raise ValueError("Canonical constraint equivalence not established")
    canonical=read(repo/"benchmark_tools/results/ob_canonical_candidate_readback_22323.json",
        "cccd9bb7eab5a14ee23ddc9530406e8f9b879a72442a4d19dbde690595569a53")
    installed=read(repo/"benchmarks/work/publication_installed_orthobench_20260926/plan.json",
        "5fd8dc70c337951706fe9246c2bf08d1191196939da066c5b74bba0bd3dd2b12")
    records.extend(installed["checked_records"])
    row=next(r for r in canonical["rows"] if r["label"]=="fresh_full_self_control")
    candidate=row["prediction"]
    constraints=next(r for r in row["readback"]["original_records"] if r["path"].endswith("/phylogeny_candidate_merges.json"))
    inputs=repo/"benchmarks/work/publication_installed_orthobench_20260926/input"
    python=repo/"benchmarks/work/ob_dependency_replay_v2_20260926/leiden011_venv/bin/python"
    pth=python.parent.parent/"lib/python3.10/site-packages/frozen_diagnostic.pth"
    known={r["path"]:r for r in canonical["checked_records"]}
    records.extend([known[str(pth)],candidate,constraints])
    for name in ("run_canonical_ob_phylogeny.py","replay_phylogeny.py","orthobench_stage_diagnostics.py",
                 "prepare_ob_candidate_neighborhood.py","probe_installed_ob_clustering.py"):
        records.append(record(repo/"benchmark_tools"/name))
    if record(repo/"benchmark_tools/replay_phylogeny.py")["sha256"]!="e4de6f39c66c0f7e7f63229188a76ece44dbdc9b7771ab1c3b1465cec63aa539":
        raise ValueError("Changed phylogeny replay harness")
    environment={"HOME":os.environ["HOME"],"LANG":"C.UTF-8",**installed["environment_overrides"]}
    code=("import orthohmm,sys,json;sys.path.insert(0,"+repr(str(repo))+");"
          "from benchmark_tools.run_canonical_ob_phylogeny import runtime_inventory;print(json.dumps(runtime_inventory()))")
    runtime=json.loads(subprocess.check_output([str(python),"-I","-c",code],env=environment,text=True,timeout=60))
    for item in runtime["scientific_sources"]:
        if known.get(item["path"])!=item: raise ValueError("Installed scientific source differs")
    records.extend(runtime["scientific_sources"]+runtime["dendropy_sources"])
    command=installed["command"]
    aligner=command[command.index("--aligner")+1]
    tree_builder=command[command.index("--tree_builder")+1]
    if fixture:
        fixture_report=read(repo/"benchmark_tools/results/publication_frozen_phylogeny_smoke_20260926.json",
            "8d258122281bed99d26e4c66bdc03aa3b8bfa3367db3de1bdb5daf86c815e4e0")
        fixture_root=repo/"benchmarks/work/publication_frozen_phylogeny_smoke_20260926"
        inputs=fixture_root/"input"
        candidate=record(fixture_root/"inference/orthohmm_working_res/phylogeny_candidate_superfamilies.txt")
        if candidate["sha256"]!=fixture_report["result"]["manifest"]["input_cluster_sha256"]:
            raise ValueError("Fixture candidate identity differs")
        records.extend([candidate,*[record(p) for p in sorted(inputs.iterdir())]])
    for item in records: check(item)
    directory.mkdir(parents=True)
    if fixture:
        genes=Path(candidate["path"]).read_text().splitlines()[-1].split()
        write_json(directory/"fixture_constraints.json",[dict(source_genes=genes[:1],target_genes=genes[1:])])
        constraints=record(directory/"fixture_constraints.json")
        records.append(constraints)
    plan=dict(repo=str(repo),output=str(directory),inference=str(directory/"inference"),inputs=str(inputs),
        candidate=candidate,constraints=constraints,aligner=aligner,tree_builder=tree_builder,
        cpu=1 if fixture else 32,memory_gib=128,timeout_seconds=180 if fixture else 28800,
        attempts=1,fixture=fixture,python=str(python),environment=environment,runtime=runtime,
        checked_records=records,source=record(__file__),checkpoint_reuse=False,scoring=False)
    write_json(directory/"plan.json",plan)


def run(path,sha):
    import orthohmm
    raw=path.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=sha: raise ValueError("Changed or unpinned plan")
    plan=json.loads(raw)
    sys.path.insert(0,plan["repo"])
    from benchmark_tools.prepare_ob_candidate_neighborhood import record,check
    from benchmark_tools.probe_installed_ob_clustering import write_json
    from benchmark_tools import replay_phylogeny as replay
    if plan["source"]!=record(__file__): raise ValueError("Changed phylogeny driver")
    if plan["checkpoint_reuse"] is not False or plan["scoring"] is not False or plan["attempts"]!=1:
        raise ValueError("Wrong execution scope")
    directory=Path(plan["output"])
    if (directory/"started.json").exists() or Path(plan["inference"]).exists():
        raise FileExistsError("Existing attempt; no retry/resume")
    if any(os.environ.get(k)!=v for k,v in plan["environment"].items()):
        raise ValueError("Environment differs")
    runtime=runtime_inventory()
    if runtime!=plan["runtime"]: raise ValueError("Installed runtime differs")
    if Path(sys.modules["orthohmm.phylogeny_pipeline"].__file__).parent!=Path(orthohmm.__file__).parent:
        raise ValueError("Scientific import escaped installation")
    for item in plan["checked_records"]: check(item)
    args=replay_arguments(plan)
    write_json(directory/"started.json",dict(plan=record(path),runtime=runtime,arguments=args,command=sys.argv))
    start=time.monotonic()
    try:
        rc=replay.main(args)
        if rc: raise RuntimeError(f"Replay returned {rc}")
        for item in plan["checked_records"]: check(item)
        write_json(directory/"native_complete.json",dict(status="fresh_canonical_phylogeny_complete_pending_readback",
            plan=record(path),replay=record(directory/"replay.json"),wall_seconds_descriptive=time.monotonic()-start,
            accuracy_evaluated=False))
    except BaseException as error:
        write_json(directory/"failure.json",dict(error_type=type(error).__name__,error=str(error),retry=False,
                   wall_seconds=time.monotonic()-start))
        raise


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo",type=Path)
    p.add_argument("--prepare",type=Path)
    p.add_argument("--fixture",action="store_true")
    p.add_argument("--run",type=Path)
    p.add_argument("--plan-sha256")
    a=p.parse_args()
    if a.prepare: prepare(a.repo.resolve(),a.prepare.resolve(),a.fixture)
    elif a.run: run(a.run.resolve(),a.plan_sha256)
    else: p.error("Require --prepare or --run")
