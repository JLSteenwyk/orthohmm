"""Fixed-seed candidate-only score/order factorial with a self-hit control."""

import argparse
import json
import os
from pathlib import Path
import shutil
import sys
import time

LABELS = ("historical_order_historical_scores", "historical_order_fresh_scores",
          "fresh_order_historical_scores", "fresh_order_fresh_scores", "fresh_full_self_control")


def prepare(repo, output):
    sys.path.insert(0, str(repo))
    from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
    from benchmark_tools.probe_installed_ob_clustering import write_json
    if output.exists():
        raise FileExistsError(output)
    parent = repo / "benchmark_tools/results/ob_dependency_replay_readback_22320.json"
    if record(parent)["sha256"] != "973d1bf10bb0e4e0028e0762cbee36def2c57e4897aa1308c25e231b87f55df4":
        raise ValueError("Changed dependency evidence")
    previous = json.loads(parent.read_text())
    pinned = {r["path"]: r for r in previous["checked_records"]}
    search_path = repo / "benchmark_tools/results/installed_ob_search_comparison_20260926.json"
    if record(search_path)["sha256"] != "247133956257e25cfbf5b2b79a19544ea382de07f9a7393f275d707ee8a658c5":
        raise ValueError("Changed historical cache evidence")
    search = json.loads(search_path.read_text())
    for item in search["checked_records"]:
        if item["path"] in pinned and pinned[item["path"]] != item:
            raise ValueError("Conflicting retained input identities")
        pinned[item["path"]] = item
    prior = repo / "benchmarks/work/ob_dependency_replay_v2_20260926"
    runtime_path = prior / "leiden011/runtime.json"
    check(pinned[str(runtime_path)])
    runtime = json.loads(runtime_path.read_text())
    checkpoint = repo / "benchmarks/work/publication_installed_orthobench_20260926/inference/orthohmm_working_res/high_sensitivity_checkpoint"
    seed = prior / "leiden011/replay/orthogroups_profiles_refined.txt"
    cache = repo / "benchmarks/results/hits_BLOSUM62_mc100.pkl"
    records = [record(parent), record(search_path), pinned[str(runtime_path)], *runtime["sources"], *runtime["native"],
               pinned[str(seed)], pinned[str(cache)],
               pinned[str(prior / "leiden011_venv/lib/python3.10/site-packages/frozen_diagnostic.pth")]]
    records.extend(pinned[str(p)] for p in sorted(checkpoint.iterdir()))
    for name in ("probe_ob_candidate_order_scores.py", "probe_installed_ob_graph.py",
                 "replay_high_sensitivity.py", "prepare_ob_candidate_neighborhood.py",
                 "audit_ob_dependency_replay.py", "audit_candidate_arm.py",
                 "audit_historical_profile_ablation.py", "compare_installed_ob_search.py",
                 "audit_installed_orthobench.py", "probe_installed_ob_clustering.py"):
        records.append(record(repo / "benchmark_tools" / name))
    for item in records:
        check(item)
    output.mkdir(parents=True)
    write_json(output / "plan.json", dict(repo=str(repo), output=str(output), labels=list(LABELS),
        checked_records=records, runtime=runtime, source=record(__file__), seed=record(seed), cache=record(cache),
        checkpoint=str(checkpoint), attempts=1, scoring=False, phylogeny=False,
        env=dict(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", PYTHONHASHSEED="0")))


def factorial_arrays(n, old, fresh):
    import numpy as np
    from benchmark_tools.probe_installed_ob_graph import align_scores
    q, t, s = fresh
    keep = q != t
    nonself = q[keep], t[keep], s[keep]
    old_in_new, new_in_old = align_scores(n, *old, *nonself)
    return dict(zip(LABELS, (old, (old[0], old[1], new_in_old),
                (nonself[0], nonself[1], old_in_new), nonself, fresh))), int(np.count_nonzero(~keep))


def run(plan_path, sha):
    # Import the frozen installation before making benchmark helpers visible.
    import orthohmm
    import numpy as np
    import importlib.metadata as metadata
    raw = plan_path.read_bytes()
    import hashlib
    if hashlib.sha256(raw).hexdigest() != sha:
        raise ValueError("Changed or unpinned plan")
    plan = json.loads(raw)
    sys.path.insert(0, plan["repo"])
    from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
    from benchmark_tools.probe_installed_ob_clustering import write_json
    from benchmark_tools.replay_high_sensitivity import load_replay_input
    from benchmark_tools.audit_ob_dependency_replay import candidate_readback
    from benchmark_tools.compare_installed_ob_search import partition
    from benchmark_tools.audit_installed_orthobench import compare_partitions
    from orthohmm.orthohmm import _expand_phylogeny_candidates
    from orthohmm.accuracy import load_accuracy_checkpoint
    import igraph._igraph
    import leidenalg._c_leiden
    output = Path(plan["output"])
    if (output / "started.json").exists():
        raise FileExistsError("Existing attempt; no retry/resume")
    if plan["labels"] != list(LABELS) or any(os.environ.get(k) != v for k, v in plan["env"].items()):
        raise ValueError("Changed arms or thread/hash settings")
    expected = plan["runtime"]
    observed = dict(versions={n:metadata.version(n) for n in ("numpy", "igraph", "leidenalg")},
        sources=[record(p) for p in sorted(Path(orthohmm.__file__).parent.rglob("*.py"))],
        native=[record(igraph._igraph.__file__), record(leidenalg._c_leiden.__file__)])
    if sys.executable != expected["executable"] or any(observed[k] != expected[k] for k in observed):
        raise ValueError("Runtime differs from fixed dependency arm")
    for item in plan["checked_records"]:
        check(item)
    write_json(output / "started.json", dict(plan=record(plan_path), runtime=observed,
               command=sys.argv, job_id=os.environ.get("SLURM_JOB_ID")))
    started = time.monotonic()
    try:
        # load_replay_input is permitted only after the admitted pickle hash check.
        names, species, q, t, s, _ = load_replay_input(pickle_path=Path(plan["cache"]["path"]))
        new_names, new_species, nq, nt, ns = load_accuracy_checkpoint(Path(plan["checkpoint"]))
        if names != new_names:
            raise ValueError("Gene indexing differs")
        classes = set(zip(species.tolist(), new_species.tolist()))
        if len(classes) != len(set(species)) or len(classes) != len(set(new_species)):
            raise ValueError("Species memberships differ")
        arrays, self_hits = factorial_arrays(len(names), (q,t,s), (nq,nt,ns))
        universe = set(names)
        partition(Path(plan["seed"]["path"]), universe)
        rows, groups = [], {}
        for label in LABELS:
            root = output / label
            target = root / "replay"
            work = target / "orthohmm_working_res"
            work.mkdir(parents=True)
            shutil.copyfile(plan["seed"]["path"], work / "orthohmm_edges_clustered.txt")
            tick = time.monotonic()
            details = _expand_phylogeny_candidates(str(target), names, species, arrays[label], profile="satellite_v2")
            details.pop("_membership_constraints", None)
            evidence = candidate_readback(root, dict(candidates=details), plan["seed"], universe)
            pred = record(work / "phylogeny_candidate_superfamilies.txt")
            groups[label] = partition(Path(pred["path"]), universe)
            row = dict(label=label, prediction=pred, candidate_summary=details, readback=evidence,
                       wall_seconds_descriptive=time.monotonic()-tick)
            write_json(root / "result.json", row)
            rows.append(row)
        comparisons = {a+"__vs__"+b:compare_partitions(groups[a],groups[b])
                       for i,a in enumerate(LABELS) for b in LABELS[i+1:]}
        for item in plan["checked_records"]:
            check(item)
        write_json(output / "report.json", dict(status="candidate_score_order_factorial_complete",
            plan=record(plan_path), source=record(__file__), rows=rows, comparisons=comparisons,
            genes=len(names), nonself_hits=len(q), removed_self_hits=self_hits,
            accuracy_evaluated=False, phylogeny_run=False, wall_seconds_descriptive=time.monotonic()-started,
            limitations=["One execution per arm in a fixed runtime; no general determinism or performance ranking.",
                         "Candidate-only diagnostic; not a final-phylogeny score comparison."]))
    except BaseException as error:
        write_json(output / "failure.json", dict(error_type=type(error).__name__, error=str(error),
                   wall_seconds=time.monotonic()-started, retry=False))
        raise


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path)
    p.add_argument("--prepare", type=Path)
    p.add_argument("--run", type=Path)
    p.add_argument("--plan-sha256")
    a = p.parse_args()
    if a.prepare:
        prepare(a.repo.resolve(), a.prepare.resolve())
    elif a.run:
        run(a.run.resolve(), a.plan_sha256)
    else:
        p.error("Require --prepare or --run")
