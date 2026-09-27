"""Diagnostic-only canonical input adapter around the retained factorial driver."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np


def array_identity(array):
    value = np.ascontiguousarray(array)
    return dict(shape=list(value.shape), dtype=value.dtype.str,
                sha256=hashlib.sha256(memoryview(value).cast("B")).hexdigest())


def canonical_factorial(names, arms):
    from benchmark_tools.candidate_hit_order_policy import canonical_hit_order_v1, POLICY
    from benchmark_tools.probe_ob_candidate_order_scores import LABELS
    if tuple(arms) != LABELS:
        raise ValueError("Wrong factorial arm inventory")
    result = {label:canonical_hit_order_v1(names,*arms[label]) for label in LABELS}
    equal = lambda a,b: all(np.array_equal(x,y) for x,y in zip(a,b))
    checks = dict(historical_scores_order_invariant=equal(result[LABELS[0]],result[LABELS[2]]),
                  fresh_scores_order_invariant=equal(result[LABELS[1]],result[LABELS[3]]))
    full = result[LABELS[4]]
    keep = full[0] != full[1]
    checks["self_control_nonself_equal"] = equal(result[LABELS[3]],tuple(a[keep] for a in full))
    if not all(checks.values()):
        raise ValueError("Canonical factorial input identities disagree")
    evidence = dict(policy=POLICY, checks=checks,
        canonical_arrays={label:[array_identity(a) for a in result[label]] for label in LABELS},
        self_hits=int(np.count_nonzero(~keep)), scores_rounded=False, scoring=False)
    return result,evidence


def prepare(repo, directory, fixture=False):
    sys.path.insert(0,str(repo))
    from benchmark_tools.prepare_ob_candidate_neighborhood import record,check
    from benchmark_tools.probe_installed_ob_clustering import write_json
    from benchmark_tools.candidate_hit_order_policy import POLICY
    if directory.exists():
        raise FileExistsError(directory)
    parent = repo / ("benchmarks/work/ob_candidate_order_smoke_20260926/plan.json" if fixture else
                     "benchmarks/work/ob_candidate_order_scores_v2_20260926/plan.json")
    expected = ("af1722e1a1fb2ca667289b6eb9859a469d15dfaa554a18329a70a653c55f5b14" if fixture else
                "cc3bdcef846c68b5d6a2e524a976d68325cfd17e184fdb2934eac618011a9ae8")
    if record(parent)["sha256"] != expected:
        raise ValueError("Changed parent plan")
    plan = json.loads(parent.read_text())
    driver = repo / "benchmark_tools/probe_ob_candidate_order_scores.py"
    policy = repo / "benchmark_tools/candidate_hit_order_policy.py"
    if (record(driver)["sha256"] != "e23aa1e8de1ae368624a1efd036924bb42d1d634092922e7eab782cddbc02ff9"
            or record(policy)["sha256"] != "cc70b884e40fe73a3c25ef9ae60a2133508127323a11382a8c610e0d192ec811"):
        raise ValueError("Changed driver or experimental policy")
    for item in plan["checked_records"]:
        check(item)
    plan.update(output=str(directory), parent_plan=record(parent), wrapper=record(__file__),
                intervention=POLICY, fixture=fixture)
    plan["checked_records"].extend([record(parent),record(__file__),record(policy)])
    directory.mkdir(parents=True)
    write_json(directory / "plan.json",plan)


def run(path, sha):
    # Preserve installed-package selection before exposing any audit modules.
    import orthohmm
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != sha:
        raise ValueError("Changed or unpinned canonical plan")
    plan = json.loads(raw)
    sys.path.insert(0,plan["repo"])
    from benchmark_tools.prepare_ob_candidate_neighborhood import record,check
    from benchmark_tools.probe_installed_ob_clustering import write_json
    from benchmark_tools.candidate_hit_order_policy import POLICY
    from benchmark_tools import probe_ob_candidate_order_scores as base
    if plan["wrapper"] != record(__file__) or plan["intervention"] != POLICY:
        raise ValueError("Wrong adapter or policy")
    for item in plan["checked_records"]:
        check(item)
    directory = Path(plan["output"])
    names = (Path(plan["checkpoint"]) / "gene_names.txt").read_text().splitlines()
    original = base.factorial_arrays
    def adapted(n, old, fresh):
        if n != len(names):
            raise ValueError("Changed gene indexing")
        arrays,self_hits = original(n,old,fresh)
        canonical,evidence = canonical_factorial(names,arrays)
        if evidence["self_hits"] != self_hits:
            raise ValueError("Self-hit counts disagree")
        write_json(directory / "canonical_inputs.json",evidence)
        return canonical,self_hits
    # Only the diagnostic input constructor changes; scientific functions do not.
    base.factorial_arrays = adapted
    try:
        base.run(path,sha)
    finally:
        base.factorial_arrays = original
    write_json(directory / "wrapper_complete.json",dict(status="canonical_candidate_native_complete_pending_readback",
               plan=record(path), wrapper=record(__file__), native_report=record(directory / "report.json"),
               canonical_inputs=record(directory / "canonical_inputs.json"), production_changed=False,
               installed_package=str(Path(orthohmm.__file__).resolve())))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo",type=Path)
    p.add_argument("--prepare",type=Path)
    p.add_argument("--fixture",action="store_true")
    p.add_argument("--run",type=Path)
    p.add_argument("--plan-sha256")
    a = p.parse_args()
    if a.prepare:
        prepare(a.repo.resolve(),a.prepare.resolve(),a.fixture)
    elif a.run:
        run(a.run.resolve(),a.plan_sha256)
    else:
        p.error("Require --prepare or --run")
