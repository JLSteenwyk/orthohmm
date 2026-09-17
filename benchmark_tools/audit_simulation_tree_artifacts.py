"""Compare retained upstream artifacts across admitted tree arms without scoring truth."""

import argparse
from collections import Counter
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_simulation_mode_pilot import retained_equivalence
from benchmark_tools.audit_mode_partitions import ADMISSION_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.run_simulation_tree_experiment import PORTABLE_SHA
from benchmark_tools.run_simulation_tree_mode_control import METHODS, compare_inventory
from benchmark_tools.verify_ygob_validation import require_completed_job

CONTRASTS = (("generating", "inferred"), ("nni1", "generating"), ("nni2", "generating"))
VARIANTS = {"generating", "nni1", "nni2"}


def index_panel(admission, trees):
    if admission["status"] != "tree_panel_verified_unscored" or admission["accuracy_evaluated"] is not False:
        raise ValueError("Require independent complete-panel admission before artifact contrasts")
    expected = {(t["condition"], t["seed"], m, t["variant"]) for t in trees for m in METHODS}
    rows = admission["records"]
    indexed = {(r["condition"], r["seed"], r["method"], r["variant"]): r for r in rows}
    if len(expected) != 420 or len(rows) != 420 or set(indexed) != expected:
        raise ValueError("Incomplete or duplicated method/tree inventory")
    if any(r["status"] not in {"admitted", "failed", "supplied_tree_not_retained"} for r in rows):
        raise ValueError("Unknown admitted panel outcome")
    return indexed


def contrast(method, target, reference):
    statuses = {"target": target["status"], "reference": reference["status"]}
    if any(r["status"] != "admitted" for r in (target, reference)):
        return {"status": "unavailable", "arm_statuses": statuses,
                "reason": "At least one arm lacks admitted output; not interpreted as upstream equivalence"}
    before, after = reference["retained_artifacts"], target["retained_artifacts"]
    if not before or not after:
        raise ValueError("Empty upstream artifact inventory")
    for item in [*before.values(), *after.values()]:
        check(item)
    comparison = compare_inventory(before, after)
    equivalent = retained_equivalence(method, before, after)
    return {"status": "retained_upstream_equivalent" if equivalent else "retained_upstream_different",
            "arm_statuses": statuses, "byte_comparison": comparison,
            "representation_exception_used": equivalent and not comparison["identical"],
            "scope": "Retained candidate/graph, alignment and raw gene-tree files only; not all internal computation"}


def inferred_inventory(row):
    if row["status"] == "unavailable":
        return {"status": "unavailable", "baseline_failure": row["baseline_failure"]}
    if row["status"] != "equivalent":
        raise ValueError("Original mode-control admission is not equivalent")
    if row.get("reused_pilot"):
        check(row["admission"])
        link = json.loads(Path(row["admission"]["path"]).read_text())["pilot_report"]
    else:
        link = row["native_report"]
    check(link)
    native = json.loads(Path(link["path"]).read_text())
    baseline = native["provenance"]["baseline"][row["method"]]
    if baseline["admission"]["status"] != "admitted":
        raise ValueError("Original baseline not natively admitted")
    for item in baseline["retained_artifacts"].values():
        check(item)
    return {"status": "admitted", "retained_artifacts": baseline["retained_artifacts"], "native_report": link}


def audit(root, output, digest):
    if output.exists():
        raise FileExistsError(output)
    path = root / "benchmarks/results/simulation_tree_panel_admission_v1/results.json"
    admission = read_frozen(path, digest)
    accounting = subprocess.check_output(["sacct", "-j", "21435", "--parsable2",
                                         "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = require_completed_job(accounting, 21435)
    executor = root / "benchmarks/work/publication_simulation_tree_admission_v1"
    if admission["source"] != record(executor / "benchmark_tools/admit_simulation_tree_panel.py"):
        raise ValueError("Unexpected tree-panel auditor source")
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != "c13b17002b22e5918d2e26883da2b3ab03efb5d8":
        raise ValueError("Tree-panel auditor revision changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    for item in [admission["source"], *admission["helper_sources"], *admission["gates"].values()]:
        check(item)
    portable_path = root / "benchmark_tools/results/simulation_portable_trees_prepared_20260917.json"
    portable = read_frozen(portable_path, PORTABLE_SHA)
    indexed = index_panel(admission, portable["trees"])
    mode_path = root / "benchmarks/results/simulation_mode_panel_admission_v1/results.json"
    mode = read_frozen(mode_path, ADMISSION_SHA)
    baselines = {(r["condition"], r["seed"], r["method"]): r for r in mode["records"]}
    datasets = sorted({(t["condition"], t["seed"]) for t in portable["trees"]})
    expected_baselines = {(c, s, m) for c, s in datasets for m in METHODS}
    if len(mode["records"]) != 140 or set(baselines) != expected_baselines:
        raise ValueError("Incomplete original-baseline inventory")
    contrasts = []
    for condition, seed in datasets:
        for method in METHODS:
            arms = {v: indexed[condition, seed, method, v] for v in VARIANTS}
            arms["inferred"] = inferred_inventory(baselines[condition, seed, method])
            for arm in arms.values():
                for key in ("native_report", "execution", "preflight"):
                    if key in arm:
                        check(arm[key])
            for target, reference in CONTRASTS:
                contrasts.append({"condition": condition, "seed": seed, "method": method,
                                  "target": target, "reference": reference,
                                  **contrast(method, arms[target], arms[reference])})
    if len(contrasts) != 420:
        raise ValueError("Incomplete planned artifact contrasts")
    read_frozen(path, digest)
    read_frozen(mode_path, ADMISSION_SHA)
    result = {"status": "tree_artifact_contrasts_checked", "accuracy_evaluated": False,
              "source": record(__file__), "admission": record(path), "mode_admission": record(mode_path),
              "portable_manifest": record(portable_path), "admission_scheduler": scheduler, "contrasts": contrasts,
              "counts": dict(Counter(r["status"] for r in contrasts)),
              "helper_sources": [record(Path(__file__).with_name(n)) for n in
                                 ("admit_simulation_mode_pilot.py", "run_simulation_tree_mode_control.py")],
              "limitations": ["Artifact differences are not an exclusion rule for prespecified accuracy contrasts.",
                              "Equivalence covers retained upstream artifacts, not a proof of a topology-only mechanism.",
                              "Unavailable inferred baselines and failed supplied arms remain explicit."]}
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--admission-sha256", required=True)
    args = parser.parse_args()
    audit(args.root.resolve(), args.output.resolve(), args.admission_sha256)
