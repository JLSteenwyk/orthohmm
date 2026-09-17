"""Freeze all simulation mode controls, retaining unavailable inferred baselines."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_simulation_tree_controls import CONDITIONS
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.run_simulation_tree_mode_control import METHOD_SHA, RESULT_SHA, METHODS
from benchmark_tools.admit_simulation_mode_pilot import PILOT_SHA


def panel_rows(manifest, evidence):
    expected = {(condition, seed) for condition in CONDITIONS for seed in range(20261101, 20261111)}
    datasets = manifest["datasets"]
    if len(datasets) != 70 or {(d["condition"], d["seed"]) for d in datasets} != expected:
        raise ValueError("Incomplete or duplicate simulation inventory")
    rows = []
    for dataset in sorted(datasets, key=lambda d: (d["seed"], CONDITIONS.index(d["condition"]))):
        label = f"{dataset['condition']}_{dataset['seed']}"
        if dataset["label"] != label:
            raise ValueError("Dataset label differs from condition and seed")
        available, unavailable = [], {}
        for method in METHODS:
            matches = [r for r in evidence["records"] if r["condition"] == dataset["condition"]
                       and r["seed"] == dataset["seed"] and r["method"] == method]
            if len(matches) != 1 or matches[0]["status"] not in {"complete", "failed", "inapplicable"}:
                raise ValueError("Missing, duplicate or nonterminal baseline outcome")
            source = matches[0]
            if source["status"] == "complete":
                available.append(method)
            else:
                unavailable[method] = {key: source.get(key) for key in ("status", "failure_stage", "reason")}
        reuse = label == "baseline_20261101"
        if reuse and tuple(available) != METHODS:
            raise ValueError("Pilot method inventory differs")
        rows.append({"index": len(rows), "label": label, "condition": dataset["condition"], "seed": dataset["seed"],
                     "methods": available, "unavailable": unavailable, "reuse_pilot": reuse})
    return rows


def prepare(root, output):
    if output.exists():
        raise FileExistsError(output)
    results = root / "benchmark_tools/results"
    method_path = results / "publication_variable_native_methods_20260916.json"
    evidence_path = results / "simulation_variable_native_results_20260916.json"
    admission_path = root / "benchmarks/results/simulation_mode_pilot_admission_v1/results.json"
    admission = json.loads(admission_path.read_text())
    if (admission["status"] != "pilot_equivalence_verified" or admission["pilot_report"]["sha256"] != PILOT_SHA
            or admission["accuracy_evaluated"] is not False):
        raise ValueError("Pilot not independently admitted")
    check(admission["source"])
    check(admission["pilot_report"])
    for item in admission["executor_sources"]:
        check(item)
    rows = panel_rows(read_frozen(method_path, METHOD_SHA), read_frozen(evidence_path, RESULT_SHA))
    result = {"status": "mode_panel_frozen_not_executed", "accuracy_evaluated": False, "source": record(__file__),
              "method_manifest": record(method_path), "baseline_results": record(evidence_path),
              "pilot_admission": record(admission_path), "rows": rows,
              "new_method_runs": sum(len(r["methods"]) for r in rows if not r["reuse_pilot"]),
              "unavailable_method_controls": sum(len(r["unavailable"]) for r in rows),
              "reused_method_runs": sum(len(r["methods"]) for r in rows if r["reuse_pilot"]),
              "scope": "All70datasets; absent inferred baselines are unavailable controls, not removed from the future oracle panel"}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.output.resolve())
