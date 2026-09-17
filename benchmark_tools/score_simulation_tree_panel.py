"""Score complete admitted tree arms against frozen evolutionary truth, retaining every failure."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_simulation_tree_artifacts import index_panel
from benchmark_tools.assemble_simulation_results import input_universe, verify_scoring_dependencies
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen, verify_inputs
from benchmark_tools.run_simulation_tree_mode_control import METHOD_SHA, RESULT_SHA, METHODS
from benchmark_tools.run_simulation_tree_experiment import PORTABLE_SHA
from benchmark_tools.simulation_method_outputs import load_predictions
from benchmark_tools.simulation_conditions import score_pairs
from benchmark_tools.summarize_simulation_tree_panel import summarize, index_records
from benchmark_tools.verify_ygob_validation import require_completed_job


def check_artifact_gate(artifacts, admission_record, indexed):
    if (artifacts["status"] != "tree_artifact_contrasts_checked" or artifacts["accuracy_evaluated"] is not False
            or artifacts["admission"] != admission_record):
        raise ValueError("Artifact audit does not reference the admitted panel")
    expected = {(c, s, m, a, b) for c, s, m, v in indexed for a, b in
                (("generating", "inferred"), ("nni1", "generating"), ("nni2", "generating"))}
    rows = artifacts["contrasts"]
    if (len(rows) != 420 or {(r["condition"], r["seed"], r["method"], r["target"], r["reference"]) for r in rows} != expected
            or any(r["status"] not in {"retained_upstream_equivalent", "retained_upstream_different", "unavailable"} for r in rows)):
        raise ValueError("Incomplete artifact interpretation inventory")


def score(root, output, admission_sha, artifact_sha):
    if output.exists():
        raise FileExistsError(output)
    path = root / "benchmarks/results/simulation_tree_panel_admission_v1/results.json"
    artifact_path = root / "benchmarks/results/simulation_tree_artifacts_v1.json"
    admission, artifacts = read_frozen(path, admission_sha), read_frozen(artifact_path, artifact_sha)
    accounting = subprocess.check_output(["sacct", "-j", "21435", "--parsable2", "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = require_completed_job(accounting, 21435)
    results = root / "benchmark_tools/results"
    portable = read_frozen(results / "simulation_portable_trees_prepared_20260917.json", PORTABLE_SHA)
    indexed = index_panel(admission, portable["trees"])
    check_artifact_gate(artifacts, record(path), indexed)
    for item in [admission["source"], *admission["helper_sources"], *admission["gates"].values(),
                 artifacts["source"], *artifacts["helper_sources"]]:
        check(item)
    manifest = read_frozen(results / "publication_variable_native_methods_20260916.json", METHOD_SHA)
    original = read_frozen(results / "simulation_variable_native_results_20260916.json", RESULT_SHA)
    baseline = {(r["condition"], r["seed"], r["method"]): r for r in original["records"] if r["method"] in METHODS}
    if len(baseline) != 140:
        raise ValueError("Incomplete original baseline inventory")
    gen = manifest["generation_manifest"]
    generation = read_frozen(Path(gen["absolute_path"]), gen["sha256"])
    verify_scoring_dependencies(manifest)
    output.mkdir(parents=True)
    report = {"status": "scoring", "source": record(__file__), "admission": record(path), "artifact_audit": record(artifact_path),
              "admission_scheduler": scheduler, "records": [],
              "helper_sources": [record(Path(__file__).with_name(n)) for n in
                 ("simulation_conditions.py", "simulation_method_outputs.py", "summarize_simulation_tree_panel.py",
                  "assemble_simulation_results.py", "run_simulation_methods.py", "audit_simulation_tree_artifacts.py")]}
    try:
        for dataset in manifest["datasets"]:
            verified = verify_inputs(dataset, generation, root / "benchmarks/work/publication_variable_simulation_panel_v2", gen["sha256"])
            truth = json.loads(Path(dataset["truth"]).read_text())
            owners, species = input_universe(verified["inputs"], truth)
            common = {"condition": dataset["condition"], "seed": dataset["seed"], "truth_sha256": verified["truth"]["sha256"]}
            for method in METHODS:
                original_row = baseline[dataset["condition"], dataset["seed"], method]
                rows = [("inferred", original_row)] + [(arm, indexed[dataset["condition"], dataset["seed"], method, arm])
                                                       for arm in ("generating", "nni1", "nni2")]
                for arm, native in rows:
                    row = {**common, "method": method, "arm": arm}
                    if native["status"] not in {"complete", "admitted"}:
                        row.update(status="failed", reason=native.get("reason") or native.get("admission", {}).get("reason") or native["status"],
                                   native_outcome=native)
                    else:
                        if arm == "inferred":
                            directory = Path(dataset["methods"][method]["output"])
                            expected_files = native["prediction_artifacts"]
                        else:
                            check(native["native_report"])
                            execution = json.loads(Path(native["native_report"]["path"]).read_text())
                            directory = Path(execution["provenance"]["configured"]["methods"][method]["output"])
                            expected_files = native["prediction_files"]
                        predictions, files = load_predictions(method, directory, owners, species)
                        actual = [record(p) for p in files]
                        expected = [{"path": f.get("absolute_path", f["path"]), "bytes": f["bytes"], "sha256": f["sha256"]}
                                    for f in expected_files]
                        if actual != expected:
                            raise ValueError("Native prediction artifacts differ from admission")
                        computed = score_pairs(predictions, truth["ortholog_pairs"], owners)
                        if arm == "inferred" and (computed != native["score"] or native["truth_sha256"] != common["truth_sha256"]):
                            raise ValueError("Recomputed inferred baseline differs from frozen original")
                        row.update(status="complete", score=computed, prediction_files=actual)
                    report["records"].append(row)
        index_records(report["records"])
        read_frozen(path, admission_sha)
        read_frozen(artifact_path, artifact_sha)
        verify_scoring_dependencies(manifest)
        for item in [report["source"], *report["helper_sources"]]:
            check(item)
        report["status"] = "complete_tree_panel_scored"
        summary = summarize(report["records"])
        summary.update(source=record(Path(__file__).with_name("summarize_simulation_tree_panel.py")), artifact_audit=record(artifact_path))
        (output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    except Exception as error:
        report.update(status="scoring_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        (output / "results.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--admission-sha256", required=True)
    parser.add_argument("--artifact-sha256", required=True)
    args = parser.parse_args()
    score(args.root.resolve(), args.output.resolve(), args.admission_sha256, args.artifact_sha256)
