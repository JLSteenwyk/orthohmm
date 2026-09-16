"""Verify supplied-tree execution, native semantics and baseline equivalence."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import Phylo
from benchmark_tools.audit_historical_profile_ablation import verify_file
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.prepare_species_tree_robustness import topology
from benchmark_tools.recover_factorial_postflight import PREPARED_HASH, ENVIRONMENT_HASH
from benchmark_tools.run_orthobench_factorial_cell import select_cell, verify_prepared
from benchmark_tools.run_simulation_methods import read_frozen, verify_environment
from benchmark_tools.run_species_tree_control import TREE_MANIFEST_SHA, control_cell
from benchmark_tools.score_ygob_groups import read_predictions
from benchmark_tools.validate_factorial_native import validate as validate_baseline, validate_native_cell
from benchmark_tools.validate_simulation_outputs import verify_process
from benchmark_tools.verify_ygob_validation import require_completed_job


def check_status(status, postflight, cell):
    if (status["status"] != "finished_pending_native_validation" or status["failed_methods"] != []
            or status["accuracy_evaluated"] is not False or status["native_outputs_validated"] is not False):
        raise ValueError("Control is not successfully completed and unscored")
    if status["dataset"] != cell["label"] or set(status["methods"]) != {cell["label"]}:
        raise ValueError("Control method inventory differs")
    if status["provenance"]["job_id"] != "21298" or status["provenance"]["cell"] != cell:
        raise ValueError("Wrong control scheduler or cell provenance")
    if postflight != {"status": "complete_pending_native_equivalence", "accuracy_evaluated": False,
                      "baseline_source_unchanged": True, "cell": cell}:
        raise ValueError("Control postflight did not pass")


def validate(root):
    root = root.resolve()
    if Path.cwd().resolve() != root:
        raise ValueError("Run from the original repository verification directory")
    accounting = subprocess.check_output(["sacct", "-j", "21298", "--parsable2",
                                         "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = require_completed_job(accounting, 21298)
    results = root / "benchmark_tools/results"
    prepared_path = results / "orthobench_factorial_prepared_20260916.json"
    env_path = results / "publication_variable_native_methods_20260916.json"
    trees_path = results / "ob_species_tree_robustness_prepared_20260916.json"
    prepared = read_frozen(prepared_path, PREPARED_HASH)
    environment = read_frozen(env_path, ENVIRONMENT_HASH)
    trees = read_frozen(trees_path, TREE_MANIFEST_SHA)
    baseline = validate_baseline(root, 3)
    original, _, launcher = select_cell(prepared, 3)
    output = root / "benchmarks/results/ob_supplied_tree_control_v1"
    control = next(t for t in trees["variants"] if t["label"] == "supplied_control")
    cell = control_cell(original, control, output)
    executor = root / "benchmarks/work/publication_ob_tree_control_v1"
    revision = subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip()
    expected_revision = subprocess.check_output(["git", "-C", str(root), "rev-parse", "d9ea049^{commit}"], text=True).strip()
    if revision != expected_revision:
        raise ValueError("Control executor revision changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    records = [file_provenance(output / name) for name in ("preflight.json", "postflight.json", "execution/status.json")]
    preflight, postflight, status = [json.loads(Path(record["path"]).read_text()) for record in records]
    check_status(status, postflight, cell)
    if status["provenance"] != preflight:
        raise ValueError("Control execution provenance differs from preflight")
    for key, path in {"source": executor / "benchmark_tools/run_species_tree_control.py",
                      "tree_manifest": trees_path, "prepared_manifest": prepared_path, "environment_manifest": env_path}.items():
        if preflight[key] != file_provenance(path):
            raise ValueError("Control input or executor provenance differs: " + key)
    for key in ("native_manifest", "native_metrics", "species_tree", "partition"):
        if baseline[key] != preflight["baseline_native_admission"][key]:
            raise ValueError("Original native checkpoint source changed")
    verify_prepared(prepared, original, launcher, 21161)
    verify_environment(environment)
    verify_file(Path(control["tree"]["path"]), control["tree"])
    method = {"argv": cell["argv"], "output": str(output / "output"), "metrics": str(output / "metrics.json")}
    artifacts = verify_process(method, status["methods"][cell["label"]])
    native = validate_native_cell(prepared, environment, cell, output, launcher,
                                  {"scheduler": scheduler, "verified_artifacts": len(artifacts)})
    source_tree = Phylo.read(trees["source_tree"]["path"], "newick")
    supplied = Phylo.read(native["species_tree"]["path"], "newick")
    if topology(source_tree) != topology(supplied):
        raise ValueError("Unchanged supplied control changed rooted topology")
    before = {frozenset(g) for g in read_predictions(Path(original["prediction"]), "root_hogs").values()}
    after = {frozenset(g) for g in read_predictions(Path(cell["prediction"]), "root_hogs").values()}
    for record in records:
        verify_file(Path(record["path"]), record)
    return {"status": "equivalent" if before == after else "not_equivalent", "accuracy_evaluated": False,
            "native_validation": native, "baseline_validation": baseline, "execution": records,
            "expected_only_groups": len(before - after), "observed_only_groups": len(after - before),
            "expected_groups": len(before), "observed_groups": len(after),
            "verifier": file_provenance(Path(__file__))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = validate(args.root)
    with args.output.open("x") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    if report["status"] != "equivalent":
        raise SystemExit("Supplied-tree control is not equivalent; report preserved")


if __name__ == "__main__":
    main()
