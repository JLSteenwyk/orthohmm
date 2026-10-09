"""Bind retained controls and prepare one prospective fragment observation panel."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from benchmark_tools.assemble_simulation_results import (
    admit_method, input_universe, terminal_tasks, verify_execution_status, verify_scoring_dependencies,
)
from benchmark_tools.benchmark_production import file_record
from benchmark_tools.controlled_fragment_observations import (
    CONDITION, SEEDS, METHODS, INFERENCE_METHODS, observe, verify_observation, fresh_methods, score_strata,
)
from benchmark_tools.run_simulation_generation import verify_file
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.simulation_conditions import canonical_pairs, score_pairs
from benchmark_tools.simulation_method_outputs import load_predictions

PINS = {
    "simulation_variable_native_results_20260916.json": "cc99fc31c3433809098212d3c8dc12f4ad829b4cfeb66dabcb13aede4e95258f",
    "publication_variable_native_methods_20260916.json": "bf677728cf9312edd9755d9ac1ceefbce3abb8c3848d43413631dc00ca00eb2f",
    "publication_variable_simulation_manifest_20260916.json": "806aa1e5f6976c323ff2f6641dd88e25264e7417749e7d77565294666aee776b",
    "publication_variable_methods_20260916.json": "f68b0caf508fde8f40d311b5d303f7569a118bd15e8cf960d227a7d1cfee5984",
    "publication_native_runtime_20260916.json": "aebea83807356b02307473506fa30c2dbd2c511d7ba75a0655eb12180a474d74",
}
PROTOCOL = "CONTROLLED_FRAGMENT_OBSERVATION_PROTOCOL_20261009.md"
CORE = "7f3a9e40dd7e79f842cc2c11fb8b548f9a802806"
EXECUTORS = {
    "native": ("publication_native_simulation_v3", "b66225dc5fc355702575aebe34b644916894f236"),
    "comparator": ("publication_variable_methods_v2", "f5f4e1b662118e7624746546e81505ddeabf9bb8"),
}


def record(path):
    path = Path(path).resolve()
    return dict(file_record(path, path.parent), absolute_path=str(path))


def baseline_rows(report):
    rows = [row for row in report["records"] if row["condition"] == "baseline"]
    expected = {(seed, method) for seed in SEEDS for method in METHODS}
    keys = [(row["seed"], row["method"]) for row in rows]
    if len(keys) != len(expected) or set(keys) != expected or any(row["status"] != "complete" for row in rows):
        raise ValueError("All 40 prespecified baseline outcomes must be complete and unique")
    return {(row["seed"], row["method"]): row for row in rows}


def read_input_sequences(verified, truth, input_directory):
    owners, species = input_universe(verified["inputs"], truth)
    input_directory = Path(input_directory).resolve()
    paths = {Path(item["absolute_path"]).resolve() for item in verified["inputs"]}
    if paths != set(input_directory.iterdir()) or any(path.parent != input_directory or path.suffix != ".fasta" for path in paths):
        raise ValueError("Baseline input directory differs from its retained inventory")
    prepared = {item["path"]: {k: item[k] for k in ("path", "bytes", "sha256")} for item in truth["prepared_inputs"]}
    actual = {}
    sequences = {}
    for item in verified["inputs"]:
        path = Path(item["absolute_path"])
        relative = "input/" + path.name
        actual[relative] = {"path": relative, "bytes": item["bytes"], "sha256": item["sha256"]}
        for sequence in SeqIO.parse(path, "fasta"):
            sequences[sequence.id] = str(sequence.seq)
    if len(prepared) != len(truth["prepared_inputs"]) or prepared != actual:
        raise ValueError("Truth's prepared input pins differ from actual retained inputs")
    genes = [gene for members in truth["families"].values() for gene in members]
    if len(set(genes)) != len(genes) or set(genes) != set(owners):
        raise ValueError("Truth family membership is not an exact gene partition")
    _, duplicates = canonical_pairs(truth["ortholog_pairs"], owners)
    if duplicates:
        raise ValueError("Duplicate evolutionary truth pairs")
    return sequences, owners, species


def bind_baseline(dataset, rows, manifest, previous, tasks, executors):
    """Revalidate actual admission, command, scheduler and artifact bindings."""
    common_verified, predictions, bindings = None, {}, []
    for method in METHODS:
        row = rows[dataset["seed"], method]
        reused = method.startswith("orthofinder_")
        active = previous if reused else manifest
        manifest_hash = PINS["publication_variable_methods_20260916.json" if reused else "publication_variable_native_methods_20260916.json"]
        if row["inference_method_manifest_sha256"] != manifest_hash or row.get("reused_comparator") is not reused:
            raise ValueError("Baseline result has inconsistent method provenance")
        if row["scheduler"] != tasks["comparator" if reused else "native"]:
            raise ValueError("Baseline result differs from retained terminal accounting")
        expected_status = Path(dataset["methods"]["orthofinder_full" if reused else "orthohmm_high_sensitivity"]["output"]).parent / "execution/status.json"
        item = row["execution_evidence"]
        if Path(item["absolute_path"]).resolve() != expected_status.resolve():
            raise ValueError("Baseline execution evidence belongs to a different output")
        verify_file(expected_status, item)
        status = json.loads(expected_status.read_text())
        verified = status["verified_inputs"]
        if verified["status"] != "ready" or Path(verified["truth"]["absolute_path"]).resolve() != Path(dataset["truth"]).resolve():
            raise ValueError("Baseline truth path or readiness differs")
        verify_file(Path(dataset["truth"]), verified["truth"])
        if row["truth_sha256"] != verified["truth"]["sha256"]:
            raise ValueError("Baseline scored truth identity differs")
        if common_verified is not None and common_verified != verified:
            raise ValueError("Methods used different baseline inputs or truth")
        common_verified = verified
        verify_execution_status(status, dataset, row["scheduler"], verified, manifest_hash,
                                manifest["generation_manifest"]["sha256"], executors["comparator" if reused else "native"])
        admission = admit_method(method, dataset, status, verified, active)
        if admission["status"] != "admitted" or admission["native_validation"] != row["native_validation"]:
            raise ValueError("Fresh baseline native admission differs from retained result")
        truth = json.loads(Path(dataset["truth"]).read_text())
        sequences, owners, species = read_input_sequences(verified, truth, dataset["input"])
        pairs, artifacts = load_predictions(method, Path(dataset["methods"][method]["output"]), owners, species)
        if [record(path) for path in artifacts] != row["prediction_artifacts"]:
            raise ValueError("Baseline prediction artifacts differ from retained admission")
        parent = "orthofinder_full" if method == "orthofinder_sequence_only" else method
        inventory = {Path(r["absolute_path"]).resolve() for r in status["methods"][parent]["outputs"]}
        if not {path.resolve() for path in artifacts} <= inventory:
            raise ValueError("Baseline prediction artifact not in execution inventory")
        if score_pairs(pairs, truth["ortholog_pairs"], owners) != row["score"]:
            raise ValueError("Recomputed baseline counts or ratios differ")
        predictions[method] = pairs
        bindings.append({"method": method, "execution": item, "scheduler": row["scheduler"],
                         "admission": admission, "predictions": row["prediction_artifacts"],
                         "retained_score": row["score"]})
    return {"dataset": dataset, "verified": common_verified, "truth": truth,
            "sequences": sequences, "owners": owners, "predictions": predictions, "bindings": bindings}


def check_baselines(root):
    results = root / "benchmark_tools/results"
    data = {name: read_frozen(results / name, sha) for name, sha in PINS.items()}
    report = data["simulation_variable_native_results_20260916.json"]
    manifest = data["publication_variable_native_methods_20260916.json"]
    previous = data["publication_variable_methods_20260916.json"]
    if (report["panel"] != "variable_length_v2" or report["method_manifest_sha256"] != PINS["publication_variable_native_methods_20260916.json"]
            or report["generation_manifest_sha256"] != PINS["publication_variable_simulation_manifest_20260916.json"]
            or manifest["core_commit"] != CORE or manifest["core_commit"] != previous["core_commit"]
            or manifest["generation_manifest"] != previous["generation_manifest"]
            or manifest["generation_manifest"] != record(results / "publication_variable_simulation_manifest_20260916.json")
            or manifest["reused_comparator_manifest"] != record(results / "publication_variable_methods_20260916.json")
            or manifest["native_runtime"] != record(results / "publication_native_runtime_20260916.json")):
        raise ValueError("Retained panel and manifest identity bindings differ")
    rows = baseline_rows(report)
    executors = {}
    for kind, (directory, commit) in EXECUTORS.items():
        path = root / "benchmarks/work" / directory
        actual = subprocess.check_output(["git", "-C", str(path), "rev-parse", "HEAD"], text=True).strip()
        if actual != commit:
            raise ValueError("Retained executor revision differs")
        executors[kind] = path
    if (report["inference_executor_commit"] != EXECUTORS["native"][1]
            or report["reused_comparator_execution"]["executor_commit"] != EXECUTORS["comparator"][1]):
        raise ValueError("Retained assembly executor identity differs")
    for item in (manifest, previous):
        verify_scoring_dependencies(item)
    accounting = {"native": terminal_tasks(report["accounting_raw"], 21142, 70),
                  "comparator": terminal_tasks(report["reused_comparator_execution"]["accounting_raw"], 21010, 70)}
    old_datasets = {d["label"]: d for d in previous["datasets"]}
    if len(old_datasets) != len(previous["datasets"]):
        raise ValueError("Duplicate comparator dataset identity")
    controls = []
    for index, dataset in enumerate(manifest["datasets"]):
        if dataset["condition"] != "baseline":
            continue
        if dataset["label"] != f"baseline_{dataset['seed']}" or dataset["seed"] not in SEEDS:
            raise ValueError("Unknown baseline dataset")
        old = old_datasets[dataset["label"]]
        if any(dataset[k] != old[k] for k in ("seed", "condition", "input", "truth")):
            raise ValueError("Comparator dataset or truth differs")
        if any(dataset["methods"][method] != old["methods"][method] for method in METHODS[2:]):
            raise ValueError("Reused comparator configuration differs")
        tasks = {kind: values[index] for kind, values in accounting.items()}
        controls.append(bind_baseline(dataset, rows, manifest, previous, tasks, executors))
    if len(controls) != len(SEEDS) or {c["dataset"]["seed"] for c in controls} != set(SEEDS):
        raise ValueError("Incomplete baseline dataset inventory")
    return sorted(controls, key=lambda c: c["dataset"]["seed"]), manifest


def prepare(root, output, inference_root):
    if output.exists() or inference_root.exists():
        raise FileExistsError("Existing fragment inputs or inference artifacts; no restart")
    controls, manifest = check_baselines(root)
    prepared = []
    # Check every planned command and transform before writing any dataset.
    for control in controls:
        baseline = control["dataset"]
        label = f"{CONDITION}_{baseline['seed']}"
        inputs = output / label / "input"
        methods = fresh_methods(baseline, inputs, inference_root / label)
        sequences, coordinates = observe(control["sequences"], control["owners"], baseline["seed"])
        verification = verify_observation(control["sequences"], sequences, control["owners"], coordinates, baseline["seed"])
        prepared.append((control, label, inputs, methods, sequences, coordinates, verification))
    output.mkdir(parents=True, exist_ok=False)
    report = {"schema": "controlled_fragment_observations_v1", "status": "preparing",
              "condition": CONDITION, "protocol": record(root / "benchmark_tools/results" / PROTOCOL),
              "baseline_pins": [record(root / "benchmark_tools/results" / name) for name in PINS],
              "source": record(__file__), "transform_source": record(Path(__file__).with_name("controlled_fragment_observations.py")),
              "inference_order": list(INFERENCE_METHODS), "inference_launched": False,
              "new_accuracy_evaluated": False, "datasets": [],
              "runtime_scope": "Parent runtime bound for retained admission; current launch environment requires separate preflight",
              "limitations": ["Development-exposed synthetic observation test, not natural fragment truth or independent biology.",
                              "Shared-host resources are descriptive; cached control costs are not paired timings."]}
    try:
        for control, label, inputs, methods, sequences, coordinates, verification in prepared:
            inputs.mkdir(parents=True, exist_ok=False)
            for original in control["verified"]["inputs"]:
                species = Path(original["absolute_path"]).stem
                genes = [g for g in sorted(sequences) if control["owners"][g] == species]
                SeqIO.write([SeqRecord(Seq(sequences[g]), id=g, description="") for g in genes], inputs / (species + ".fasta"), "fasta")
            truth_path = inputs.parent / "parent_truth.json"
            shutil.copyfile(control["dataset"]["truth"], truth_path)
            verify_file(truth_path, control["verified"]["truth"])
            coordinate_path = inputs.parent / "coordinates.json"
            coordinate_path.write_text(json.dumps(coordinates, indent=2, sort_keys=True) + "\n")
            flags = {row["gene"]: row["fragment"] for row in coordinates}
            baseline = [{"arm": "baseline", "seed": control["dataset"]["seed"], "method": method,
                         "status": "complete", **score_strata(control["predictions"][method], control["truth"]["ortholog_pairs"], control["owners"], flags)}
                        for method in METHODS]
            report["datasets"].append({"label": label, "condition": CONDITION, "seed": control["dataset"]["seed"],
                                       "parent": control["dataset"]["label"], "input": str(inputs), "truth": str(truth_path),
                                       "verified_inputs": {"status": "ready", "truth": record(truth_path),
                                                           "inputs": [record(p) for p in sorted(inputs.iterdir())]},
                                       "parent_inputs": control["verified"], "baseline_admission": control["bindings"],
                                       "coordinates": record(coordinate_path), "verification": verification,
                                       "methods": methods, "baseline_records": baseline})
        report.update(status="prepared_unexecuted", inference_identities=len(prepared) * len(INFERENCE_METHODS))
    except Exception as error:
        report.update(status="preparation_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        (output / "manifest.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--inference-root", type=Path)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    root = args.root.resolve()
    if args.check_only:
        controls, _ = check_baselines(root)
        print(json.dumps({"status": "retained_baselines_verified", "datasets": len(controls),
                          "methods": len(controls) * len(METHODS), "new_inference": False,
                          "genes": {str(c['dataset']['seed']): len(c['owners']) for c in controls}}, sort_keys=True))
    else:
        if args.output is None or args.inference_root is None:
            parser.error("--output and --inference-root are required for preparation")
        report = prepare(root, args.output.resolve(), args.inference_root.resolve())
        print(json.dumps({"status": report["status"], "datasets": len(report["datasets"]),
                          "inference_identities": report["inference_identities"]}))


if __name__ == "__main__":
    main()
