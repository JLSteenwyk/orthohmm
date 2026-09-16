"""Freeze method commands/provenance before scientific simulation evaluation."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.benchmark_production import file_record


CORE_COMMIT = "7f3a9e40dd7e79f842cc2c11fb8b548f9a802806"
GENERATION_HASH = "ee31ea38d3b5c047abf80f04636f06838959f941d6704649a216bb93958343b2"


def commands(dataset, output, frozen, python, orthofinder):
    inputs = Path(dataset["input"])
    base = [str(python), str(frozen / "benchmark_tools/benchmark_production.py"), str(inputs)]
    settings = ["--cpu", "4", "--threads-per-worker", "4", "--matrix", "BLOSUM62",
                "--evalue", "0.0001", "--clustering", "leiden", "--cpm-resolution", "0.1",
                "--refinement-profile", "default", "--accuracy-profile", "high_sensitivity"]
    records = {}
    for method in ("orthohmm_high_sensitivity", "orthohmm_satellite_v2"):
        directory = output / method
        argv = [*base, str(directory), str(output / f"{method}.json"), *settings]
        if method == "orthohmm_satellite_v2":
            argv += ["--phylogeny", "reconcile", "--species-tree-mode", "infer",
                     "--aligner", "mafft", "--tree-builder", "FastTree",
                     "--phylogeny-candidates", "satellite_v2", "--phylogeny-root-rule", "species_overlap",
                     "--phylogeny-pair-rule", "positive_paralogy", "--species-tree-rooting", "min_variance"]
        records[method] = {"argv": argv, "output": str(directory), "metrics": str(output / f"{method}.json"),
                           "semantics": "group-derived cross-species pairs" if method.endswith("sensitivity") else "native phylogenetic ortholog pairs"}
    of_output = output / "orthofinder_full"
    copy = of_output / "input"
    records["orthofinder_full"] = {"argv": [str(orthofinder), "-f", str(copy), "-t", "4", "-a", "4", "-S", "diamond"],
        "output": str(of_output), "copy_inputs_from": str(inputs), "copy_inputs_to": str(copy),
        "semantics": "native phylogenetic ortholog pairs"}
    records["orthofinder_sequence_only"] = {"parent_method": "orthofinder_full", "output": str(of_output),
        "semantics": "MCL checkpoint group-derived pairs; diagnostic, no independent timing"}
    return records


def prepare(generation_path, frozen, python, orthofinder, output_root, destination, generation_hash=GENERATION_HASH):
    raw = generation_path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != generation_hash:
        raise ValueError("Wrong frozen simulation generation manifest")
    generation = json.loads(raw)
    if generation.get("status") != "materialized_not_executed" or len(generation.get("datasets", [])) != 70:
        raise ValueError("Expected complete frozen generation panel")
    if destination.exists() or output_root.exists():
        raise FileExistsError("Refusing existing method manifest or output root")
    frozen = frozen.resolve()
    commit = subprocess.check_output(["git", "-C", str(frozen), "rev-parse", "HEAD"], text=True).strip()
    if commit != CORE_COMMIT:
        raise ValueError("Unexpected OrthoHMM source revision")
    subprocess.run(["git", "-C", str(frozen), "diff", "--exit-code", "HEAD", "--", "orthohmm", "benchmark_tools/benchmark_production.py"], check=True, capture_output=True)
    tracked = subprocess.check_output(["git", "-C", str(frozen), "ls-files", "orthohmm"], text=True).splitlines()
    core_sources = [frozen / p for p in tracked] + [frozen / "benchmark_tools/benchmark_production.py"]
    python = python.absolute()
    orthofinder = orthofinder.absolute()
    inventory_code = "import importlib.metadata as m,json,sys; names=sorted({d.metadata['Name'] for d in m.distributions() if d.metadata['Name']}); print(json.dumps({'python':sys.version,'packages':{n:m.version(n) for n in names}}))"
    environments = {"orthohmm": json.loads(subprocess.check_output([str(python), "-c", inventory_code], text=True)),
                    "orthofinder": json.loads(subprocess.check_output([str(orthofinder.parent / "python"), "-c", inventory_code], text=True))}
    of_version = environments["orthofinder"]["packages"].get("orthofinder")
    if of_version != "3.1.5":
        raise ValueError(f"Expected OrthoFinder 3.1.5, found {of_version}")
    query = "import importlib.metadata as m,json; d=m.distribution('orthofinder'); print(json.dumps([str(d.locate_file(p).resolve()) for p in d.files if not str(p).endswith('.pyc')]))"
    of_files = json.loads(subprocess.check_output([str(orthofinder.parent / "python"), "-c", query], text=True))
    tools = {"orthofinder": orthofinder, "orthohmm_python": python}
    for name in ("mafft", "FastTree", "diamond"):
        path = shutil.which(name)
        if path is None:
            raise ValueError(f"Missing required tool: {name}")
        tools[name] = Path(path).absolute()
    workflow = Path(__file__).resolve().parent
    adapters = [workflow / name for name in ("simulation_method_outputs.py", "simulation_conditions.py",
        "orthofinder_to_pairwise.py", "orthofinder_mcl_to_orthogroups.py", "report_ygob_validation.py",
        "score_ygob_groups.py", "summarize_simulation_panel.py", "prepare_simulation_methods.py", "benchmark_production.py",
        "validate_simulation_outputs.py", "verify_ygob_validation.py", "run_simulation_generation.py")]
    report = {"schema_version": 1, "status": "frozen_not_executed", "core_commit": commit,
        "core_root": str(frozen), "generation_manifest": dict(file_record(generation_path, generation_path.parent), absolute_path=str(generation_path.resolve())),
        "core_sources": [dict(file_record(p, frozen), absolute_path=str(p)) for p in core_sources],
        "adapter_sources": [dict(file_record(p, workflow), absolute_path=str(p)) for p in adapters],
        "tool_entrypoints": {name: dict(file_record(p, p.parent), absolute_path=str(p)) for name, p in tools.items()},
        "orthofinder_distribution": [dict(file_record(Path(p), Path(p).parent), absolute_path=p) for p in sorted(of_files)],
        "environments": environments,
        "environment_overrides": {"PYTHONPATH": str(frozen), "PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1",
                                  "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"},
        "prepend_path": list(dict.fromkeys(str(tools[n].parent) for n in ("mafft", "FastTree", "diamond"))),
        "execution_order": ["orthohmm_high_sensitivity", "orthohmm_satellite_v2", "orthofinder_full"],
        "datasets": [], "limitations": ["Entrypoint hashes do not inventory all external-tool auxiliary files.",
            "Package inventories do not establish a portable environment rebuild.",
            "Generation completion, input hashes and native history equivalence require verification before inference.",
            "Actual resolved subprocess executables must be recorded; prepend_path alone does not prove runtime resolution.",
            "Truth and scientific scoring remain separate from inference inputs."]}
    for dataset in generation["datasets"]:
        label = f"{dataset['condition']}_{dataset['seed']}"
        report["datasets"].append({**dataset, "label": label,
                                  "methods": commands(dataset, output_root / label, frozen, python, orthofinder)})
    if len(report["datasets"]) != 70 or len({r["label"] for r in report["datasets"]}) != 70:
        raise ValueError("Not the complete seventy-dataset panel")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generation-manifest", type=Path, required=True)
    parser.add_argument("--generation-manifest-sha256", default=GENERATION_HASH)
    parser.add_argument("--frozen-root", type=Path, required=True)
    parser.add_argument("--python", type=Path, required=True)
    parser.add_argument("--orthofinder", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    report = prepare(args.generation_manifest.resolve(), args.frozen_root, args.python, args.orthofinder,
                     args.output_root.resolve(), args.manifest.resolve(), args.generation_manifest_sha256)
    print(f"Froze {len(report['datasets'])} method configurations; no inference")


if __name__ == "__main__":
    main()
