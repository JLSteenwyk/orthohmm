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
from benchmark_tools.build_publication_runtime import verify_runtime
from benchmark_tools.validate_profile_runtime import require_profile_runtime


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


def reuse_comparators(report, previous):
    """Reuse configurations only; native completion/admission remain scoring gates."""
    if previous["generation_manifest"]["sha256"] != report["generation_manifest"]["sha256"]:
        raise ValueError("Comparator reuse requires identical generated datasets")
    if previous["core_commit"] != report["core_commit"]:
        raise ValueError("Runtime correction must retain the frozen scientific revision")
    for field in ("tool_entrypoints", "orthofinder_distribution", "environments"):
        if previous[field] != report[field]:
            raise ValueError("Comparator environment changed: " + field)
    old = {d["label"]: d for d in previous["datasets"]}
    if len(old) != len(previous["datasets"]) or set(old) != {d["label"] for d in report["datasets"]}:
        raise ValueError("Comparator dataset set differs")
    for dataset in report["datasets"]:
        source = old[dataset["label"]]
        if {k: v for k, v in dataset.items() if k != "methods"} != {k: v for k, v in source.items() if k != "methods"}:
            raise ValueError("Comparator dataset metadata differs")
        for name in ("orthohmm_high_sensitivity", "orthohmm_satellite_v2"):
            before, after = source["methods"][name]["argv"], dataset["methods"][name]["argv"]
            if before[0] != after[0] or before[2] != after[2] or before[5:] != after[5:]:
                raise ValueError("Runtime correction changed scientific arguments")
        old_output = Path(source["methods"]["orthofinder_full"]["output"]).parent
        expected = commands(dataset, old_output, Path(report["core_root"]),
                            Path(report["tool_entrypoints"]["orthohmm_python"]["absolute_path"]),
                            Path(report["tool_entrypoints"]["orthofinder"]["absolute_path"]))
        for name in ("orthofinder_full", "orthofinder_sequence_only"):
            if source["methods"][name] != expected[name]:
                raise ValueError("Comparator configuration differs from frozen settings")
            dataset["methods"][name] = source["methods"][name]
    report["execution_order"] = ["orthohmm_high_sensitivity", "orthohmm_satellite_v2"]


def prepare(generation_path, frozen, python, orthofinder, output_root, destination, generation_hash=GENERATION_HASH,
            runtime_manifest=None, reuse_manifest=None, reuse_hash=None):
    raw = generation_path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != generation_hash:
        raise ValueError("Wrong frozen simulation generation manifest")
    generation = json.loads(raw)
    if generation.get("status") != "materialized_not_executed" or len(generation.get("datasets", [])) != 70:
        raise ValueError("Expected complete frozen generation panel")
    if destination.exists() or output_root.exists():
        raise FileExistsError("Refusing existing method manifest or output root")
    if runtime_manifest is None:
        raise ValueError("A validated native runtime manifest is required")
    if (reuse_manifest is None) != (reuse_hash is None):
        raise ValueError("Comparator manifest and expected hash must be supplied together")
    frozen = frozen.resolve()
    verify_runtime(runtime_manifest, frozen)
    profile_probe = require_profile_runtime(frozen, str(python))
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
        "validate_simulation_outputs.py", "verify_ygob_validation.py", "run_simulation_generation.py",
        "build_publication_runtime.py", "validate_profile_runtime.py")]
    report = {"schema_version": 1, "status": "frozen_not_executed", "core_commit": commit,
        "core_root": str(frozen), "native_runtime": dict(file_record(runtime_manifest, runtime_manifest.parent), absolute_path=str(runtime_manifest.resolve())),
        "profile_probe": profile_probe,
        "generation_manifest": dict(file_record(generation_path, generation_path.parent), absolute_path=str(generation_path.resolve())),
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
    if reuse_manifest is not None:
        raw = reuse_manifest.read_bytes()
        if hashlib.sha256(raw).hexdigest() != reuse_hash:
            raise ValueError("Comparator reuse manifest changed")
        reuse_comparators(report, json.loads(raw))
        report["reused_comparator_manifest"] = dict(file_record(reuse_manifest, reuse_manifest.parent),
                                                   absolute_path=str(reuse_manifest.resolve()))
        report["limitations"].append("Reused comparator outputs require their original scheduler, executor and native admission evidence; reuse is not automatic admission.")
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
    parser.add_argument("--runtime-manifest", type=Path, required=True)
    parser.add_argument("--reuse-comparator-manifest", type=Path)
    parser.add_argument("--reuse-comparator-sha256")
    args = parser.parse_args()
    report = prepare(args.generation_manifest.resolve(), args.frozen_root, args.python, args.orthofinder,
                     args.output_root.resolve(), args.manifest.resolve(), args.generation_manifest_sha256,
                     args.runtime_manifest.resolve(), args.reuse_comparator_manifest, args.reuse_comparator_sha256)
    print(f"Froze {len(report['datasets'])} method configurations; no inference")


if __name__ == "__main__":
    main()
