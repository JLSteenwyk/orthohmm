"""Run a bounded mechanical fixture; never authorize publication timing runs."""

import argparse
import importlib
import json
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.native_factorial_adapter import entrypoint, factors
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def diagnose(core, baseline_path, baseline_sha, fasta, output, cell):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    if record(baseline_path)["sha256"] != baseline_sha:
        raise ValueError("Baseline checksum differs")
    baseline = json.loads(baseline_path.read_text())
    if baseline["core_commit"] != "7f3a9e40dd7e79f842cc2c11fb8b548f9a802806" or baseline["core_root"] != str(core):
        raise ValueError("Wrong frozen core root or commit")
    expected_env = {"PYTHONHASHSEED": "0", "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    if any(os.environ.get(k) != v for k, v in expected_env.items()) or not sys.dont_write_bytecode:
        raise ValueError("Require fixed diagnostic environment and disabled bytecode writes")
    if not sys.pycache_prefix or Path(sys.pycache_prefix).exists() or Path(sys.pycache_prefix).is_symlink():
        raise ValueError("Require an absent, separate diagnostic bytecode cache prefix")
    files = sorted(fasta.glob("*.fa"))
    if len(files) != 4 or set(fasta.iterdir()) != set(files):
        raise ValueError("Require exactly four diagnostic FASTAs and no other files")
    inputs, genes = [], set()
    for path in files:
        inputs.append(record(path))
        for line in path.read_text().splitlines():
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in genes:
                    raise ValueError("Duplicate diagnostic gene")
                genes.add(name)
    if not 16 <= len(genes) <= 64:
        raise ValueError("Diagnostic only: require 16-64 proteins")
    sources = [{"path": r["absolute_path"], "bytes": r["bytes"], "sha256": r["sha256"]}
               for r in baseline["core_sources"]]
    for pin in sources:
        check(pin)
    tools = {name: {"path": baseline["tool_entrypoints"][name]["absolute_path"],
                   "bytes": baseline["tool_entrypoints"][name]["bytes"],
                   "sha256": baseline["tool_entrypoints"][name]["sha256"]}
             for name in ("mafft", "FastTree")}
    for pin in tools.values():
        check(pin)
    sys.path.insert(0, str(core))
    module = importlib.import_module("orthohmm.orthohmm")
    if Path(module.__file__).resolve() != core / "orthohmm/orthohmm.py":
        raise ValueError("Wrong frozen pipeline import")
    adapted, flags = entrypoint(module, cell)
    config = factors(cell)
    output.mkdir(parents=True, exist_ok=False)
    metrics_path = output / "metrics.json"
    result = dict(status="diagnostic_running", cell=cell, diagnostic_only=True,
        scientific_execution_authorized=False, publication_ready=False,
        baseline=record(baseline_path), source=record(__file__),
        adapter_source=record(Path(__file__).with_name("native_factorial_adapter.py")),
        core_sources=sources, inputs=inputs, tools=tools, factors=flags,
        python=sys.version, cpu=2, process_cpu_affinity=sorted(os.sched_getaffinity(0)),
        environment=expected_env, accuracy_evaluated=False)
    report = output / "diagnostic.json"
    report.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    try:
        try:
            adapted(fasta_directory=str(fasta), output_directory=str(output / "native"),
                phmmer="unused", cpu=2, threads_per_worker=2, single_copy_threshold=.5,
                mcl="unused", inflation_value=1.5, start=None, stop=module.StopStep.infer,
                substitution_matrix=module.SubstitutionMatrix.blosum62, evalue_threshold=.0001,
                search_mode="builtin", clustering="leiden", cpm_resolution=.1,
                refinement_profile="default", accuracy_profile="high_sensitivity", metrics_json=str(metrics_path),
                phylogeny="reconcile" if config["reconciliation"] else "off",
                phylogeny_candidates="satellite_v2" if config["candidate_expansion"] else "seed",
                species_tree_mode="infer", species_tree_rooting="min_variance",
                phylogeny_root_rule="species_overlap", phylogeny_pair_rule="positive_paralogy",
                aligner=tools["mafft"]["path"], tree_builder=tools["FastTree"]["path"])
        except SystemExit as error:
            if error.code not in (None, 0):
                raise
        metrics = json.loads(metrics_path.read_text())
        if metrics["status"] != "complete" or metrics["metadata"]["native_factorial"] != flags:
            raise ValueError("Diagnostic metrics incomplete or factors differ")
        expected_stages = {"search", "edge_thresholds", "network_edges", "clustering", "refinement", "orthogroup_materialization"}
        if config["profile_expansion"]:
            expected_stages.add("profile_expansion")
        if config["candidate_expansion"]:
            expected_stages.add("phylogeny_candidates")
        if config["reconciliation"]:
            expected_stages.add("phylogeny")
        if set(metrics["stages"]) != expected_stages or metrics["counts"]["genes"] != len(genes):
            raise ValueError("Diagnostic stage/gene inventory differs")
        if config["reconciliation"] and (metrics["counts"]["phylogeny_checkpoint_hits"] != 0
                or metrics["counts"]["phylogeny_species_tree_checkpoint_hit"] is not False):
            raise ValueError("Diagnostic reused reconciliation checkpoint")
        for pin in sources + inputs + list(tools.values()) + [result["adapter_source"], result["source"], result["baseline"]]:
            check(pin)
        result.update(status="native_factorial_mechanical_diagnostic_complete",
            genes=len(genes), metrics=record(metrics_path), counts=metrics["counts"], stages=sorted(metrics["stages"]),
            checkpoint=record(output / "native/orthohmm_working_res/high_sensitivity_checkpoint/manifest.json"),
            native_groups=record(output / "native/orthohmm_orthogroups.txt"),
            limitations=["Four-species fixture, not an evolutionary accuracy simulation or publication timing observation.",
                "Stage presence/checkpoints/outputs do not establish full-dataset equivalence or sensitivity.",
                "The production source and globals remain unchanged; factor metadata makes the in-memory adapter explicit.",
                "Full-dataset runs require their own prospective resource protocol, actual safe-capacity check and valid whole-run accounting."])
    except BaseException as error:
        result.update(status="diagnostic_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagnostic", action="store_true", required=True)
    for name in ("core-root", "baseline", "fasta-directory", "output-directory"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--baseline-sha256", required=True)
    parser.add_argument("--cell", required=True)
    args = parser.parse_args()
    diagnose(args.core_root.resolve(), args.baseline.resolve(), args.baseline_sha256,
             args.fasta_directory.resolve(), args.output_directory.absolute(), args.cell)


if __name__ == "__main__":
    main()
