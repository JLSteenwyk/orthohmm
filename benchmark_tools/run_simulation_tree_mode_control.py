"""Rerun a frozen simulation with each method's own species tree supplied."""

import argparse
from copy import deepcopy
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import Phylo
from benchmark_tools.assemble_simulation_results import admit_method, input_universe
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_species_tree_robustness import topology
from benchmark_tools.run_simulation_methods import (
    read_frozen, verify_environment, verify_inputs, execution_environment, execute,
)
from benchmark_tools.simulation_method_outputs import load_predictions, unique_path
from benchmark_tools.simulation_supplied_commands import fresh_supplied_method
from benchmark_tools.run_simulation_generation import verify_file

METHOD_SHA = "bf677728cf9312edd9755d9ac1ceefbce3abb8c3848d43413631dc00ca00eb2f"
RESULT_SHA = "cc99fc31c3433809098212d3c8dc12f4ad829b4cfeb66dabcb13aede4e95258f"
METHODS = ("orthohmm_satellite_v2", "orthofinder_full")


def native_tree(method, output):
    if method == "orthohmm_satellite_v2":
        return output / "orthohmm_phylogeny/species_tree.rooted.nwk"
    if method == "orthofinder_full":
        return unique_path(output, "**/Species_Tree/SpeciesTree_rooted.txt")
    raise ValueError("Unsupported method")


def artifact_inventory(method, output):
    if method == "orthohmm_satellite_v2":
        base = output
        paths = [base / "orthohmm_working_res/phylogeny_candidate_superfamilies.txt"]
        paths += sorted(base.glob("orthohmm_phylogeny/gene_trees/*.raw.nwk"))
        paths += sorted(base.glob("orthohmm_phylogeny/alignments/*.faa"))
    elif method == "orthofinder_full":
        base = unique_path(output, "**/WorkingDirectory")
        paths = [base / name for name in ("OrthoFinder_graph.txt", "SequenceIDs.txt", "SpeciesIDs.txt")]
        paths += sorted(base.glob("clusters_OrthoFinder*"))
        paths += sorted(base.glob("Trees_ids/*"))
        paths += sorted(base.glob("Alignments_ids/*"))
    else:
        raise ValueError("Unsupported method")
    if not paths or any(not p.is_file() for p in paths):
        raise ValueError("Incomplete retained artifact inventory")
    return {str(p.relative_to(base)): record(p) for p in paths}


def compare_inventory(before, after):
    missing, extra = sorted(set(before) - set(after)), sorted(set(after) - set(before))
    changed = sorted(k for k in set(before) & set(after)
                     if (before[k]["bytes"], before[k]["sha256"]) != (after[k]["bytes"], after[k]["sha256"]))
    return {"identical": not (missing or extra or changed), "missing": missing, "extra": extra,
            "changed": changed, "before_files": len(before), "after_files": len(after)}


def compare_pairs(before, after):
    a = {tuple(sorted(pair)) for pair in before}
    b = {tuple(sorted(pair)) for pair in after}
    return {"identical": a == b, "before_pairs": len(a), "after_pairs": len(b),
            "removed": [list(p) for p in sorted(a - b)], "added": [list(p) for p in sorted(b - a)]}


def baseline_records(dataset, evidence, manifest, verified):
    baseline = {}
    for method in METHODS:
        rows = [r for r in evidence["records"] if r["condition"] == dataset["condition"]
                and r["seed"] == dataset["seed"] and r["method"] == method]
        if len(rows) != 1 or rows[0]["status"] != "complete":
            raise ValueError("Mode control needs a completed inferred baseline: " + method)
        source = rows[0]["execution_evidence"]
        verify_file(Path(source["absolute_path"]), source)
        status = json.loads(Path(source["absolute_path"]).read_text())
        admission = admit_method(method, dataset, status, verified, manifest)
        if admission["status"] != "admitted":
            raise ValueError("Original output admission failed: " + method)
        output = Path(dataset["methods"][method]["output"])
        tree = native_tree(method, output)
        inventoried = {r["absolute_path"] for r in status["methods"][method]["outputs"]}
        if str(tree) not in inventoried:
            raise ValueError("Original species tree was not inventoried")
        baseline[method] = {"execution_evidence": source, "admission": admission,
                            "tree": record(tree), "retained_artifacts": artifact_inventory(method, output)}
    return baseline


def parser_checks(manifest, baseline, species, env):
    python = manifest["tool_entrypoints"]["orthohmm_python"]["absolute_path"]
    of_python = str(Path(manifest["tool_entrypoints"]["orthofinder"]["absolute_path"]).parent / "python")
    commands = [
        [python, "-c", "import json,sys; from orthohmm.phylogeny import parse_species_tree; "
         "parse_species_tree(sys.argv[1],json.loads(sys.argv[2])); print('accepted')",
         baseline["orthohmm_satellite_v2"]["tree"]["path"], json.dumps([s + ".fasta" for s in species])],
        [of_python, "-c", "import json,sys; from orthofinder.run.species_info import CheckUserSpeciesTree; "
         "CheckUserSpeciesTree(sys.argv[1],json.loads(sys.argv[2])); print('accepted')",
         baseline["orthofinder_full"]["tree"]["path"], json.dumps(species)],
    ]
    evidence = []
    for argv in commands:
        result = subprocess.run(argv, cwd=manifest["core_root"], env=env, text=True, capture_output=True, check=True)
        if result.stdout.strip() != "accepted":
            raise ValueError("Unexpected native parser output")
        evidence.append({"argv": argv, "cwd": manifest["core_root"], "stdout": result.stdout,
                         "stderr": result.stderr, "exit_code": result.returncode})
    return evidence


def run(root, label, output):
    if Path.cwd().resolve() != root:
        raise ValueError("Verify environment from original repository directory")
    if output.exists():
        raise FileExistsError(output)
    results = root / "benchmark_tools/results"
    manifest_path = results / "publication_variable_native_methods_20260916.json"
    evidence_path = results / "simulation_variable_native_results_20260916.json"
    manifest = read_frozen(manifest_path, METHOD_SHA)
    evidence = read_frozen(evidence_path, RESULT_SHA)
    matches = [d for d in manifest["datasets"] if d["label"] == label]
    if len(matches) != 1:
        raise ValueError("Unknown or duplicate simulation dataset")
    dataset = matches[0]
    generation_record = manifest["generation_manifest"]
    generation = read_frozen(Path(generation_record["absolute_path"]), generation_record["sha256"])
    verify_environment(manifest)
    verified = verify_inputs(dataset, generation, root / "benchmarks/work/publication_variable_simulation_panel_v2",
                             generation_record["sha256"])
    if verified["status"] != "ready":
        raise ValueError("Dataset is not applicable")
    baseline = baseline_records(dataset, evidence, manifest, verified)
    owners, species = input_universe(verified["inputs"], json.loads(Path(dataset["truth"]).read_text()))
    configured = deepcopy(dataset)
    configured["methods"] = {m: fresh_supplied_method(m, dataset["methods"][m], baseline[m]["tree"]["path"],
                                                     output / m) for m in METHODS}
    env, resolved = execution_environment(manifest)
    parser_evidence = parser_checks(manifest, baseline, species, env)
    sources = [record(Path(__file__).with_name(name)) for name in
               ("run_simulation_tree_mode_control.py", "simulation_supplied_commands.py", "run_simulation_methods.py",
                "assemble_simulation_results.py", "validate_simulation_outputs.py", "simulation_method_outputs.py")]
    provenance = {"source": record(__file__), "sources": sources, "manifest": record(manifest_path),
                  "baseline_results": record(evidence_path), "baseline": baseline, "configured": configured,
                  "native_parser_checks": parser_evidence,
                  "resolved_tools": resolved, "job_id": os.environ.get("SLURM_JOB_ID"),
                  "executor_commit": subprocess.check_output(["git", "-C", str(Path(__file__).resolve().parent.parent),
                                                               "rev-parse", "HEAD"], text=True).strip(),
                  "scope": "Fresh supplied-own-tree mode control; no accuracy scoring or matched runtime claim"}
    output.mkdir(parents=True)
    (output / "preflight.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    report = {"status": "running", "provenance": provenance, "accuracy_evaluated": False, "methods": {}}
    try:
        execution = execute(configured, list(METHODS), env, output / "execution", verified, provenance, manifest)
        report["execution"] = record(output / "execution/status.json")
        for method in METHODS:
            admission = admit_method(method, configured, execution, verified, manifest)
            row = {"admission": admission}
            report["methods"][method] = row
            if admission["status"] != "admitted":
                continue
            before = Path(dataset["methods"][method]["output"])
            after = Path(configured["methods"][method]["output"])
            old_pairs, _ = load_predictions(method, before, owners, species)
            new_pairs, _ = load_predictions(method, after, owners, species)
            row["pairs"] = compare_pairs(old_pairs, new_pairs)
            original_tree = Phylo.read(baseline[method]["tree"]["path"], "newick")
            resulting_tree = Phylo.read(native_tree(method, after), "newick")
            row["rooted_tree_identical"] = (topology(original_tree) == topology(resulting_tree)
                and sorted(t.name for t in original_tree.get_terminals()) == sorted(t.name for t in resulting_tree.get_terminals()))
            row["retained_artifacts"] = artifact_inventory(method, after)
            row["artifact_comparison"] = compare_inventory(baseline[method]["retained_artifacts"], row["retained_artifacts"])
        verify_environment(manifest)
        if baseline_records(dataset, evidence, manifest, verified) != baseline:
            raise ValueError("Original baseline changed during execution")
        for source in sources:
            check(source)
        read_frozen(manifest_path, METHOD_SHA)
        read_frozen(evidence_path, RESULT_SHA)
        for item in verified["inputs"]:
            verify_file(Path(item["absolute_path"]), item)
        report["status"] = "complete_pending_independent_admission"
    except Exception as error:
        report.update(status="failed", error=str(error))
        raise
    finally:
        (output / "results.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.root.resolve(), args.dataset, args.output.resolve())
