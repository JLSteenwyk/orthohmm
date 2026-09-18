"""Stage corrected FastOMA inputs only after corrected OrthoFinder admission."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_qfo_corrected_high_sensitivity import PLAN_SHA
from benchmark_tools.admit_wgd_comparators import check_root_tree
from benchmark_tools.fastoma_to_pairwise import input_owners
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_qfo_corrected_orthofinder_pairs import ADMITTER, validate_admission
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.validate_scaling_outputs import input_universe
from benchmark_tools.verify_ygob_validation import require_completed_job

ASSETS_SHA = "bdb6878dbab9396c6cb5360d35625ebe44149cf422b1d18dbf9b0e693a2b7975"
PROBE_SHA = "f2f2e8d62acfee09311536ceb458746de8824936f4a10f07640b1753c7c5e761"


def validate_tree_binding(content, admission, primary, inputs):
    tree = content["species_tree"]
    expected_root = Path(primary["methods"]["orthofinder_full"]["output"])
    results = Path(content["results_directory"])
    if not results.is_relative_to(expected_root):
        raise ValueError("Tree results outside corrected OrthoFinder output")
    if tree["path"] != str(results / "Species_Tree/SpeciesTree_rooted_node_labels.txt"):
        raise ValueError("Unexpected admitted species tree path")
    if tree not in admission["checked_records"] or any(r not in admission["checked_records"] for r in inputs):
        raise ValueError("Tree or inputs missing from admitted records")
    if len(inputs) != 78 or len({Path(r["path"]).stem for r in inputs}) != 78:
        raise ValueError("Require 78 distinct corrected species")
    if any(Path(r["path"]).parent != Path(primary["input_directory"])
           or Path(r["path"]).suffix != ".fasta" for r in inputs):
        raise ValueError("Inputs are not corrected primary FASTAs")
    return tree


def copy_inputs(inputs, tree, directory):
    """Create fresh copies and independently compare their content checksums."""
    directory.mkdir(parents=True, exist_ok=False)
    proteomes = directory / "proteome"
    proteomes.mkdir()
    copied = []
    targets = [proteomes / (Path(r["path"]).stem + ".fa") for r in inputs]
    if len(set(targets)) != len(targets):
        raise ValueError("Duplicate staged proteome filename")
    for source, target in zip([*inputs, tree], [*targets, directory / "species_tree.nwk"]):
        check(source)
        with Path(source["path"]).open("rb") as reader, target.open("xb") as writer:
            shutil.copyfileobj(reader, writer, length=1024 * 1024)
        observed = record(target)
        if any(observed[k] != source[k] for k in ("bytes", "sha256")):
            raise ValueError("Staged content differs from admitted source")
        check(source)
        copied.append({"source": source, "staged": observed})
    return copied


def prepare(root, admission_path, admission_sha, admission_job, directory, destination):
    if directory.exists() or destination.exists():
        raise FileExistsError("Require fresh stage and manifest paths")
    accounting = subprocess.check_output(["sacct", "-j", str(admission_job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = require_completed_job(accounting, admission_job)
    if scheduler["NodeList"] != "bizon" or scheduler["AllocCPUS"] != "2":
        raise ValueError("Wrong native admission allocation")
    admission = read_frozen(admission_path, admission_sha)
    content = validate_admission(admission)
    executor = root / "benchmarks/work/publication_qfo_corrected_orthofinder_admission_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != ADMITTER:
        raise ValueError("Native admission executor changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    if admission["source"] != record(executor / "benchmark_tools/admit_qfo_corrected_orthofinder.py"):
        raise ValueError("Wrong native admission source")
    results = root / "benchmark_tools/results"
    primary_path = results / "qfo_corrected_primary_commands_20260918.json"
    primary = read_frozen(primary_path, PLAN_SHA)
    assets_path = results / "qfo_corrected_fastoma_assets_20260918.json"
    assets = read_frozen(assets_path, ASSETS_SHA)
    probe_path = results / "fastoma_corrected_resource_probe_20260918.json"
    probe = read_frozen(probe_path, PROBE_SHA)
    if record(primary_path) not in admission["checked_records"]:
        raise ValueError("Admission not bound to corrected primary plan")
    inputs = assets["input_fastas"]
    if any(r not in primary["inputs"] for r in inputs):
        raise ValueError("FastOMA asset inputs differ from corrected primary plan")
    tree = validate_tree_binding(content, admission, primary, inputs)
    checked = [record(admission_path), admission["source"], record(primary_path), record(assets_path),
               record(probe_path), *inputs, tree, *probe["checked_records"]]
    for item in checked:
        check(item)
    _, species = input_universe({"inputs": inputs, "proteins": 984137, "proteomes": 78})
    if len(input_owners([Path(row["path"]) for row in inputs])) != 984137:
        raise ValueError("FastOMA accession transformation changes input universe")
    check_root_tree(Path(tree["path"]), species)
    copied = copy_inputs(inputs, tree, directory)
    check_root_tree(directory / "species_tree.nwk", species)
    for item in checked + [row["staged"] for row in copied]:
        check(item)
    report = {"status": "corrected_fastoma_inputs_staged_pending_launch_freeze",
              "source": record(__file__), "admission": record(admission_path),
              "admission_scheduler": scheduler, "admission_accounting": accounting,
              "checked_records": checked, "copies": copied, "input_directory": str(directory),
              "input_proteins": 984137, "input_species": 78,
              "execution_authorized": False, "accuracy_evaluated": False, "publication_ready": False,
              "limitations": ["Supplied corrected OrthoFinder tree, not independent FastOMA tree inference.",
                  "Only admitted tree and input assets rechecked; this does not repeat the full OrthoFinder output audit.",
                  "No inference launched; fresh command/runtime freeze and native-output admission remain required."]}
    with destination.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "admission", "directory", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--admission-sha256", required=True)
    parser.add_argument("--admission-job", type=int, required=True)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.admission.resolve(), args.admission_sha256,
            args.admission_job, args.directory.resolve(), args.output.resolve())
