"""Check frozen YGOB commands and native group conversions without scoring."""

import argparse
import json
from pathlib import Path
import shlex
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.report_ygob_validation import read_checkpoint
from benchmark_tools.score_ygob_groups import membership, read_predictions
from benchmark_tools.verify_ygob_reference import verify as verify_reference
from benchmark_tools.verify_ygob_validation import verify_run

PYTHON = "/home/bizon/anaconda3/bin/python"
OF = "/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/orthofinder_3.1.5/bin/orthofinder"
CHECKPOINT_PATTERN = "**/WorkingDirectory/clusters_OrthoFinder_I1.2.txt_id_pairs.txt"


def expected_command(root, method):
    output = root / "benchmarks/results/ygob_validation_v1"
    args = [str(root / "benchmarks/work/ygob_validation_v1/input"), "-o", str(output / method),
            "-c", "32", "--threads_per_worker", "8", "-x", "BLOSUM62", "-e", "0.0001",
            "--clustering", "leiden", "--cpm_resolution", "0.1", "--refinement_profile", "default",
            "--accuracy_profile", "high_sensitivity", "--metrics_json", str(output / (method + ".json")),
            "--stop", "infer"]
    if method == "satellite_v2":
        args += ["--phylogeny", "reconcile", "--species_tree_mode", "infer", "--aligner", "mafft",
                 "--tree_builder", "FastTree", "--phylogeny_candidates", "satellite_v2",
                 "--phylogeny_root_rule", "species_overlap", "--phylogeny_pair_rule", "positive_paralogy",
                 "--species_tree_rooting", "min_variance"]
    elif method != "high_sensitivity":
        raise ValueError("Unknown frozen method")
    return args


def verify_command(metrics, root, method):
    args = expected_command(root, method)
    native = [PYTHON, str(root / "benchmarks/work/publication_method_native_v2/orthohmm/__main__.py"), *args]
    harness = [PYTHON, "-m", "orthohmm", *args]
    if metrics["command"] != native or metrics["harness"]["command"] != harness:
        raise ValueError("Native or harness command differs from frozen specification")
    return {"native": native, "harness": harness}


def check_groups(groups, universe):
    found = set(membership(groups))
    if not found <= universe:
        raise ValueError("Native prediction contains foreign IDs")
    return {"groups": len(groups), "genes": len(found), "missing_input_genes": len(universe - found),
            "singletons": sum(len(g) == 1 for g in groups.values())}


def unique_path(root, pattern):
    paths = list(root.glob(pattern))
    if len(paths) != 1:
        raise ValueError("Missing or ambiguous native output: " + pattern)
    return paths[0]


def verify(root):
    root = root.resolve()
    output = root / "benchmarks/results/ygob_validation_v1"
    accounting = subprocess.check_output(["sacct", "-j", "21192", "--parsable2",
                                         "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    files = verify_run(root, root / "benchmarks/work/publication_method_native_v2", 21192, accounting)
    reference = verify_reference(root)
    from Bio import SeqIO
    universe = {r.id for item in reference["prepared_inputs"] for r in SeqIO.parse(item["path"], "fasta")}
    commands = {}
    for method in ("high_sensitivity", "satellite_v2"):
        commands[method] = verify_command(json.loads((output / (method + ".json")).read_text()), root, method)
    environment = json.loads((output / "orthofinder_environment.json").read_text())
    if [v for name, v in environment["packages"] if name.lower() == "orthofinder"] != ["3.1.5"]:
        raise ValueError("OrthoFinder recorded version differs")
    timing = (output / "orthofinder.time.log").read_text().splitlines()
    timed = [line.strip().removeprefix('Command being timed: "').removesuffix('"')
             for line in timing if line.strip().startswith("Command being timed:")]
    expected_of = [OF, "-f", str(output / "orthofinder/input"), "-t", "32", "-a", "8", "-S", "diamond"]
    if len(timed) != 1 or shlex.split(timed[0]) != expected_of:
        raise ValueError("OrthoFinder executed command differs")
    log = (output / "orthofinder.log").read_text()
    if "Starting OrthoFinder v3.1.5" not in log or "Running with the recommended MSA tree inference by default" not in log:
        raise ValueError("OrthoFinder log does not confirm full default MSA inference")
    commands["orthofinder"] = expected_of
    paths = {"orthohmm_high_sensitivity": (output / "high_sensitivity/orthohmm_orthogroups.txt", "named_groups"),
             "orthohmm_satellite_v2": (output / "satellite_v2/orthohmm_phylogeny/orthohmm_root_hogs.tsv", "root_hogs"),
             "orthofinder_full": (unique_path(output / "orthofinder", "**/Orthogroups/Orthogroups.txt"), "named_groups")}
    converted = {}
    for method, (path, format) in paths.items():
        converted[method] = {**check_groups(read_predictions(path, format), universe),
                             "native": file_provenance(path), "format": format}
    checkpoint = unique_path(output / "orthofinder", CHECKPOINT_PATTERN)
    ids = checkpoint.parent / "SequenceIDs.txt"
    groups = read_checkpoint(checkpoint, ids, universe)
    converted["orthofinder_sequence_only"] = {**check_groups(groups, universe), "native": file_provenance(checkpoint),
                                              "sequence_ids": file_provenance(ids), "format": "native_mcl_checkpoint"}
    return {"schema_version": 1, "status": "native_commands_and_conversion_verified", "accuracy_evaluated": False,
            "all_scoring_gates_verified": False, "file_verification": files, "reference_reconstruction": reference,
            "commands": commands, "native_groups": converted, "orthofinder_environment": file_provenance(output / "orthofinder_environment.json"),
            "verifier": file_provenance(Path(__file__)),
            "remaining": ["Overlap/reference-resource admission audit", "Independent score arithmetic and frozen uncertainty assembly"],
            "limitations": ["Native group membership, not resolved pairwise orthology.",
                            "Recorded package versions and entrypoint hashes do not inventory all external executable dependencies.",
                            "Shared-machine timings are not controlled efficiency comparisons."]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = verify(args.root)
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    main()
