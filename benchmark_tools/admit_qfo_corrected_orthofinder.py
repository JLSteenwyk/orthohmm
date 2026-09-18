"""Admit corrected OrthoFinder native evidence before conversion or scoring."""

import argparse
import json
from pathlib import Path
import shlex
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_qfo_corrected_high_sensitivity import PLAN_SHA, EXECUTOR
from benchmark_tools.admit_wgd_comparators import check_root_tree
from benchmark_tools.audit_orthofinder_input_parity import read_maps, compare_species
from benchmark_tools.audit_orthofinder_pair_tables import audit_tables
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.report_ygob_validation import read_checkpoint
from benchmark_tools.run_qfo_corrected_primary import verify
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.score_ygob_groups import membership
from benchmark_tools.validate_scaling_outputs import input_universe, unique_path
from benchmark_tools.validate_simulation_outputs import validate_graph_weights
from benchmark_tools.verify_ygob_validation import require_completed_job

METHOD = "orthofinder_full"


def validate_execution(plan, execution, scheduler, plan_record, runner_record):
    config = plan["methods"][METHOD]
    if (execution.get("status") != "process_succeeded_pending_native_admission"
            or execution.get("exit_code") != 0):
        raise ValueError("Native execution has not succeeded")
    if (scheduler["State"] != "COMPLETED" or scheduler["ExitCode"] != "0:0"
            or scheduler["NodeList"] != "bizon" or scheduler["AllocCPUS"] != "32"
            or execution["job_id"] != scheduler["JobIDRaw"] or execution["node"] != "bizon"):
        raise ValueError("Wrong scheduler identity/allocation")
    if (execution["method"] != METHOD or type(execution["index"]) is not int or execution["index"] != 1
            or execution["array_task_id"] != "1" or not execution["array_job_id"]):
        raise ValueError("Wrong primary method/task")
    if execution["source"] != runner_record or execution["manifest"] != plan_record:
        raise ValueError("Source/plan binding differs")
    if execution["native_argv"] != config["native_argv"] or execution["cwd"] != config["cwd"]:
        raise ValueError("Native command/cwd differs")
    if not execution["finished_epoch"] >= execution["started_epoch"] > 0:
        raise ValueError("Invalid execution timestamps")
    if execution["accuracy_admitted"] is not False or execution["native_outputs_validated"] is not False:
        raise ValueError("Unexpected native admission state")
    directory = Path(plan["output_root"]) / "execution" / METHOD
    for key, name in (("log", "native.log"), ("timing", "time.txt")):
        if execution[key]["path"] != str(directory / name) or execution[key]["bytes"] <= 0:
            raise ValueError("Invalid execution log/timing")


def validate_content(config, inputs, proteins=984137, proteomes=78):
    owners, species = input_universe({"inputs": inputs, "proteins": proteins, "proteomes": proteomes})
    output = Path(config["output"])
    log = unique_path(output, "**/Log.txt")
    lines = log.read_text().splitlines()
    starts = [i for i, s in enumerate(lines) if "Started OrthoFinder version " in s]
    ends = [i for i, s in enumerate(lines) if s.endswith(" : OrthoFinder run completed")]
    commands = [shlex.split(s.removeprefix("Command Line: ")) for s in lines if s.startswith("Command Line: ")]
    if (len(starts) != 1 or not lines[starts[0]].endswith("Started OrthoFinder version 3.1.5")
            or len(ends) != 1 or ends[0] <= starts[0] or commands != [config["native_argv"]]):
        raise ValueError("Native OrthoFinder version/command/completion differs")
    expected = {Path(r["path"]).name: r for r in inputs}
    copied = [record(p) for p in sorted(Path(config["copy_inputs_to"]).iterdir()) if p.is_file()]
    if ({Path(r["path"]).name: (r["bytes"], r["sha256"]) for r in copied}
            != {name: (r["bytes"], r["sha256"]) for name, r in expected.items()}):
        raise ValueError("Native input copies differ")
    working = log.parent / "WorkingDirectory"
    species_map, sequences = read_maps(working / "SpeciesIDs.txt", working / "SequenceIDs.txt")
    if set(species_map.values()) != set(expected):
        raise ValueError("Internal species inventory differs")
    if {p.name for p in working.glob("Species*.fa")} != {f"Species{i}.fa" for i in species_map}:
        raise ValueError("Internal FASTA inventory differs")
    sequence_rows = []
    for index, name in species_map.items():
        row = compare_species(Path(expected[name]["path"]), working / f"Species{index}.fa", sequences[index])
        if row["differences"]:
            raise ValueError("OrthoFinder internal sequences differ from corrected inputs")
        sequence_rows.append({"species_index": index, "filename": name, **row})
    graph = working / "OrthoFinder_graph.txt"
    validate_graph_weights(graph)
    clusters = unique_path(working, "clusters_OrthoFinder_I*.txt_id_pairs.txt")
    groups = read_checkpoint(clusters, working / "SequenceIDs.txt", owners)
    tree = log.parent / "Species_Tree/SpeciesTree_rooted_node_labels.txt"
    check_root_tree(tree, species)
    pairs = audit_tables(log.parent, owners, species, membership(groups))
    return {"results_directory": str(log.parent), "genes": len(owners), "species": len(species),
            "internal_sequences": sequence_rows, "checkpoint_groups": len(groups),
            "checkpoint": record(clusters), "sequence_ids": record(working / "SequenceIDs.txt"),
            "species_tree": record(tree), "native_pairs": pairs, "copied_inputs": copied}


def admit(root, job, destination):
    if destination.exists():
        raise FileExistsError(destination)
    accounting = subprocess.check_output(["sacct", "-j", str(job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = require_completed_job(accounting, job)
    plan_path = root / "benchmark_tools/results/qfo_corrected_primary_commands_20260918.json"
    plan = read_frozen(plan_path, PLAN_SHA)
    config = plan["methods"][METHOD]
    executor = root / "benchmarks/work/publication_qfo_corrected_primary_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != EXECUTOR:
        raise ValueError("Inference executor revision changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    status_path = Path(plan["output_root"]) / "execution" / METHOD / "status.json"
    status_record = record(status_path)
    execution = json.loads(status_path.read_text())
    validate_execution(plan, execution, scheduler, record(plan_path),
                       record(executor / "benchmark_tools/run_qfo_corrected_primary.py"))
    _, inputs = verify(plan)
    output = Path(config["output"])
    observed = [record(p) for p in sorted(output.rglob("*")) if p.is_file()]
    if not observed or observed != execution["outputs"]:
        raise ValueError("Native output inventory/content differs")
    content = validate_content(config, inputs)
    if content["copied_inputs"] != execution["copied_inputs"]:
        raise ValueError("Copied inputs differ from execution inventory")
    checked = [status_record, record(plan_path), execution["source"], execution["log"], execution["timing"],
               *plan["inputs"], *inputs, *observed]
    for item in checked:
        check(item)
    if {str(p) for p in output.rglob("*") if p.is_file()} != {r["path"] for r in observed}:
        raise ValueError("Output inventory changed during admission")
    report = {"status": "corrected_orthofinder_native_evidence_admitted", "source": record(__file__),
              "scheduler": scheduler, "accounting": accounting, "checked_records": checked,
              "content": content, "accuracy_evaluated": False, "publication_ready": False,
              "helpers": [record(Path(__file__).with_name(n)) for n in (
                  "run_qfo_corrected_primary.py", "audit_orthofinder_input_parity.py",
                  "audit_orthofinder_pair_tables.py", "admit_wgd_comparators.py",
                  "validate_scaling_outputs.py", "report_ygob_validation.py", "score_ygob_groups.py")],
              "limitations": ["Native integrity, not accuracy or a validated biological species tree.",
                  "Full native pairs and pre-phylogenetic MCL cliques must be converted and scored separately.",
                  "Shared-host resources are not matched dedicated timing."]}
    with destination.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--job", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    admit(args.root.resolve(), args.job, args.output.resolve())
