"""Validate a fresh DGX full-pipeline fixture and compare admitted x86 pairs."""

import argparse
import json
from pathlib import Path
import sys

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_dgx_orthohmm_smoke import compare
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.report_ygob_validation import read_checkpoint
from benchmark_tools.simulation_method_outputs import unique_path, orthofinder_pairs
from benchmark_tools.validate_scaling_outputs import validate_orthofinder

SIMULATION_SHA = "cc99fc31c3433809098212d3c8dc12f4ad829b4cfeb66dabcb13aede4e95258f"


def audit(arm, inputs, simulation):
    source = record(simulation)
    if source["sha256"] != SIMULATION_SHA:
        raise ValueError("Changed x86 admission")
    data = json.loads(simulation.read_text())
    selected = [r for r in data["records"] if r["condition"] == "missing20" and r["seed"] == 20261101
                and r["method"] == "orthofinder_full"]
    if len(selected) != 1 or selected[0]["status"] != "complete":
        raise ValueError("Missing admitted reference")
    artifacts = selected[0]["prediction_artifacts"]
    for artifact in artifacts:
        check({"path": artifact["absolute_path"], "bytes": artifact["bytes"], "sha256": artifact["sha256"]})
    reference = Path(artifacts[0]["absolute_path"]).parents[2]
    paths = sorted(inputs.glob("*.fasta"))
    owners = {}
    for path in paths:
        for seq in SeqIO.parse(path, "fasta"):
            if seq.id in owners or not seq.seq:
                raise ValueError("Invalid input universe")
            owners[seq.id] = path.stem
    species = [p.stem for p in paths]
    if len(species) != 8 or len(owners) != 645:
        raise ValueError("Wrong smoke fixture")
    transferred = {}
    for line in (arm / "input.before.sha256").read_text().splitlines():
        digest, path = line.split(maxsplit=1)
        if Path(path).name in transferred:
            raise ValueError("Duplicate input identity")
        transferred[Path(path).name] = digest
    if transferred != {p.name: record(p)["sha256"] for p in paths}:
        raise ValueError("Transferred fixture differs")
    before, after = arm / "runtime.before.json", arm / "runtime.after.json"
    if before.read_bytes() != after.read_bytes() or (arm / "exit.txt").read_text() != "exit=0\n":
        raise ValueError("Runtime changed or native command failed")
    remote = "/home/jlsteenwyk/projects/orthohmm-publication"
    output = remote + "/orthofinder_full_smoke_v1"
    run = {"configuration": {"output": str(arm / "results"), "copy_inputs_to": str(arm / "input")},
           "dataset": {"inputs": [record(p) for p in paths]},
           "native_argv": [remote + "/envs/orthofinder/bin/orthofinder", "-f", output + "/input",
                           "-t", "4", "-a", "4", "-S", "diamond", "-o", output + "/results"]}
    validation = validate_orthofinder(run, owners, species)
    comparisons = {}
    groups, pairs = [], []
    for root in (reference, arm / "results"):
        partition = read_checkpoint(unique_path(root, "**/clusters_OrthoFinder_I*.txt_id_pairs.txt"),
                                    unique_path(root, "**/SequenceIDs.txt"), owners)
        groups.append({tuple(sorted(genes)) for genes in partition.values()})
        rows = orthofinder_pairs(unique_path(root, "**/Orthologues").parent, owners, species)
        if len(rows) != len(set(rows)):
            raise ValueError("Converter emitted duplicate canonical pairs")
        pairs.append(set(rows))
    comparisons["mcl_checkpoint"] = compare(*groups)
    comparisons["native_pairs"] = compare(*pairs)
    check(source)
    for artifact in artifacts:
        check({"path": artifact["absolute_path"], "bytes": artifact["bytes"], "sha256": artifact["sha256"]})
    return {"status": "full_pipeline_fixture_audited", "publication_ready": False, "source": record(Path(__file__)),
            "reference_admission": source, "inputs": [record(p) for p in paths], "native_validation": validation,
            "runtime": record(before), "trace": record(arm / "execve.log"), "comparisons": comparisons,
            "limitations": ["One development-exposed fixture, not general equivalence or independent accuracy.",
                            "strace overhead present; not a scientific timing result.",
                            "Intermediate/final tree and HOG equality are not established by this comparison.",
                            "This fixture did not execute FastME; separate direct command probes cover finite matrices only."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("arm", "inputs", "simulation", "output"):
        parser.add_argument("--" + flag, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = audit(args.arm.resolve(), args.inputs.resolve(), args.simulation.resolve())
    with args.output.open("x") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
