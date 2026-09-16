"""Verify retained trace rows against native partitions and correct branch labels."""

import argparse
from collections import Counter
import csv
import itertools
import json
import math
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
from benchmark_tools.assemble_orthobench_factorial import load_reference_snapshot
from benchmark_tools.audit_historical_profile_ablation import verify_file
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.score_orthobench_partition import score_partition
from benchmark_tools.trace_ob_families import STAGES, TRANSITIONS, family_group_summary, partition, render, transitions
from benchmark_tools.verify_ygob_validation import require_completed_job

ORIGINAL_SHA = "82883d3d3b32a1c053775513d56987aefe943739fd532bcad4d16dea82a4f34f"


def verify_pairs(records, references, indices, owners):
    expected = {name: set(itertools.combinations(sorted(genes), 2)) for name, genes in references.items()}
    observed = {name: [] for name in references}
    fields = {"refog", "left", "right", "same_species", "forward_normalized_hit", "reverse_normalized_hit", *STAGES}
    for raw in records:
        if set(raw) != fields or raw["refog"] not in expected:
            raise ValueError("Changed pair-table columns or reference family")
        name, pair = raw["refog"], (raw["left"], raw["right"])
        if pair not in expected[name]:
            raise ValueError("Duplicate, reversed or unexpected reference pair")
        expected[name].remove(pair)
        row = {"left": pair[0], "right": pair[1]}
        for field in ("same_species", *STAGES):
            if raw[field] not in {"True", "False"}:
                raise ValueError("Invalid pair membership boolean")
            row[field] = raw[field] == "True"
            value = owners[pair[0]] == owners[pair[1]] if field == "same_species" else indices[field][pair[0]] == indices[field][pair[1]]
            if row[field] != value:
                raise ValueError("Pair table differs from native membership or species")
        for field in ("forward_normalized_hit", "reverse_normalized_hit"):
            value = None if raw[field] == "NA" else float(raw[field])
            if value is not None and (not math.isfinite(value) or value <= 0):
                raise ValueError("Invalid cached-hit value")
            row[field] = value
        observed[name].append(row)
    if any(expected.values()):
        raise ValueError("Missing reference pairs")
    return observed


def admit(root, output):
    if output.exists():
        raise FileExistsError(output)
    path = root / "benchmarks/results/ob_family_trace_v2/results.json"
    report = read_frozen(path, ORIGINAL_SHA)
    accounting = subprocess.check_output(["sacct", "-j", "21313", "--parsable2", "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = require_completed_job(accounting, 21313)
    if report["status"] != "retained_stage_trace_complete" or report["job_id"] != "21313":
        raise ValueError("Unexpected extraction identity")
    records = [report["source"], *report["inputs"], *report["frozen_snapshots"], report["pair_trace"]]
    for item in records:
        verify_file(Path(item["path"]), item)
    results = root / "benchmark_tools/results"
    references, uncertain, _, _ = load_reference_snapshot(results / "orthobench_paired_uncertainty_20260916.json")
    snapshots = {Path(item["path"]).name: json.loads(Path(item["path"]).read_text()) for item in report["frozen_snapshots"]}
    prepared = snapshots["orthobench_factorial_prepared_20260916.json"]
    factorial = snapshots["orthobench_factorial_results_20260916.json"]
    replay = json.loads(Path(prepared["replay_verification"]["path"]).read_text())
    verify_file(Path(replay["preflight"]["path"]), replay["preflight"])
    preflight = json.loads(Path(replay["preflight"]["path"]).read_text())
    verify_file(Path(preflight["replay_source"]["path"]), preflight["replay_source"])
    if preflight["replay_source"]["sha256"] != "852c3e4fc1a53de7e6046aa78324da587376c0df6a0db1cd8265af65c75bea0f":
        raise ValueError("Replay source differs from reviewed branch semantics")
    sources = {stage: replay["stages"][stage]["output"] for stage in STAGES[:4]}
    sources.update(candidates=factorial["predictions"]["p1_c1_r0"], root_hogs=factorial["predictions"]["p1_c1_r1"])
    owners = {}
    for item in prepared["fasta_inputs"]:
        for protein in SeqIO.parse(item["path"], "fasta"):
            if protein.id in owners:
                raise ValueError("Duplicate FASTA identifier")
            owners[protein.id] = Path(item["path"]).name
    groups, indices = {}, {}
    for stage, item in sources.items():
        verify_file(Path(item["path"]), item)
        groups[stage], indices[stage] = partition(Path(item["path"]), "root_hogs" if stage == "root_hogs" else "plain", set(owners))
        score = score_partition(list(groups[stage].values()), references, uncertain)
        if score != report["scores"][stage]:
            raise ValueError("Fresh full-reference stage scores differ")
    with Path(report["pair_trace"]["path"]).open() as handle:
        rows = verify_pairs(csv.DictReader(handle, delimiter="\t"), references, indices, owners)
    reference_universe = set().union(*references.values())
    for name, pairs in rows.items():
        row = report["families"][name]
        if row["possible_pairs"] != len(pairs):
            raise ValueError("Family pair total differs")
        for stage in STAGES:
            summary = family_group_summary(references[name], groups[stage], indices[stage], reference_universe)
            if summary != row["stages"][stage]:
                raise ValueError("Family group summary differs from native checkpoint")
        for before, after in zip(STAGES, STAGES[1:]):
            counts = Counter("retained" if p[before] and p[after] else "lost" if p[before] else "gained" if p[after] else "absent_both" for p in pairs)
            if {key: counts[key] for key in ("retained", "lost", "gained", "absent_both")} != row["transitions"][before + "_to_" + after]:
                raise ValueError("Original descriptive transition counts differ")
        row["transitions"] = transitions(pairs)
    for item in records:
        verify_file(Path(item["path"]), item)
    report.update(status="native_pair_trace_verified_branch_labels_corrected", original_trace=file_provenance(path),
                  admission_source=file_provenance(Path(__file__)), scheduler=scheduler,
                  admission_helpers=[file_provenance(Path(__file__).with_name(name)) for name in
                      ("trace_ob_families.py", "score_orthobench_partition.py", "score_ygob_groups.py")],
                  transition_specification=TRANSITIONS, stage_sources=sources,
                  lineage_source=preflight["replay_source"], pair_rows_verified=sum(map(len, rows.values())))
    report["limitations"].append("The original extraction ordered checkpoints as a linear chain. This admitted report corrects transitions to the source-defined branches; raw pair membership is unchanged.")
    report["limitations"].append("Cached-hit values retain the hash-verified extraction evidence; this admission independently rechecks native group membership and complete pair coverage, not search execution.")
    output.mkdir(parents=True)
    (output / "results.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    (output / "results.md").write_text(render(report))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    admit(args.root.resolve(), args.output.resolve())


if __name__ == "__main__":
    main()
