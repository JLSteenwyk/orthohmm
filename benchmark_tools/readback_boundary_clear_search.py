"""Reconstruct complete native search outputs and admit only raw paired timings."""

import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_membership_hillclimb import write_json

ARMS = ("full_clear", "boundary_only")
ORDERS = (ARMS, ARMS[::-1], ARMS, ARMS[::-1])


def timing_summary(case):
    if len(case["repeats"]) != 4 or [tuple(row["order"]) for row in case["repeats"]] != list(ORDERS):
        raise ValueError("Incomplete/reordered native search timing schedule")
    for index, row in enumerate(case["repeats"]):
        if row["repeat"] != index or set(row["arms"]) != set(ARMS):
            raise ValueError("Invalid native search arm inventory")
        if any(not math.isfinite(arm[field]) or arm[field] <= 0
               for arm in row["arms"].values() for field in ("wall_s", "process_cpu_s")):
            raise ValueError("Invalid native search cost")
    ratios = {field:statistics.median(row["arms"][ARMS[1]][field]/row["arms"][ARMS[0]][field]
              for row in case["repeats"]) for field in ("wall_s", "process_cpu_s")}
    if ratios != case["ratios"]:
        raise ValueError("Native search ratios do not reproduce from raw repeats")
    return ratios


def run(repo, report_path, output):
    if output.exists() or os.environ.get("CUDA_VISIBLE_DEVICES") != "":
        raise ValueError("Require new CPU-only native-search readback namespace")
    report = json.loads(report_path.read_text())
    scheduler = subprocess.check_output(["sacct", "-j", report["job_id"], "-X",
        "--format=JobID,State,ExitCode,AllocCPUS", "--parsable2"], text=True)
    if scheduler.splitlines() != ["JobID|State|ExitCode|AllocCPUS", report["job_id"] + "|COMPLETED|0:0|4"]:
        raise ValueError("Native search producer is not successfully terminal")
    if (report["status"] != "complete_boundary_clear_native_search_output_equivalence" or len(report["cases"]) != 12
            or any(report[key] for key in ("truth_labels_loaded", "defaults_promoted", "goal_completed"))):
        raise ValueError("Incomplete native search validation or exceeded scope")
    full, boundary = (report["compiler_commands"][arm] for arm in ARMS)
    if (len(full) != 9 or len(boundary) != 10 or full[:6] != boundary[:6]
            or boundary[6] != "-DHMM_VITERBI_BOUNDARY_CLEAR_ONLY=1" or boundary[7] != full[6]):
        raise ValueError("Native search builds differ beyond the planned flag")
    source = Path(full[6]).parents[3]
    checked = [record(report_path), *report["checked"], *report["sources"],
               *report["runtime"]["files"], *report["kernels"].values(), report["compiler"]]
    for pin in checked:
        check(pin)
    for pin in report["sources"]:
        relative = str(Path(pin["path"]).relative_to(source))
        blob = subprocess.check_output(["git", "show", report["source_revision"] + ":" + relative], cwd=repo)
        if hashlib.sha256(blob).hexdigest() != pin["sha256"]:
            raise ValueError("Native source differs from its frozen revision")
    sys.path.insert(0,str(source))
    from orthohmm.search import engine
    from orthohmm.search.sequences import SpeciesSequences
    if Path(engine.__file__).resolve() != source / "orthohmm/search/engine.py":
        raise ValueError("Readback imported the wrong scientific search engine")
    helper_path = source / "benchmark_tools/validate_boundary_clear_search.py"
    spec = importlib.util.spec_from_file_location("frozen_native_search_validation",helper_path)
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    if helper.runtime_records() != report["runtime"] or report["settings"] != helper.SETTINGS:
        raise ValueError("Native search runtime/settings changed")
    live_path = next(Path(pin["path"]) for pin in report["checked"] if Path(pin["path"]).name == "execution_plan.json")
    live = json.loads(live_path.read_text())
    inputs = next(case["inputs"] for case in live["cases"] if case["dataset"] == "QfO" and case["label"] == "control100")
    pins = {Path(pin["path"]).name:pin for pin in inputs}
    names = sorted(pins,key=lambda name:(hashlib.sha256(name.encode()).hexdigest(),name))[:3]
    expected = [(name,names[(index+1)%3],order,cap) for index,name in enumerate(names)
                for order in ("fasta","reversed") for cap in (100,400)]
    if [(case["query"],case["target"],case["query_order"],case["cap"]) for case in report["cases"]] != expected:
        raise ValueError("Native search case inventory differs from the prespecified cohort")
    loaded = {name:SpeciesSequences.from_fasta(pins[name]["path"],name) for name in names}
    libraries = {arm:helper.bind_library(report["kernels"][arm]["path"]) for arm in ARMS}
    rows, canonical = [], {}
    for case in report["cases"]:
        query, indices = helper.subset_species(loaded[case["query"]],512)
        if case["query_order"] == "reversed":
            query = helper.reordered_species(query,range(511,-1,-1))
            indices = indices[::-1]
        target = loaded[case["target"]]
        if query.ids != case["query_ids"] or indices != case["query_indices"] or target.num_sequences != case["target_genes"]:
            raise ValueError("Native search selected gene universe/order changed")
        arrays = (query.flat_sequences,query.offsets,query.lengths,target.flat_sequences,target.offsets,target.lengths)
        if helper.array_records(arrays) != case["input_sequences"]:
            raise ValueError("Native search input sequence arrays changed")
        reference = helper.output_identity(case["warm"][ARMS[0]])
        for arm in ARMS:
            if (helper.output_identity(case["warm"][arm]) != reference
                    or case["warm"][arm]["kernel_inputs"] != case["warm"][ARMS[0]]["kernel_inputs"]):
                raise ValueError("Recorded warm search outputs/inputs differ")
            for repeat in case["repeats"]:
                if helper.output_identity(repeat["arms"][arm]) != reference or repeat["arms"][arm]["kernel_calls"] != 1:
                    raise ValueError("Recorded timed search output/route differs")
            actual = helper.search_with_library(query,target,case["cap"],libraries[arm],capture_inputs=True)
            if helper.output_identity(actual) != reference or actual["kernel_inputs"] != case["warm"][arm]["kernel_inputs"]:
                raise ValueError("Fresh native search does not reproduce all recorded outputs/input arrays")
        if helper.array_records(arrays) != case["input_sequences"]:
            raise ValueError("Fresh native search changed sequence inputs")
        key = (case["query"],case["cap"])
        value = reference["canonical_gene_pair_sha256"]
        if key in canonical and canonical[key] != value:
            raise ValueError("Fresh canonical gene-pair identity differs across query orders")
        canonical[key] = value
        rows.append({"query":case["query"],"target":case["target"],"query_order":case["query_order"],"cap":case["cap"],
            "candidate_records":reference["fingerprint"]["candidate_count"],"significant_hits":reference["significant_hits"],
            "ratios":timing_summary(case),"all_native_outputs_and_inputs_reconstructed":True})
    for pin in checked:
        check(pin)
    if helper.runtime_records() != report["runtime"]:
        raise ValueError("Scientific runtime changed during independent native search replay")
    write_json(output,{"status":"all_boundary_clear_search_outputs_reconstructed_and_timings_verified",
        "producer":record(report_path),"reader":record(__file__),"source_revision":report["source_revision"],
        "scheduler":scheduler,"cases":rows,"candidate_records":sum(row["candidate_records"] for row in rows),
        "all_paired_median_cpu_and_wall_ratios_below_one":all(ratio<1 for row in rows for ratio in row["ratios"].values()),
        "scope":"Fresh replay of both private libraries on every prespecified sampled search; independent raw-timing median calculation. Warm shared-host species-pair search costs, not cold-process/full-pipeline timings or biological accuracy. Candidate records sum case instances including overlapping caps and query orders, not unique predicted ortholog pairs.",
        "defaults_promoted":False,"goal_completed":False})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo",required=True,type=Path)
    parser.add_argument("--report",required=True,type=Path)
    parser.add_argument("--output",required=True,type=Path)
    args = parser.parse_args()
    run(args.repo.resolve(),args.report.resolve(),args.output.resolve())
