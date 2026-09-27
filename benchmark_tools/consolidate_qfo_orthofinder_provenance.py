"""Consolidate retained QfO OrthoFinder provenance without rerunning inference."""

import argparse
import json
from pathlib import Path

from benchmark_tools.audit_ob_orthofinder_provenance import log_command
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.summarize_matched_resources import verbose_time

MANIFEST_SHA = "042aa221d01554114a1b0b413ca3e2ff56dad5f8c06362c69a1bf49f2de642fc"
METHODS = ("orthofinder_3_1_5_full", "orthofinder_3_1_5_sequence_only")


def bind_conversion(row, conversion):
    if (conversion["semantics"] != row["prediction_semantics"]
            or conversion["retained_pairs"] != row["retained_pairs"]
            or conversion["total_pairs"] != row["submitted_pairs"]
            or conversion["removed_mapping_pairs"] != 0
            or conversion["pairs"]["sha256"] != conversion["filtered_pairs"]["sha256"]
            or conversion["finished_epoch"] < conversion["started_epoch"]):
        raise ValueError("Conversion identity, pair count, mapping or time differs")
    return dict(key=row["key"], semantics=row["prediction_semantics"],
        retained_pairs=row["retained_pairs"], scores=row["scores"], secondary_mean=row["secondary_mean"],
        pair_file=conversion["filtered_pairs"], native_admission=conversion["admission"],
        conversion_job=conversion["job_id"],
        conversion_wall_seconds=conversion["finished_epoch"]-conversion["started_epoch"],
        standalone_sequence_only_inference_seconds=None)


def consolidate(repo):
    checked = []
    def read(item):
        check(item)
        checked.append(item)
        return json.loads(Path(item["path"]).read_text())

    source = record(repo / "benchmark_tools/results/qfo_corrected_comparison_20260926_v7/manifest.json")
    if source["sha256"] != MANIFEST_SHA:
        raise ValueError("Changed eight-method comparison")
    manifest = read(source)
    rows, conversions = [], []
    for method in METHODS:
        row = next(r for r in manifest["methods"] if r["key"] == method)
        check(row["admission"])
        checked.append(row["admission"])
        conversion = read(row["conversion"])
        for key in ("pairs", "filtered_pairs", "mapping", "source"):
            check(conversion[key])
            checked.append(conversion[key])
        rows.append(bind_conversion(row, conversion))
        conversions.append(conversion)
    if conversions[0]["admission"] != conversions[1]["admission"]:
        raise ValueError("Methods no longer share native execution")
    native = read(conversions[0]["admission"])
    if native["status"] != "corrected_orthofinder_native_evidence_admitted":
        raise ValueError("Native evidence was not admitted")
    known = {r["path"]:r for r in native["checked_records"]}
    status_path = repo / "benchmarks/results/qfo_corrected_primary_v1/execution/orthofinder_full/status.json"
    execution = read(known[str(status_path)])
    plan = read(execution["manifest"])
    config = plan["methods"]["orthofinder_full"]
    if (execution["native_argv"] != config["native_argv"] or execution["exit_code"] != 0
            or execution["job_id"] != native["scheduler"]["JobIDRaw"]
            or native["scheduler"]["State"] != "COMPLETED" or native["scheduler"]["ExitCode"] != "0:0"):
        raise ValueError("Native execution binding differs")
    for item in (execution["source"], execution["log"], execution["timing"], *plan["inputs"]):
        check(item)
        checked.append(item)
    text = Path(execution["timing"]["path"]).read_text()
    if log_command(text, "Command being timed: ") != config["native_argv"]:
        raise ValueError("Timing command differs")
    resources = verbose_time(text)
    content = native["content"]
    if (content["genes"] != 984137 or content["species"] != 78
            or len(content["copied_inputs"]) != 78
            or len(content["internal_sequences"]) != 78
            or any(r["differences"] or r["sequences"] != r["identical_sequences"]
                   for r in content["internal_sequences"])
            or sum(r["sequences"] for r in content["internal_sequences"]) != 984137
            or content["native_pairs"]["distinct_pairs"] != rows[0]["retained_pairs"]):
        raise ValueError("Retained native coverage differs")
    for item in (*content["copied_inputs"], content["checkpoint"], content["sequence_ids"], content["species_tree"]):
        check(item)
        checked.append(item)
    log = Path(content["results_directory"]) / "Log.txt"
    check(known[str(log)])
    checked.append(known[str(log)])
    if (log_command(log.read_text(), "Command Line: ") != config["native_argv"]
            or "Started OrthoFinder version 3.1.5" not in log.read_text()):
        raise ValueError("Native version/command differs")
    return dict(status="retained_qfo_orthofinder_provenance_consolidated", rows=rows,
        version="3.1.5", native_command=config["native_argv"], native_scheduler=native["scheduler"],
        full_native_resources=resources, input_genes=content["genes"], species=content["species"],
        prephylogeny_checkpoint_groups=content["checkpoint_groups"], native_pairs=content["native_pairs"]["distinct_pairs"],
        checked_records=checked, source=record(__file__), historical_scores_replaced=False, publication_ready=False,
        limitations=["Consolidates retained admissions; does not rerun all native integrity or score arithmetic checks.",
            "Copied inputs are rehashed; internal sequence parity is inherited from the pinned native admission.",
            "Sequence-only pairs come from the full run checkpoint, not an independently timed sequence-only run.",
            "GNU-time full inference and conversion epoch intervals have different scopes; shared-host times are descriptive.",
            "Maximum process RSS is not peak aggregate concurrent process-tree memory.",
            "Native tool runtime evidence is linked through the frozen plan; transitive executable audit is not repeated."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = consolidate(args.repo.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
