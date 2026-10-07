"""Diagnose retained identifier-join failure, without retrying localization."""

import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools import readback_native_qfo_candidate_groups as graph
from benchmark_tools import readback_native_qfo_vgnc_blocks as evidence

ROOT = Path(__file__).resolve().parent.parent
PROTOCOL = ROOT / "benchmark_tools/results/NATIVE_QFO_CANDIDATE_GROUP_ID_DIAGNOSIS_PROTOCOL_20261006.md"


def collect(failure_ref):
    checked = []
    def load(ref):
        checked.append(ref)
        return evidence.read_json(ref)
    failure = load(failure_ref)
    evidence.need(failure.get("schema") == "native_qfo_candidate_group_trace_failure_v1"
        and failure.get("error_type") == "ValueError"
        and failure.get("error") == "Duplicate, unordered or unmapped changed pair"
        and all(failure.get(k) is False for k in ("automatic_retry", "accuracy_rescored",
            "native_inference_reexecuted", "publication_ready")), "Require the retained unresolved join failure")
    directory = Path(failure_ref["path"]).parent
    evidence.need(not (directory / "report.json").exists() and not (directory / "changed_pair_groups.tsv").exists(),
                  "Failed attempt has an unexpected report/ledger")
    evidence.verify(failure["source"])
    evidence.need(failure["source"] == evidence.identity(Path(__file__).with_name("trace_native_qfo_candidate_groups.py")),
                  "Failure primary source differs")
    sources = [failure["source"], *[evidence.identity(p) for p in
        (Path(__file__), Path(graph.__file__), Path(evidence.__file__), PROTOCOL)]]
    checked.extend(sources)
    decomposition, readback = load(failure["decomposition"]), load(failure["readback"])
    evidence.need(decomposition.get("schema") == "native_qfo_candidate_vgnc_decomposition_v1"
        and decomposition.get("status") == "new_candidate_scored_rows_decomposed"
        and readback.get("schema") == "native_qfo_candidate_vgnc_readback_v1"
        and readback.get("status") == "candidate_decomposition_independently_verified"
        and readback.get("report") == failure["decomposition"]
        and all(r.get(k) is False for r in (decomposition, readback)
            for k in ("uncertainty_admitted", "new_scoring_or_admission", "publication_ready"))
        and all(r.get("failed_r1_timing_remains_ineligible") is True for r in (decomposition, readback)),
        "Changed prior decomposition scope/linkage")
    for item in (decomposition["source"], readback["source"]):
        evidence.verify(item); checked.append(item)
    methods = []
    for expected, name in (((6, "p0_c0_r0"), "baseline"), ((8, "p0_c1_r0"), "candidate")):
        admission_ref = decomposition[name]["admission"]
        admission = load(admission_ref)
        evidence.need(admission.get("schema") == "full_native_factorial_qfo_admission_v1"
            and admission.get("status") == "full_native_factorial_qfo_assessment_admitted"
            and admission.get("accuracy_admitted") is True and admission.get("publication_ready") is False
            and (admission["native_index"], admission["cell"]) == expected
            and admission["participant"] == decomposition[name]["participant"]
            and admission["pairs_manifest"] in admission["checked_records"], "Admission identity/inventory differs")
        conversion = load(admission["pairs_manifest"])
        job = admission["native_job_id"]
        evidence.need(conversion.get("schema") == "full_native_factorial_qfo_conversion_v1"
            and conversion.get("status") == "full_native_factorial_qfo_pairs_prepared_unscored"
            and conversion.get("conversion_kind") == "group"
            and (conversion["native_index"], conversion["cell"], conversion["native_job_id"]) == (*expected, job),
            "Conversion identity/semantics differs")
        terminal = load(conversion["terminal_review"])
        evidence.need(terminal.get("schema") == "native_factorial_terminal_review_v1"
            and terminal.get("status") == "native_success" and terminal.get("native_outputs_validated") is True
            and (terminal["index"], terminal["cell"], terminal["job_id"]) == (*expected, job), "Terminal review differs")
        output_ref = terminal["reviews"]["outputs_or_failure"]
        output = load(output_ref)
        evidence.need(output.get("schema") == "native_factorial_output_review_v1"
            and output.get("status") == "native_outputs_validated" and output.get("native_outputs_validated") is True
            and (output["index"], output["cell"], output["job_id"]) == (*expected, job)
            and conversion["native_input"] in output["checked_files"], "Output inventory differs")
        evidence.verify(conversion["native_input"]); checked.append(conversion["native_input"])
        methods.append(dict(index=expected[0], cell=expected[1], admission=admission_ref,
            conversion=admission["pairs_manifest"], terminal_review=conversion["terminal_review"],
            output_review=output_ref, partition=conversion["native_input"], input_genes=output["input_genes"],
            files=output["checked_files"]))
    refs = {}
    for key, suffix in (("trace", "/phylogeny_candidate_merges.json"), ("metrics", "/metrics.json"),
                        ("checkpoint", "/phylogeny_candidate_superfamilies.txt")):
        matches = {json.dumps(r, sort_keys=True): r for r in methods[1]["files"] if r["path"].endswith(suffix)}
        evidence.need(len(matches) == 1, "Ambiguous or missing original artifact")
        refs[key] = next(iter(matches.values()))
        evidence.verify(refs[key]); checked.append(refs[key])
    evidence.need(all(refs["checkpoint"][k] == methods[1]["partition"][k] for k in ("bytes", "sha256")),
                  "Candidate checkpoint differs")
    universe, initial = graph.partitions(methods[0]["partition"]["path"])
    final_universe, final = graph.partitions(methods[1]["partition"]["path"])
    trace, metrics = load(refs["trace"]), load(refs["metrics"])
    profile = metrics["metadata"]["phylogeny_candidate_profile"]
    evidence.need(metrics.get("status") == "complete" and profile["profile"] == "satellite_v2"
        and universe == final_universe and len(universe) == methods[0]["input_genes"] == methods[1]["input_genes"]
        and len(initial) == profile["seed_families"] and len(final) == profile["candidate_families"]
        and len(trace) == profile["merges"], "Complete universe/counts differ")
    reconstructed, _ = graph.traverse(initial, trace, set())
    evidence.need(reconstructed == final, "Complete candidate graph reconstruction differs")
    aliases = {}
    for gene in universe:
        tokens = gene.split("|")
        accession = tokens[1] if len(tokens) >= 2 and tokens[1] else gene
        evidence.need(accession not in aliases, "Noninjective native accession normalization")
        aliases[accession] = gene
    pair_ref = decomposition["transition_table"]
    evidence.verify(pair_ref); checked.append(pair_ref)
    rows_seen, transitions, unresolved, absent, affected = set(), Counter(), [], set(), Counter()
    changed = 0
    with Path(pair_ref["path"]).open(newline="") as stream:
        rows = csv.reader(stream, delimiter="\t")
        evidence.need(next(rows) == ["protein_left", "protein_right", "block_left", "block_right", "p0_c0_r0", "p0_c1_r0"],
                      "Unexpected scored transition header")
        for row in rows:
            evidence.need(len(row) == 6 and row[0] < row[1] and tuple(row[:2]) not in rows_seen,
                          "Malformed, duplicate or unordered transition row")
            rows_seen.add(tuple(row[:2])); transitions[row[4], row[5]] += 1
            if row[4] == row[5]:
                continue
            changed += 1
            evidence.need((row[4], row[5]) in (("FN", "TP"), ("not_scored", "FP")), "Unexpected changed state")
            missing = [gene for gene in row[:2] if gene not in aliases]
            if missing:
                unresolved.append([row[0], row[1], row[4], row[5], ",".join(missing)])
                absent.update(missing); affected[row[5]] += 1
    evidence.need(len(rows_seen) == decomposition["union_scored_pairs"] and decomposition["transition_counts"]
        == [dict(baseline=a, candidate=b, pairs=n) for (a, b), n in sorted(transitions.items())], "Prior pair inventory differs")
    evidence.need(unresolved, "Unmapped identifiers do not explain the retained failure")
    for method in methods:
        del method["files"]
    for item in checked:
        evidence.verify(item)
    return dict(schema="native_qfo_candidate_group_id_diagnosis_v1", status="retained_join_failure_independently_explained",
        source=evidence.identity(__file__), protocol=evidence.identity(PROTOCOL), failure=failure_ref,
        decomposition=failure["decomposition"], decomposition_readback=failure["readback"], methods=methods,
        native_trace=refs["trace"], metrics=refs["metrics"], candidate_checkpoint=refs["checkpoint"], scored_pair_table=pair_ref,
        checked_records=checked, genes=len(universe), baseline_groups=len(initial), candidate_groups=len(final),
        accepted_merges=len(trace), whole_candidate_partition_reconstructed=True, union_scored_pairs=len(rows_seen),
        changed_pairs=changed, pairs_with_unmapped_accessions=len(unresolved), changed_pairs_without_missing_accessions=changed-len(unresolved),
        unmapped_accessions=sorted(absent), affected_pairs_by_candidate_state=dict(sorted(affected.items())),
        pair_localization_admitted=False, primary_export_retried=False, identifiers_substituted=False,
        primary_or_accepted_replay_imported=False, accuracy_rescored=False, native_inference_reexecuted=False,
        uncertainty_admitted=False, failed_r1_timing_remains_ineligible=True, publication_ready=False,
        limitations=["Complete group replay is software compatibility, not biological truth or identical upstream hits.",
            "All changed pairs are inventoried; no approximate alias join or mapped-subset pair localization is admitted.",
            "The reason scored accessions differ from native identifiers remains unresolved; this is not an FP mechanism.",
            "Prior raw scoring/admission is reused, not repeated; current-byte checks do not prove continuous integrity."]), unresolved


def execute(failure_ref, output):
    output = Path(output).absolute()
    evidence.need(output.resolve() == output and not output.exists() and not output.is_symlink(), "Require fresh direct output")
    output.mkdir(parents=True, exist_ok=False)
    try:
        result, rows = collect(failure_ref)
        path = output / "unmapped_pairs.tsv"
        with path.open("x", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
            writer.writerow(["protein_left", "protein_right", "baseline_state", "candidate_state", "unmapped_accessions"])
            writer.writerows(rows)
        result["unmapped_pair_table"] = evidence.identity(path)
        with (output / "report.json").open("x") as stream:
            json.dump(result, stream, sort_keys=True, indent=2, allow_nan=False); stream.write("\n")
        return result
    except Exception as error:
        with (output / "failure.json").open("x") as stream:
            json.dump(dict(schema="native_qfo_candidate_group_id_diagnosis_failure_v1", source=evidence.identity(__file__),
                original_failure=failure_ref, error_type=type(error).__name__, error=str(error), automatic_retry=False,
                pair_localization_admitted=False, publication_ready=False), stream, sort_keys=True, indent=2)
            stream.write("\n")
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--failure", type=Path, required=True)
    parser.add_argument("--failure-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    ref = evidence.identity(args.failure)
    evidence.need(ref["sha256"] == args.failure_sha256, "Selected failure changed")
    result = execute(ref, args.output)
    print(json.dumps({k: result[k] for k in ("status", "genes", "baseline_groups", "candidate_groups", "accepted_merges",
        "changed_pairs", "pairs_with_unmapped_accessions", "unmapped_accessions", "affected_pairs_by_candidate_state")}))
