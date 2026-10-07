"""Independently derive original-protein bridges and all changed group paths."""

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
from pathlib import Path
import sqlite3
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools import readback_native_qfo_candidate_groups as graph
from benchmark_tools import readback_native_qfo_vgnc_blocks as evidence

ROOT = Path(__file__).resolve().parent.parent
HEADER = ["protein_left", "protein_right", "native_gene_left", "native_gene_right", "baseline_state", "candidate_state",
    "baseline_group_left", "baseline_group_right", "candidate_group", "first_connected_iteration",
    "first_direct_event", "connection_path"]


def bridges(mapping, sql_rows, native, missing):
    evidence.need(missing and len(set(missing)) == len(missing) and not set(missing).intersection(native), "Invalid missing cohort")
    numbers = {mapping["mapping"].get(a) for a in missing}
    evidence.need(all(type(n) is int and n > 0 for n in numbers), "Invalid original protein number")
    offsets = mapping["Goff"]
    evidence.need(offsets == sorted(set(offsets)) and all(type(x) is int for x in offsets), "Invalid original offsets")
    reverse = defaultdict(list)
    for accession, number in mapping["mapping"].items():
        if type(number) is int and number in numbers and accession in native:
            reverse[number].append(accession)
    result, targets = [], set()
    for scored in sorted(missing):
        number = mapping["mapping"][scored]
        position = sum(x <= number - 1 for x in offsets) - 1
        evidence.need(0 <= position < len(mapping["species"]), "Invalid species position")
        species = mapping["species"][position]
        rows = [r for r in sql_rows if r[1] == number]
        evidence.need(rows and all(len(r) == 4 and type(r[0]) is int and type(r[1]) is int
            and type(mapping["mapping"].get(r[2])) is int and mapping["mapping"].get(r[2]) == number and r[3] == species for r in rows)
            and scored in [r[2] for r in rows], "Original SQL identity/species differs")
        candidates = reverse[number]
        evidence.need(len(candidates) == 1 and candidates[0] in [r[2] for r in rows], "Ambiguous or absent native protein")
        canonical = candidates[0]
        evidence.need(canonical not in targets, "Repeated bridge target")
        targets.add(canonical)
        result.append(dict(scored_accession=scored, native_accession=canonical, native_gene=native[canonical], prot_nr=number,
            species=species, original_proteome_rows=[dict(zip(("rowid", "prot_nr", "uniprot_id", "species"), r)) for r in rows]))
    return result


def localize(initial, final, trace, changes, proven):
    native = {}
    for members in initial.values():
        for gene in members:
            tokens = gene.split("|")
            accession = tokens[1] if len(tokens) >= 2 and tokens[1] else gene
            evidence.need(accession not in native, "Noninjective native normalization")
            native[accession] = gene
    aliases = {r["scored_accession"]: r["native_accession"] for r in proven}
    evidence.need(len(aliases) == len(proven), "Repeated scored alias")
    converted, seen = [], set()
    for a, b, left, right in changes:
        aa, bb = aliases.get(a, a), aliases.get(b, b)
        pair = (min(aa, bb), max(aa, bb))
        evidence.need(a < b and aa in native and bb in native and aa != bb and pair not in seen,
                      "Unmapped or collapsed changed pair")
        seen.add(pair); converted.append((pair[0], pair[1], left, right))
    ledger, summary = graph.localize(initial, final, trace, converted)
    result = []
    for (a, b, left, right), localized in zip(changes, ledger):
        aa, bb = aliases.get(a, a), aliases.get(b, b)
        ga, gb = localized[4:6]
        if aa > bb:
            ga, gb = gb, ga
        result.append([a, b, native[aa], native[bb], left, right, ga, gb, *localized[6:]])
    return result, summary


def review(path, digest):
    ref = evidence.identity(path)
    evidence.need(ref["sha256"] == digest, "Changed alias report")
    report = evidence.read_json(ref)
    evidence.need(report.get("schema") == "native_qfo_candidate_alias_group_trace_v1"
        and report.get("status") == "original_protein_identity_join_and_all_changed_pairs_localized"
        and all(report.get(k) is True for k in ("whole_candidate_partition_reconstructed", "complete_pair_localization",
            "original_scored_identifiers_preserved", "failed_r1_timing_remains_ineligible"))
        and all(report.get(k) is False for k in ("original_failed_export_retried", "accuracy_rescored", "native_inference_reexecuted",
            "new_scoring_or_admission", "uncertainty_admitted", "publication_ready")), "Changed alias analysis scope")
    checked = report["checked_records"]
    for item in checked:
        evidence.verify(item)
    for key in ("source", "protocol", "diagnosis", "decomposition", "original_failure", "mapping", "candidate_database",
                "native_trace", "scored_pair_table"):
        evidence.need(report[key] in checked, "Unbound artifact: " + key)
    evidence.need(report["source"] == evidence.identity(Path(__file__).with_name("join_native_qfo_candidate_aliases.py"))
        and report["protocol"] == evidence.identity(ROOT / "benchmark_tools/results/NATIVE_QFO_CANDIDATE_ALIAS_GROUP_PROTOCOL_20261006.md")
        and evidence.identity(__file__) in checked and evidence.identity(graph.__file__) in checked
        and evidence.identity(evidence.__file__) in checked, "Source/protocol binding differs")
    for name in ("trace_native_qfo_candidate_groups.py", "replay_native_candidate_trace.py",
                 "diagnose_candidate_trace_variation.py", "link_factorial_scaling_resources.py", "score_ygob_groups.py",
                 "prepare_ob_candidate_neighborhood.py", "probe_dgx_step_separation.py", "validate_native_factorial_outputs.py"):
        evidence.need(evidence.identity(Path(__file__).with_name(name)) in checked, "Unbound accepted-kernel source: " + name)
    evidence.need(evidence.identity(ROOT / "qfo_benchmark/og_to_pairwise.py") in checked, "Unbound original normalization source")
    diagnosis = evidence.read_json(report["diagnosis"])
    failure = evidence.read_json(report["original_failure"])
    decomposition = evidence.read_json(report["decomposition"])
    evidence.need(diagnosis.get("schema") == "native_qfo_candidate_group_id_diagnosis_v1"
        and diagnosis.get("status") == "retained_join_failure_independently_explained"
        and diagnosis.get("whole_candidate_partition_reconstructed") is True
        and all(diagnosis.get(k) is False for k in ("pair_localization_admitted", "primary_export_retried", "identifiers_substituted",
            "accuracy_rescored", "native_inference_reexecuted", "uncertainty_admitted", "publication_ready"))
        and diagnosis.get("failed_r1_timing_remains_ineligible") is True
        and diagnosis["failure"] == report["original_failure"] and diagnosis["decomposition"] == report["decomposition"]
        and diagnosis["methods"] == report["methods"] and diagnosis["native_trace"] == report["native_trace"]
        and diagnosis["scored_pair_table"] == report["scored_pair_table"] == decomposition["transition_table"]
        and all(item in checked for item in diagnosis["checked_records"]), "Prior diagnosis identity/scope differs")
    evidence.need(failure.get("schema") == "native_qfo_candidate_group_trace_failure_v1"
        and failure.get("error") == "Duplicate, unordered or unmapped changed pair" and failure.get("automatic_retry") is False
        and failure["decomposition"] == report["decomposition"] and failure["readback"] == diagnosis["decomposition_readback"]
        and failure["source"] == evidence.identity(Path(__file__).with_name("trace_native_qfo_candidate_groups.py"))
        and not (Path(report["original_failure"]["path"]).parent / "report.json").exists(), "Original failure differs")
    previous_readback = evidence.read_json(diagnosis["decomposition_readback"])
    evidence.need(decomposition.get("schema") == "native_qfo_candidate_vgnc_decomposition_v1"
        and previous_readback.get("status") == "candidate_decomposition_independently_verified"
        and previous_readback["report"] == report["decomposition"], "Prior decomposition differs")
    evidence.need([(m["index"], m["cell"]) for m in report["methods"]] == [(6, "p0_c0_r0"), (8, "p0_c1_r0")], "Wrong cohort")
    admissions = []
    for method, name in zip(report["methods"], ("baseline", "candidate")):
        admission = evidence.read_json(method["admission"]); admissions.append(admission)
        conversion = evidence.read_json(method["conversion"])
        evidence.need(method["admission"] == decomposition[name]["admission"]
            and admission.get("schema") == "full_native_factorial_qfo_admission_v1"
            and admission.get("status") == "full_native_factorial_qfo_assessment_admitted"
            and admission.get("accuracy_admitted") is True and admission.get("publication_ready") is False
            and (admission["native_index"], admission["cell"], admission["participant"])
                == (method["index"], method["cell"], decomposition[name]["participant"])
            and admission["pairs_manifest"] == method["conversion"] and method["conversion"] in admission["checked_records"]
            and conversion.get("schema") == "full_native_factorial_qfo_conversion_v1"
            and conversion.get("status") == "full_native_factorial_qfo_pairs_prepared_unscored" and conversion.get("conversion_kind") == "group"
            and (conversion["native_index"], conversion["cell"], conversion["native_job_id"])
                == (method["index"], method["cell"], admission["native_job_id"])
            and conversion["native_input"] == method["partition"] and conversion["terminal_review"] == method["terminal_review"]
            and conversion["mapping"] == report["mapping"] and report["mapping"] in conversion["checked_records"],
            "Original mapping/admission/conversion differs")
    candidate = decomposition["candidate"]
    execution = evidence.read_json(candidate["execution"])
    evidence.need(candidate["database"] == report["candidate_database"]
        and admissions[1]["execution_report"] == candidate["execution"] and candidate["execution"] in admissions[1]["checked_records"]
        and report["candidate_database"] in admissions[1]["checked_records"] and report["candidate_database"] in execution["outputs"]
        and execution.get("status") == "process_succeeded_pending_independent_admission"
        and (execution["native_index"], execution["cell"], execution["exit_code"]) == (8, "p0_c1_r0", 0), "Original DB inventory differs")
    for name in ("map_relations.py", "vgnc_benchmark.py"):
        source = evidence.identity(ROOT / "qfo_benchmark/benchmark-webservice" / name)
        evidence.need(source in checked and source in admissions[1]["checked_records"], "Original scoring source differs")
    universe, initial = graph.partitions(report["methods"][0]["partition"]["path"])
    final_universe, final = graph.partitions(report["methods"][1]["partition"]["path"])
    evidence.need(universe == final_universe and len(universe) == report["genes"] == diagnosis["genes"]
        and len(initial) == report["baseline_groups"] == diagnosis["baseline_groups"]
        and len(final) == report["candidate_groups"] == diagnosis["candidate_groups"], "Whole partition counts differ")
    native = {}
    for gene in universe:
        tokens = gene.split("|")
        accession = tokens[1] if len(tokens) >= 2 and tokens[1] else gene
        evidence.need(accession not in native, "Noninjective native accessions")
        native[accession] = gene
    with gzip.open(report["mapping"]["path"], "rt") as stream:
        mapping = json.load(stream)
    missing = diagnosis["unmapped_accessions"]
    selected_numbers = {mapping["mapping"].get(a) for a in missing}
    evidence.need(selected_numbers and all(type(n) is int and n > 0 for n in selected_numbers), "Invalid missing accession mapping")
    numbers = sorted(selected_numbers)
    with sqlite3.connect(Path(report["candidate_database"]["path"]).resolve().as_uri() + "?mode=ro", uri=True) as connection:
        selected = [list(r) for r in connection.execute("SELECT rowid, prot_nr, uniprot_id, species FROM proteomes WHERE prot_nr IN ("
            + ",".join("?" for _ in numbers) + ") ORDER BY rowid", numbers)]
    proven = bridges(mapping, selected, native, missing)
    evidence.need(selected == report["original_database_rows"] and proven == report["bridges"], "Bridge proof differs")
    changes, counts, seen = [], Counter(), set()
    with Path(report["scored_pair_table"]["path"]).open(newline="") as stream:
        rows = csv.reader(stream, delimiter="\t")
        evidence.need(next(rows) == ["protein_left", "protein_right", "block_left", "block_right", "p0_c0_r0", "p0_c1_r0"], "Wrong transition header")
        for row in rows:
            evidence.need(len(row) == 6 and row[0] < row[1] and tuple(row[:2]) not in seen, "Bad transition row")
            seen.add(tuple(row[:2])); counts[row[4], row[5]] += 1
            if row[4] != row[5]:
                changes.append((row[0], row[1], row[4], row[5]))
    evidence.need(len(seen) == diagnosis["union_scored_pairs"] == decomposition["union_scored_pairs"]
        and len(changes) == report["changed_pairs"] == diagnosis["changed_pairs"]
        and decomposition["transition_counts"] == [dict(baseline=a, candidate=b, pairs=n) for (a, b), n in sorted(counts.items())]
        and sum(any(g not in native for g in c[:2]) for c in changes) == diagnosis["pairs_with_unmapped_accessions"], "Changed-pair inventory differs")
    trace = evidence.read_json(report["native_trace"])
    evidence.need(len(trace) == report["accepted_merges"] == diagnosis["accepted_merges"], "Merge count differs")
    ledger, summary = localize(initial, final, trace, changes, proven)
    evidence.need(summary == report["localized_summary"], "Connection-path summary differs")
    evidence.table(report["ledger"], HEADER, ledger)
    for item in checked:
        evidence.verify(item)
    evidence.verify(ref)
    return dict(schema="native_qfo_candidate_alias_group_readback_v1", status="original_protein_identity_join_and_all_paths_independently_verified",
        source=evidence.identity(__file__), graph_kernel_source=evidence.identity(graph.__file__), identity_kernel_source=evidence.identity(evidence.__file__),
        report=ref, bridges=proven, genes=len(universe), baseline_groups=len(initial), candidate_groups=len(final), accepted_merges=len(trace),
        changed_pairs=len(ledger), localized_summary=summary, checked_input_records=len(checked),
        exporter_or_primary_replay_imported=False, original_scored_identifiers_preserved=True, complete_pair_localization=True,
        whole_candidate_partition_reconstructed=True, original_failed_export_retried=False, accuracy_rescored=False,
        native_inference_reexecuted=False, new_scoring_or_admission=False, uncertainty_admitted=False,
        failed_r1_timing_remains_ineligible=True, publication_ready=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    evidence.need(not args.output.exists() and not args.output.is_symlink(), "Require fresh output")
    result = review(args.report, args.report_sha256)
    with args.output.open("x") as stream:
        json.dump(result, stream, sort_keys=True, indent=2, allow_nan=False); stream.write("\n")
    print(json.dumps({k: result[k] for k in ("status", "changed_pairs", "localized_summary")}))
