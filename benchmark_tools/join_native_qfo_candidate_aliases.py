"""Join retained scored aliases by proven original QfO protein identity."""

import argparse
from bisect import bisect_right
from collections import Counter
import csv
import gzip
import json
from pathlib import Path
import sqlite3
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools import trace_native_qfo_candidate_groups as original
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.validate_native_factorial_outputs import require
from qfo_benchmark.og_to_pairwise import _strip_to_uniprot

ROOT = Path(__file__).resolve().parent.parent
PROTOCOL = ROOT / "benchmark_tools/results/NATIVE_QFO_CANDIDATE_ALIAS_GROUP_PROTOCOL_20261006.md"
HEADER = ["protein_left", "protein_right", "native_gene_left", "native_gene_right", "baseline_state", "candidate_state",
    "baseline_group_left", "baseline_group_right", "candidate_group", "first_connected_iteration",
    "first_direct_event", "connection_path"]
DIAGNOSIS_SHA = "3bc6646e1dc7cb6488ceb6aaaf763e3a4e3cdcbd41268f2066ef2d8008f5be74"


def bridge(mapping, sql_rows, native, missing):
    accessions = mapping["mapping"]
    require(missing and len(set(missing)) == len(missing) and not set(missing) & set(native), "Invalid absent accession cohort")
    offsets, species = mapping["Goff"], mapping["species"]
    require(offsets == sorted(set(offsets)) and all(type(x) is int for x in offsets), "Invalid species offsets")
    result, targets = [], set()
    for scored in sorted(missing):
        protein = accessions.get(scored)
        require(type(protein) is int and protein > 0, "Missing or invalid original scored protein number")
        selected = [r for r in sql_rows if r[1] == protein]
        require(selected and all(len(r) == 4 and type(r[0]) is int and type(r[1]) is int
            and type(accessions.get(r[2])) is int and accessions.get(r[2]) == protein for r in selected), "SQL rows disagree with original mapping")
        position = bisect_right(offsets, protein - 1) - 1
        require(0 <= position < len(species), "Protein outside species offsets")
        expected_species = species[position]
        require(all(r[3] == expected_species for r in selected) and scored in {r[2] for r in selected},
                "Scored alias absent from original DB or inconsistent species")
        candidates = [a for a in native if type(accessions.get(a)) is int and accessions.get(a) == protein]
        require(len(candidates) == 1 and candidates[0] in {r[2] for r in selected}, "Require one original native accession for protein")
        target = candidates[0]
        require(target not in targets, "Alias bridge targets are not unique")
        targets.add(target)
        result.append(dict(scored_accession=scored, native_accession=target, native_gene=native[target],
            prot_nr=protein, species=expected_species,
            original_proteome_rows=[dict(rowid=r[0], prot_nr=r[1], uniprot_id=r[2], species=r[3]) for r in selected]))
    return result


def joined_localization(initial, final, trace, changes, bridges):
    native = {_strip_to_uniprot(g): g for group in initial for g in group}
    require(len(native) == sum(map(len, initial)), "Noninjective native normalization")
    aliases = {r["scored_accession"]: r["native_accession"] for r in bridges}
    require(len(aliases) == len(bridges), "Duplicate scored alias")
    converted, seen = [], set()
    for a, b, left, right in changes:
        aa, bb = aliases.get(a, a), aliases.get(b, b)
        pair = tuple(sorted((aa, bb)))
        require(a < b and aa in native and bb in native and aa != bb and pair not in seen,
                "Unmapped, self-collapsed or alias-collapsed changed pair")
        seen.add(pair); converted.append((*pair, left, right))
    ledger, summary, genes = original.localize(initial, final, trace, converted)
    rows = []
    for (a, b, left, right), localized in zip(changes, ledger):
        aa, bb = aliases.get(a, a), aliases.get(b, b)
        groups = localized[4:6] if aa < bb else localized[4:6][::-1]
        rows.append([a, b, native[aa], native[bb], left, right, *groups, *localized[6:]])
    require(len(rows) == len(changes), "Incomplete localization")
    return rows, summary, genes


def collect(diagnosis_ref):
    checked = []
    def load(ref):
        check(ref); checked.append(ref)
        return json.loads(Path(ref["path"]).read_text())
    diagnosis = load(diagnosis_ref)
    require(diagnosis.get("schema") == "native_qfo_candidate_group_id_diagnosis_v1"
        and diagnosis.get("status") == "retained_join_failure_independently_explained"
        and diagnosis.get("whole_candidate_partition_reconstructed") is True
        and all(diagnosis.get(k) is False for k in ("pair_localization_admitted", "primary_export_retried",
            "identifiers_substituted", "accuracy_rescored", "native_inference_reexecuted", "uncertainty_admitted", "publication_ready"))
        and diagnosis.get("failed_r1_timing_remains_ineligible") is True, "Changed prior diagnosis scope")
    for ref in diagnosis["checked_records"]:
        check(ref); checked.append(ref)
    check(diagnosis["unmapped_pair_table"]); checked.append(diagnosis["unmapped_pair_table"])
    failure = load(diagnosis["failure"])
    require(failure.get("schema") == "native_qfo_candidate_group_trace_failure_v1"
        and failure.get("error") == "Duplicate, unordered or unmapped changed pair"
        and failure.get("automatic_retry") is False and failure["source"] == record(original.__file__)
        and failure["decomposition"] == diagnosis["decomposition"] and failure["readback"] == diagnosis["decomposition_readback"]
        and not (Path(diagnosis["failure"]["path"]).parent / "report.json").exists(), "Original failed attempt differs")
    decomposition = load(diagnosis["decomposition"])
    previous_readback = load(diagnosis["decomposition_readback"])
    require(decomposition.get("schema") == "native_qfo_candidate_vgnc_decomposition_v1"
        and previous_readback.get("status") == "candidate_decomposition_independently_verified"
        and previous_readback["report"] == diagnosis["decomposition"]
        and diagnosis["scored_pair_table"] == decomposition["transition_table"], "Decomposition binding differs")
    require([(r["index"], r["cell"]) for r in diagnosis["methods"]] == [(6, "p0_c0_r0"), (8, "p0_c1_r0")], "Wrong cohort")
    mappings, admissions = [], []
    for method, name in zip(diagnosis["methods"], ("baseline", "candidate")):
        admission = load(method["admission"]); admissions.append(admission)
        conversion = load(method["conversion"])
        require(method["admission"] == decomposition[name]["admission"]
            and admission.get("schema") == "full_native_factorial_qfo_admission_v1"
            and admission.get("status") == "full_native_factorial_qfo_assessment_admitted"
            and admission.get("accuracy_admitted") is True and admission.get("publication_ready") is False
            and (admission["native_index"], admission["cell"], admission["participant"])
                == (method["index"], method["cell"], decomposition[name]["participant"])
            and admission["pairs_manifest"] == method["conversion"] and method["conversion"] in admission["checked_records"]
            and conversion.get("schema") == "full_native_factorial_qfo_conversion_v1"
            and conversion.get("status") == "full_native_factorial_qfo_pairs_prepared_unscored"
            and conversion.get("conversion_kind") == "group"
            and (conversion["native_index"], conversion["cell"], conversion["native_job_id"])
                == (method["index"], method["cell"], admission["native_job_id"])
            and conversion["native_input"] == method["partition"]
            and conversion["terminal_review"] == method["terminal_review"]
            and conversion["mapping"] in conversion["checked_records"], "Original admission/conversion/mapping differs")
        mappings.append(conversion["mapping"])
    require(mappings[0] == mappings[1], "Original mappings differ")
    mapping_ref = mappings[0]; check(mapping_ref); checked.append(mapping_ref)
    candidate = decomposition["candidate"]
    execution = load(candidate["execution"])
    require(admissions[1]["execution_report"] == candidate["execution"]
        and candidate["execution"] in admissions[1]["checked_records"]
        and candidate["database"] in admissions[1]["checked_records"] and candidate["database"] in execution["outputs"]
        and execution.get("status") == "process_succeeded_pending_independent_admission"
        and (execution["native_index"], execution["cell"], execution["exit_code"]) == (8, "p0_c1_r0", 0), "Original DB execution differs")
    database_ref = candidate["database"]; check(database_ref); checked.append(database_ref)
    for name in ("map_relations.py", "vgnc_benchmark.py"):
        path = ROOT / "qfo_benchmark/benchmark-webservice" / name
        ref = record(path)
        require(ref in admissions[1]["checked_records"], "Original scoring source not inventoried")
        checked.append(ref)
    source_paths = [Path(__file__), PROTOCOL, Path(original.__file__),
        ROOT / "benchmark_tools/replay_native_candidate_trace.py",
        ROOT / "benchmark_tools/diagnose_candidate_trace_variation.py",
        ROOT / "benchmark_tools/link_factorial_scaling_resources.py",
        ROOT / "benchmark_tools/score_ygob_groups.py", ROOT / "qfo_benchmark/og_to_pairwise.py",
        ROOT / "benchmark_tools/prepare_ob_candidate_neighborhood.py",
        ROOT / "benchmark_tools/probe_dgx_step_separation.py",
        ROOT / "benchmark_tools/validate_native_factorial_outputs.py",
        ROOT / "benchmark_tools/readback_native_qfo_candidate_groups.py",
        ROOT / "benchmark_tools/readback_native_qfo_vgnc_blocks.py",
        ROOT / "benchmark_tools/readback_native_qfo_candidate_aliases.py"]
    checked.extend(record(p) for p in source_paths)
    initial = original.accepted.partition(Path(diagnosis["methods"][0]["partition"]["path"]), "space_separated_groups")
    final = original.accepted.partition(Path(diagnosis["methods"][1]["partition"]["path"]), "space_separated_groups")
    require(initial[0] == final[0] and len(initial[0]) == diagnosis["genes"]
        and len(initial[1]) == diagnosis["baseline_groups"] and len(final[1]) == diagnosis["candidate_groups"], "Partition counts differ")
    native = {_strip_to_uniprot(g): g for g in initial[0]}
    require(len(native) == len(initial[0]), "Noninjective native normalization")
    with gzip.open(mapping_ref["path"], "rt") as stream:
        mapping = json.load(stream)
    missing = diagnosis["unmapped_accessions"]
    selected_numbers = {mapping["mapping"].get(a) for a in missing}
    require(selected_numbers and all(type(n) is int and n > 0 for n in selected_numbers), "Missing original mapped protein")
    numbers = sorted(selected_numbers)
    with sqlite3.connect(Path(database_ref["path"]).resolve().as_uri() + "?mode=ro", uri=True) as connection:
        rows = [list(r) for r in connection.execute("SELECT rowid, prot_nr, uniprot_id, species FROM proteomes WHERE prot_nr IN ("
            + ",".join("?" for _ in numbers) + ") ORDER BY rowid", numbers)]
    bridges = bridge(mapping, rows, native, missing)
    changes, counts, seen = [], Counter(), set()
    with Path(diagnosis["scored_pair_table"]["path"]).open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        require(reader.fieldnames == ["protein_left", "protein_right", "block_left", "block_right", "p0_c0_r0", "p0_c1_r0"], "Wrong transition header")
        for row in reader:
            require(None not in row and None not in row.values(), "Malformed transition row")
            a, b, left, right = (row[k] for k in ("protein_left", "protein_right", "p0_c0_r0", "p0_c1_r0"))
            require(a < b and (a, b) not in seen, "Duplicate or unordered transition pair")
            seen.add((a, b)); counts[left, right] += 1
            if left != right:
                changes.append((a, b, left, right))
    require(len(seen) == decomposition["union_scored_pairs"] == diagnosis["union_scored_pairs"]
        and len(changes) == diagnosis["changed_pairs"]
        and decomposition["transition_counts"] == [dict(baseline=a, candidate=b, pairs=n) for (a, b), n in sorted(counts.items())], "Transition inventory differs")
    trace = load(diagnosis["native_trace"])
    require(len(trace) == diagnosis["accepted_merges"], "Accepted count differs")
    ledger, summary, genes = joined_localization(initial[1], final[1], trace, changes, bridges)
    require(sum(any(g not in native for g in r[:2]) for r in changes) == diagnosis["pairs_with_unmapped_accessions"], "Affected pair count differs")
    for ref in checked:
        check(ref)
    return dict(schema="native_qfo_candidate_alias_group_trace_v1", status="original_protein_identity_join_and_all_changed_pairs_localized",
        source=record(__file__), protocol=record(PROTOCOL), diagnosis=diagnosis_ref, decomposition=diagnosis["decomposition"],
        original_failure=diagnosis["failure"], mapping=mapping_ref, candidate_database=database_ref, original_database_rows=rows,
        bridges=bridges, methods=diagnosis["methods"], native_trace=diagnosis["native_trace"], scored_pair_table=diagnosis["scored_pair_table"],
        checked_records=checked, genes=genes, baseline_groups=len(initial[1]), candidate_groups=len(final[1]),
        accepted_merges=len(trace), changed_pairs=len(ledger), localized_summary=summary,
        whole_candidate_partition_reconstructed=True, complete_pair_localization=True,
        original_scored_identifiers_preserved=True, original_failed_export_retried=False, accuracy_rescored=False,
        native_inference_reexecuted=False, new_scoring_or_admission=False, uncertainty_admitted=False,
        failed_r1_timing_remains_ineligible=True, publication_ready=False,
        limitations=["Explicit joins use original QfO protein identity, not independent biological validation of aliases.",
            "Development-exposed software paths, not true homology/duplication or general accuracy claims.",
            "Accepted-union replay does not recompute candidate eligibility, rejected alternatives or numeric support.",
            "Original failed namespace, scored identifiers, categories, metrics and admissions remain unchanged.",
            "Current-byte checks do not establish uninterrupted integrity or repeat transitive admission."]), ledger


def execute(diagnosis_ref, output):
    output = Path(output).absolute()
    require(output.resolve() == output and not output.exists() and not output.is_symlink(), "Require fresh direct output")
    output.mkdir(parents=True, exist_ok=False)
    try:
        result, ledger = collect(diagnosis_ref)
        path = output / "changed_pair_groups.tsv"
        with path.open("x", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
            writer.writerow(HEADER); writer.writerows(ledger)
        result["ledger"] = record(path)
        save(output / "report.json", result)
        return result
    except Exception as error:
        save(output / "failure.json", dict(schema="native_qfo_candidate_alias_group_trace_failure_v1", source=record(__file__),
            diagnosis=diagnosis_ref, error_type=type(error).__name__, error=str(error), automatic_retry=False,
            original_failed_export_retried=False, accuracy_rescored=False, publication_ready=False))
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagnosis", type=Path, required=True)
    parser.add_argument("--diagnosis-sha256", default=DIAGNOSIS_SHA)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    ref = record(args.diagnosis)
    require(ref["sha256"] == args.diagnosis_sha256, "Selected diagnosis changed")
    result = execute(ref, args.output)
    print(json.dumps({k: result[k] for k in ("status", "genes", "bridges", "changed_pairs", "localized_summary")}))
