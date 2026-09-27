"""Trace retained Three Kingdoms OrthoMCL input transformations without reruns."""

import argparse
import json
from pathlib import Path

from benchmark_tools.audit_ob_orthomcl_provenance import mapping
from benchmark_tools.audit_three_kingdoms_sources import sequence_inventory
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_publication_pipeline import save

SOURCE_SHA = "a79bf83e1ea28a9790e4597ea50a44d6409f4f36404a35e660d56a0b754fc1f3"


def compare(expected, actual):
    shared = expected.keys() & actual.keys()
    return dict(expected_genes=len(expected), observed_genes=len(actual),
        missing_ids=sorted(expected.keys() - actual.keys()),
        extra_ids=sorted(actual.keys() - expected.keys()),
        changed_sequences=[dict(gene=g, expected=expected[g], observed=actual[g])
                           for g in sorted(shared) if expected[g] != actual[g]],
        content_equal=expected == actual)


def audit(repo):
    source_path = repo / "benchmark_tools/results/three_kingdoms_sources_20260918.json"
    if record(source_path)["sha256"] != SOURCE_SHA:
        raise ValueError("Changed staged-panel source evidence")
    source = json.loads(source_path.read_text())
    root = repo / "three_kingdoms/results/parity_20260907/orthomcl_1_4"
    directory = root / "orthomcl"
    native_files = [directory / "Sep_7/tmp/all.fa", directory / "Sep_7/tmp/all.gg",
                    directory / "Sep_7/orthomcl.log", root / "recovery_run.log",
                    root / "input.sha256"]
    watched = [record(source_path), *map(record, native_files)]
    staged, copied, species, rows = {}, {}, {}, []
    if len(source["inputs"]) != 12:
        raise ValueError("Expected twelve staged proteomes")
    for item in source["inputs"]:
        baseline = item["files"]["staged"]
        check(baseline)
        copy_path = directory / "data" / (item["code"] + ".fa")
        copy_record = record(copy_path)
        watched.extend([baseline, copy_record])
        expected = sequence_inventory(Path(baseline["path"]))
        observed = sequence_inventory(copy_path)
        if staged.keys() & expected.keys() or copied.keys() & observed.keys():
            raise ValueError("Cross-species duplicate gene identifiers")
        staged.update(expected)
        copied.update(observed)
        species[item["code"]] = set(expected)
        rows.append(dict(species=item["code"], staged=baseline, native_copy=copy_record,
                         byte_equal=baseline["sha256"] == copy_record["sha256"],
                         comparison=compare(expected, observed)))
    if len(staged) != 443217:
        raise ValueError("Wrong full Three Kingdoms universe")
    merged = sequence_inventory(native_files[0])
    species_mapping = mapping(native_files[1], species)
    if any(k != v["input_species"] for k, v in species_mapping.items()):
        raise ValueError("Native species labels are permuted")
    for item in watched:
        check(item)
    return dict(status="retained_three_kingdoms_orthomcl_input_content_audited",
        per_species=rows, copies_vs_staged=compare(staged, copied),
        merged_vs_staged=compare(staged, merged), merged_vs_copies=compare(copied, merged),
        species_mapping=species_mapping, checked_records=watched,
        source=record(__file__), reader_sources=[record(Path(__file__).with_name(n)) for n in
            ("audit_three_kingdoms_sources.py", "audit_ob_orthomcl_provenance.py")],
        historical_consumption_proven=False, scores_changed=False, publication_ready=False,
        limitations=["Retained sequence-content and species-mapping identity, not immutable execution attestation.",
            "Logs document the original merged search input and recovered genome map; no BLAST database or BPO audit here.",
            "Does not close the historical high-sensitivity input-hash gap or establish proteome-wide accuracy."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.repo.resolve())
    save(args.output.resolve(), result)
    print(json.dumps({k: result[k] for k in ("status", "copies_vs_staged", "merged_vs_staged")}))
