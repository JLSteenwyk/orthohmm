"""Stage corrected QfO FASTAs only after explicitly pinned compatibility audits."""

import argparse
import gzip
import hashlib
import json
from pathlib import Path, PurePosixPath
import sys
import tarfile

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_qfo_archive_sequences import DATABASE_SHA, classify_matches
from benchmark_tools.inventory_swiss_sequences import PREPARED_SHA
from benchmark_tools.resolve_swiss_sequence_aliases import MAPPING_SHA
from benchmark_tools.audit_swiss_missing_input_relations import ALIASES_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check

XENOPUS = "UP000008143_8364.fasta"
REFERENCE_COUNT = 984137


def require(condition, message):
    if not condition:
        raise ValueError(message)


def validate_reports(comparison, sequences, expected_names, missing_accessions):
    require(comparison["status"] == "qfo_canonical_archive_compared" and
            sequences["status"] == "archive_native_sequence_comparison", "Unexpected audit status")
    require(comparison["inputs"][:3] == sequences["inputs"][:3], "Audit source identities disagree")
    require([x["sha256"] for x in comparison["inputs"][1:]] ==
            [PREPARED_SHA, MAPPING_SHA, ALIASES_SHA], "Changed comparison references")
    require([x["sha256"] for x in sequences["inputs"][1:]] ==
            [PREPARED_SHA, MAPPING_SHA, DATABASE_SHA], "Changed sequence references")
    require(comparison["inputs"][0]["bytes"] == 2648666198, "Unexpected corrected archive size")
    archive = sequences["archive_comparison"]
    shared_keys = ("canonical_files", "canonical_proteomes", "changed_canonical_files", "sequences",
                   "unique_mapped_numeric_ids", "unmapped_accessions", "mapping_numeric_identity_count",
                   "mapping_numeric_ids_without_canonical_accession", "missing_numeric_ids_by_species",
                   "gzip_read_to_eof")
    require(all(comparison[k] == archive[k] for k in shared_keys), "Archive audit results disagree")
    files = comparison["canonical_files"]
    require(len(expected_names) == len(files) == comparison["canonical_proteomes"] == 78 and
            set(files) == set(expected_names), "Unexpected canonical inventory")
    require(comparison["changed_canonical_files"] == [XENOPUS] and
            sorted(k for k, v in files.items() if not v["identical_to_original"]) == [XENOPUS],
            "Release changes extend beyond expected Xenopus correction")
    require(comparison["gzip_read_to_eof"] is True and sequences["inputs_modified"] is False,
            "Incomplete archive or modified source")
    require(all(comparison[k] == REFERENCE_COUNT for k in
                ("sequences", "unique_mapped_numeric_ids", "mapping_numeric_identity_count")) and
            comparison["unmapped_accessions"] == 0 and
            not comparison["mapping_numeric_ids_without_canonical_accession"] and
            not comparison["missing_numeric_ids_by_species"] and sequences["complete_mapping_coverage"] is True,
            "Incomplete reference mapping coverage")
    require(sum(r["sequences"] for r in files.values()) == REFERENCE_COUNT and
            all(r["sequences"] == r["mapped_accessions"] > 0 and not r["unmapped_accessions"]
                for r in files.values()), "Per-proteome counts disagree")
    recovered = comparison["recovered_missing_reference_accessions"]
    require(len(missing_accessions) == 14 and set(recovered) == set(missing_accessions) and
            not comparison["remaining_missing_accessions"], "Incomplete SwissTrees accession recovery")
    recovered_numbers = [r["numeric_protein_id"] for r in recovered.values()]
    require(len(set(recovered_numbers)) == 14 and all(type(n) is int and 1 <= n <= REFERENCE_COUNT
                                                   for n in recovered_numbers), "Invalid recovered identities")
    native = sequences["native_sequence_comparison"]
    require(native["native_entries"] == native["mapped_input_sequences"] == REFERENCE_COUNT and
            not native["native_entries_without_input_by_species"], "Incomplete native database coverage")
    differences = native["differences"]
    ids = [r["numeric_id"] for r in differences]
    require(len(ids) == len(set(ids)) == native["sequence_different"] and
            all(type(n) is int and 1 <= n <= REFERENCE_COUNT for n in ids), "Invalid native difference IDs")
    classes = classify_matches(native)
    require(classes == sequences["sequence_match_classes"] and
            not classes["unexplained_sequence_difference_ids"], "Unexplained or inconsistent sequence differences")
    return files


def extract_canonical(archive, files, destination):
    require(all(Path(name).name == name and name.endswith(".fasta") and name not in (".", "..")
                for name in files), "Unsafe output filename")
    destination.mkdir(parents=True, exist_ok=False)
    seen, outputs = set(), []
    # Never extract archive paths. Copy only approved regular members to flat, exclusive files.
    with gzip.open(archive, "rb") as compressed:
        with tarfile.open(fileobj=compressed, mode="r|") as stream:
            for member in stream:
                name = PurePosixPath(member.name).name
                if not name.endswith(".fasta") or name.endswith(("_DNA.fasta", "_additional.fasta")):
                    continue
                require(name in files and name not in seen and member.isfile(), "Unexpected or duplicate canonical member")
                expected = files[name]
                require(member.name == expected["archive_member"] and member.size == expected["bytes"],
                        "Archive member identity changed")
                require(not PurePosixPath(member.name).is_absolute() and ".." not in PurePosixPath(member.name).parts,
                        "Unsafe archive member path")
                digest, count = hashlib.sha256(), 0
                with stream.extractfile(member) as source, (destination / name).open("xb") as output:
                    for chunk in iter(lambda: source.read(1024 * 1024), b""):
                        digest.update(chunk)
                        count += len(chunk)
                        output.write(chunk)
                require(count == expected["bytes"] and digest.hexdigest() == expected["sha256"],
                        "Extracted FASTA checksum mismatch")
                outputs.append(record(destination / name))
                seen.add(name)
        for _ in iter(lambda: compressed.read(1024 * 1024), b""):
            pass
    require(seen == set(files), "Missing canonical member")
    return sorted(outputs, key=lambda r: r["path"])


def stage(comparison_path, comparison_sha, sequence_path, sequence_sha, destination):
    reports = [record(comparison_path), record(sequence_path)]
    require([r["sha256"] for r in reports] == [comparison_sha, sequence_sha], "Unreviewed or changed audit report")
    comparison, sequences = [json.loads(p.read_text()) for p in (comparison_path, sequence_path)]
    sources = comparison["inputs"] + sequences["inputs"]
    for identity in sources:
        check(identity)
    prepared = json.loads(Path(comparison["inputs"][1]["path"]).read_text())
    aliases = json.loads(Path(comparison["inputs"][3]["path"]).read_text())
    files = validate_reports(comparison, sequences,
                             {Path(r["path"]).name for r in prepared["input_fastas"]},
                             set(aliases["summary"]["missing_genes"]))
    outputs = extract_canonical(Path(sources[0]["path"]), files, destination)
    for identity in [*sources, *reports, *outputs]:
        check(identity)
    result = {"status": "corrected_inputs_staged_pending_inference_freeze", "source": record(__file__),
              "audits": reports, "sources": sources, "input_fastas": outputs,
              "sequence_match_classes": sequences["sequence_match_classes"],
              "total_sequences": comparison["sequences"], "inputs_normalized": False,
              "inference_authorized": False,
              "limitations": ["Explicit report hashes require prior human/agent review of scheduler, pinned executors and audit provenance.",
                              "No inference or scoring launched. Staged files require independent inventory and execution-manifest freeze.",
                              "On failure, partial directory is preserved without a success manifest; never overwrite it."]}
    with (destination / "staging_manifest.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("comparison", "sequences", "destination"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("comparison-sha256", "sequences-sha256"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    result = stage(args.comparison, args.comparison_sha256, args.sequences, args.sequences_sha256, args.destination)
    print(result["status"])
