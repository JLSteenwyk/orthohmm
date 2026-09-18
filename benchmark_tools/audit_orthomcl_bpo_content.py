"""Check every converted BPO record against its contiguous source BLAST HSPs."""

import argparse
from itertools import groupby, zip_longest
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.convert_orthomcl_blast import read_fasta_lengths
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


def source_rows(stream, lengths):
    for number, line in enumerate(stream, 1):
        fields = line.rstrip(b"\r\n").split(b"\t")
        if len(fields) != 12 or any(identifier not in lengths for identifier in fields[:2]):
            raise ValueError(f"Malformed/unknown source BLAST record at line {number}")
        yield fields


def expected_records(stream, lengths, counts):
    seen_queries, seen_subjects, current_query = set(), set(), None
    serial = 0
    for (query, subject), hsps in groupby(source_rows(stream, lengths), key=lambda row: tuple(row[:2])):
        if query != current_query:
            if query in seen_queries:
                raise ValueError("Noncontiguous source query block")
            seen_queries.add(query)
            seen_subjects = set()
            current_query = query
        if subject in seen_subjects:
            raise ValueError("Noncontiguous source subject block")
        seen_subjects.add(subject)
        counts["source_pair_blocks"] += 1
        spans, weighted_identity, total_subject_span, first_evalue = [], 0.0, 0, None
        for hsp_id, fields in enumerate(hsps, 1):
            counts["source_hsp_rows"] += 1
            aligned, mismatches, openings, qs, qe, ss, se = map(int, fields[3:10])
            literal = b"1" + fields[10] if fields[10].startswith(b"e-") else fields[10]
            evalue, identity, bits = float(literal), float(fields[2]), float(fields[11])
            if (not all(math.isfinite(v) for v in (evalue, identity, bits)) or evalue < 0 or bits < 0
                    or aligned <= 0 or not 0 <= mismatches <= aligned
                    or not 1 <= qs <= qe <= lengths[query] or not 1 <= ss <= se <= lengths[subject]):
                raise ValueError("Invalid source HSP scores/counts/coordinates")
            query_span, subject_span = qe - qs + 1, se - ss + 1
            gap_residues = (aligned - query_span) + (aligned - subject_span)
            identities = aligned - mismatches - gap_residues
            if (query_span > aligned or subject_span > aligned or identities < 0
                    or not 0 <= openings <= gap_residues or (openings == 0) != (gap_residues == 0)
                    or not 0 <= identity <= 100 or abs(identity - identities / aligned * 100) > .011):
                raise ValueError("Inconsistent source HSP alignment accounting")
            if first_evalue is None:
                first_evalue = literal
            weighted_identity += (identities / aligned * 100) * subject_span
            total_subject_span += subject_span
            spans.append(str(hsp_id).encode() + b":" + fields[6] + b"-" + fields[7]
                         + b":" + fields[8] + b"-" + fields[9] + b".")
        if float(first_evalue) > 1e-5:
            counts["cutoff_excluded_pair_blocks"] += 1
            continue
        serial += 1
        yield [str(serial).encode(), query, str(lengths[query]).encode(), subject,
               str(lengths[subject]).encode(), first_evalue,
               str(int(weighted_identity / total_subject_span)).encode(), b"".join(spans)]


def validate(blast, fasta, bpo):
    lengths = read_fasta_lengths(fasta)
    if not lengths or any(value <= 0 for value in lengths.values()):
        raise ValueError("Require nonempty FASTA sequences")
    counts = {"source_hsp_rows": 0, "source_pair_blocks": 0, "cutoff_excluded_pair_blocks": 0,
              "bpo_pair_records": 0, "input_proteins": len(lengths)}
    with blast.open("rb") as source, bpo.open("rb") as converted:
        expected = expected_records(source, lengths, counts)
        for number, (fields, line) in enumerate(zip_longest(expected, converted), 1):
            if fields is None or line is None:
                raise ValueError(f"BPO/source record counts differ at record {number}")
            observed = line.rstrip(b"\r\n").split(b";")
            if observed != fields:
                raise ValueError(f"BPO fields differ from source HSPs at record {number}")
            counts["bpo_pair_records"] += 1
    if counts["bpo_pair_records"] == 0:
        raise ValueError("Empty retained BPO")
    return counts


def audit(blast, fasta, bpo, output):
    if output.exists():
        raise FileExistsError(output)
    checked = [record(p) for p in (blast, fasta, bpo, Path(__file__),
                                   Path(__file__).with_name("convert_orthomcl_blast.py"))]
    content = validate(blast, fasta, bpo)
    for item in checked:
        check(item)
    report = {"status": "bpo_content_matches_source_hsps", "content": content, "checked_records": checked,
              "accuracy_admitted": False, "publication_ready": False,
              "limitations": [
                  "Complete source-to-BPO field verification, not terminal execution or biological validation.",
                  "Recomputes the reviewed legacy first-HSP E-value and subject-span-weighted identity semantics.",
                  "Uses the same binary floating arithmetic convention as the native-compatible converter; not a proof against shared algorithmic errors.",
                  "The independent native BioPerl fixture comparison is complementary, not a full native parse of production inputs.",
                  "Search provenance/query diagnostics, native index validation and final groups remain separate gates."]}
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("blast", "fasta", "bpo", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.blast.resolve(), args.fasta.resolve(), args.bpo.resolve(), args.output.resolve())
    print(json.dumps({"status": result["status"], "content": result["content"]}))
