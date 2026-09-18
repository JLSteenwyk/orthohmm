"""Stream-check protein BLAST m8 structure and query coverage before BPO conversion."""

import argparse
from collections import Counter
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_orthomcl_blast import parse_diagnostics
from benchmark_tools.convert_orthomcl_blast import read_fasta_lengths, _numeric_evalue
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


def audit_table(blast, fasta, log):
    lengths = read_fasta_lengths(fasta)
    if not lengths or any(length <= 0 for length in lengths.values()):
        raise ValueError("Require nonempty input sequences")
    diagnostics = parse_diagnostics(log)
    if any(gene.encode() not in lengths for gene in diagnostics):
        raise ValueError("Diagnostic query absent from input")
    queries, subjects, self_hits = set(), set(), set()
    current_query, current_subject = None, None
    query_subjects = set()
    rows = pair_blocks = over_cutoff = 0
    with blast.open("rb") as stream:
        for number, line in enumerate(stream, 1):
            fields = line.rstrip(b"\r\n").split(b"\t")
            if len(fields) != 12:
                raise ValueError(f"Expected 12 BLAST columns at line {number}")
            query, subject = fields[:2]
            if query not in lengths or subject not in lengths:
                raise ValueError(f"Unknown BLAST identifier at line {number}")
            try:
                identity = float(fields[2])
                aligned, mismatches, openings, qstart, qend, sstart, send = map(int, fields[3:10])
                evalue, bits = _numeric_evalue(fields[10]), float(fields[11])
            except (ValueError, OverflowError) as error:
                raise ValueError(f"Invalid BLAST number at line {number}") from error
            if (not all(math.isfinite(v) for v in (identity, evalue, bits))
                    or not 0 <= identity <= 100 or evalue < 0 or bits < 0
                    or aligned <= 0 or not 0 <= mismatches <= aligned
                    or not 0 <= openings <= aligned):
                raise ValueError(f"Invalid BLAST value range at line {number}")
            if not (1 <= qstart <= qend <= lengths[query] and 1 <= sstart <= send <= lengths[subject]):
                raise ValueError(f"Invalid protein alignment coordinates at line {number}")
            qspan, sspan = qend - qstart + 1, send - sstart + 1
            gaps = 2 * aligned - qspan - sspan
            matches = aligned - mismatches - gaps
            if (qspan > aligned or sspan > aligned or matches < 0
                    or openings > gaps or (openings == 0) != (gaps == 0)
                    or abs(identity - 100 * matches / aligned) > .011):
                raise ValueError(f"Inconsistent alignment accounting at line {number}")
            # BPO conversion combines adjacent HSPs; revisiting a closed block
            # would silently create multiple entries for the same sequence pair.
            if query != current_query:
                if query in queries:
                    raise ValueError(f"Noncontiguous query block at line {number}")
                current_query, current_subject = query, None
                query_subjects = set()
                queries.add(query)
            if subject != current_subject:
                if subject in query_subjects:
                    raise ValueError(f"Noncontiguous subject block at line {number}")
                query_subjects.add(subject)
                current_subject = subject
                pair_blocks += 1
            subjects.add(subject)
            if query == subject:
                self_hits.add(query)
            rows += 1
            over_cutoff += evalue > 1e-5
    if rows == 0:
        raise ValueError("Empty BLAST hit table")
    failed = {gene.encode() for gene, value in diagnostics.items() if value["query_failed"]}
    absent = lengths.keys() - queries
    details = []
    for gene, value in sorted(diagnostics.items()):
        key = gene.encode()
        details.append({**value, "length": lengths[key], "has_query_hits": key in queries,
                        "has_subject_hits": key in subjects, "has_self_hit": key in self_hits})
    return {"input_proteins": len(lengths), "hsp_rows": rows, "distinct_directed_pairs": pair_blocks,
            "queries_with_hits": len(queries), "subjects_with_hits": len(subjects),
            "proteins_with_self_hits": len(self_hits), "hsp_rows_above_1e_minus_5": over_cutoff,
            "failed_queries": len(failed), "failed_queries_with_query_hits": len(failed & queries),
            "failed_queries_with_subject_hits": len(failed & subjects),
            "queries_without_hits": len(absent), "queries_without_hits_and_without_logged_failure": len(absent - failed),
            "query_ids_without_hits": [gene.decode() for gene in sorted(absent)],
            "diagnostic_lines_by_category": dict(Counter(
                message["category"] for value in diagnostics.values() for message in value["messages"])),
            "diagnostics": details}


def audit(blast, fasta, log, output):
    if output.exists():
        raise FileExistsError(output)
    checked = [record(path) for path in (blast, fasta, log, Path(__file__),
        Path(__file__).with_name("audit_orthomcl_blast.py"),
        Path(__file__).with_name("convert_orthomcl_blast.py"))]
    content = audit_table(blast, fasta, log)
    for item in checked:
        check(item)
    report = {"status": "blast_table_structure_and_query_coverage_verified",
              "checked_records": checked, "content": content, "search_admitted": False,
              "accuracy_admitted": False, "publication_ready": False,
              "limitations": [
                  "This checks a hit table and diagnostic log, not terminal execution or database identity.",
                  "No-hit queries are not automatically failed searches; masking or statistical thresholds can suppress hits.",
                  "Incoming subject hits are distinguished from successful outgoing query searches.",
                  "HSPs above the cutoff are counted, not silently removed or interpreted as native failure.",
                  "BPO/index validation, final clustering, failure impact and reference scoring remain separate checks."]}
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("blast", "fasta", "log", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.blast.resolve(), args.fasta.resolve(), args.log.resolve(), args.output.resolve())
    print(json.dumps({"status": result["status"], "hsp_rows": result["content"]["hsp_rows"]}))
