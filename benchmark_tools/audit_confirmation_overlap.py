"""Audit candidate protein ZIPs for sequence exposure without reading truth labels."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import io
import json
from pathlib import Path, PurePosixPath
import resource
import sys
import time
from zipfile import ZipFile

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
from benchmark_tools.audit_ygob_overlap import fingerprint
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_membership_hillclimb import write_json

SUFFIXES = frozenset((".fa", ".faa", ".fasta", ".fsa", ".fas"))
MAX_BYTES = 2 ** 31
MAX_MEMBERS = 256


def fasta_members(archive):
    members = archive.infolist()
    names = [member.filename for member in members]
    if len(members) > MAX_MEMBERS or len(names) != len(set(names)):
        raise ValueError("Too many or duplicate ZIP member names")
    if sum(member.file_size for member in members) > MAX_BYTES:
        raise ValueError("ZIP uncompressed size exceeds the audit bound")
    admitted = []
    for member in members:
        path = PurePosixPath(member.filename)
        if (path.is_absolute() or ".." in path.parts or "\\" in member.filename
                or not path.parts or ":" in path.parts[0] or member.flag_bits & 1):
            raise ValueError("Unsafe or encrypted ZIP member")
        if member.is_dir() or "__MACOSX" in path.parts or any(part.startswith(".") for part in path.parts):
            continue
        if path.suffix.lower() in SUFFIXES:
            admitted.append(member)
    if not admitted:
        raise ValueError("No candidate FASTA ZIP members")
    return sorted(admitted, key=lambda member: member.filename)


def proteins(handle):
    seen = set()
    count = 0
    for protein in SeqIO.parse(handle, "fasta"):
        if not protein.id or protein.id in seen or not str(protein.seq).upper().removesuffix("*"):
            raise ValueError("Empty sequence or duplicate/empty FASTA identifier")
        seen.add(protein.id)
        count += 1
        yield protein.id, fingerprint(protein.seq)
    if not count:
        raise ValueError("Empty FASTA input")


def audit(archives, development):
    started, cpu_started = time.monotonic(), time.process_time()
    paths = [Path(path).resolve() for path in archives]
    if not paths or len(set(paths)) != len(paths) or not development:
        raise ValueError("Require unique candidate archives and development inputs")
    identities = [record(path) for path in paths]
    source_paths = [Path(__file__), *[Path(__file__).with_name(name) for name in (
        "audit_ygob_overlap.py", "prepare_ob_candidate_neighborhood.py",
        "orthobench_stage_diagnostics.py", "run_membership_hillclimb.py")]]
    sources = [record(path) for path in source_paths]
    index, candidates, members, seen_members = defaultdict(list), [], [], set()
    for path in paths:
        with ZipFile(path) as archive:
            for member in fasta_members(archive):
                key = str(path) + "::" + member.filename
                if member.filename in seen_members:
                    raise ValueError("Duplicate FASTA member name across archives")
                seen_members.add(member.filename)
                count = 0
                with archive.open(member) as binary, io.TextIOWrapper(binary, encoding="ascii") as handle:
                    for identifier, digest in proteins(handle):
                        index[digest].append(len(candidates))
                        candidates.append((key, identifier))
                        count += 1
                members.append({"archive": str(path), "member": member.filename, "uncompressed_bytes": member.file_size,
                                "proteins": count})
    candidate_id_counts = Counter(identifier for _, identifier in candidates)
    overlaps, development_inputs = {}, {}
    for label, files in development.items():
        files = [Path(path).resolve() for path in files]
        if not label or not files or len(set(files)) != len(files):
            raise ValueError("Require nonempty development labels and unique FASTAs")
        inputs = [record(path) for path in files]
        development_inputs[label] = inputs
        identities.extend(inputs)
        matched, total, matching_records = set(), 0, 0
        for path in files:
            with path.open(encoding="ascii") as handle:
                for _, digest in proteins(handle):
                    total += 1
                    matching = index.get(digest, ())
                    if matching:
                        matching_records += 1
                        matched.update(matching)
        overlaps[label] = {"development_proteins": total, "matching_development_records": matching_records,
            "matching_candidate_proteins": len(matched), "matching_candidate_fraction": len(matched) / len(candidates),
            "matching_candidate_records_by_member": dict(sorted(Counter(candidates[i][0] for i in matched).items()))}
    for identity in [*identities, *sources]:
        check(identity)
    return {"status": "candidate_sequence_exposure_audited", "scoring_performed": False,
        "inference_launched": False, "reference_classifications_read": False, "independence_established": False,
        "candidate_proteins": len(candidates), "candidate_unique_normalized_sequences": len(index),
        "candidate_fasta_members": members,
        "candidate_identifier_collisions_across_members": sum(count - 1 for count in candidate_id_counts.values()),
        "overlaps": overlaps, "strict_sequence_disjoint_confirmation_ruled_out": any(
            row["matching_candidate_proteins"] > 0 for row in overlaps.values()),
        "archives": identities[:len(paths)], "development_inputs": development_inputs, "sources": sources,
        "sequence_normalization": "Uppercase and remove one terminal stop only; preserve internal stops, ambiguous residues and gaps.",
        "resources": {"scope": "Single-process exposure audit, not tool inference timing; shared-host wall time may be confounded.",
                      "wall_s": time.monotonic() - started, "cpu_s": time.process_time() - cpu_started,
                      "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},
        "limitations": ["Sequence overlap disproves strict sequence disjointness; absent exact matches do not establish family or taxonomic independence.",
            "Matching records count every candidate occurrence, not only distinct sequences; identifier renaming cannot hide sequence reuse.",
            "No labels or predictions are examined; this does not validate the reference classification or score orthology accuracy."]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    identity = record(args.manifest)
    manifest = json.loads(args.manifest.read_text())
    expected = manifest["expected_development_proteins"]
    if (set(expected) != set(manifest["development_inputs"]) or not expected
            or any(isinstance(count, bool) or not isinstance(count, int) or count < 1 for count in expected.values())):
        raise ValueError("Development inventory/count schema differs")
    count = manifest["expected_candidate_proteomes"]
    if isinstance(count, bool) or not isinstance(count, int) or count < 1:
        raise ValueError("Invalid expected candidate proteome count")
    for item in [manifest["protocol"], *manifest["archives"]]:
        check(item)
    development = {label: [Path(item["path"]) for item in rows] for label, rows in manifest["development_inputs"].items()}
    for rows in manifest["development_inputs"].values():
        for item in rows:
            check(item)
    result = audit([Path(item["path"]) for item in manifest["archives"]], development)
    if len(result["candidate_fasta_members"]) != manifest["expected_candidate_proteomes"]:
        raise ValueError("Candidate proteome inventory differs")
    for label, expected in manifest["expected_development_proteins"].items():
        if result["overlaps"][label]["development_proteins"] != expected:
            raise ValueError("Development protein universe differs")
    check(identity)
    result.update(generated_at_utc=datetime.now(timezone.utc).isoformat(), manifest=identity,
                  protocol=manifest["protocol"], upstream=manifest["upstream"], goal_completed=False,
                  command=[sys.executable, *sys.argv])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_json(args.output, result)
    print(json.dumps({"candidate_proteins": result["candidate_proteins"],
        "overlaps": {label: row["matching_candidate_proteins"] for label, row in result["overlaps"].items()},
        "independence_established": False}))


if __name__ == "__main__":
    main()
