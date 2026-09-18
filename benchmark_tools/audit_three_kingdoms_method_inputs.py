"""Compare retained method inputs without equating commands with consumed bytes."""

import argparse
import csv
import json
from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_simulation_methods import read_frozen

SOURCE_SHA = "a79bf83e1ea28a9790e4597ea50a44d6409f4f36404a35e660d56a0b754fc1f3"
PARITY = "three_kingdoms/results/parity_20260907/"
METHODS = {
    "orthohmm_high_sensitivity": ("three_kingdoms/results/orthohmm_high_sensitivity", None, None),
    "orthohmm_phylogeny_satellite_v2": (PARITY + "orthohmm_phylogeny_satellite_v2", None, None),
    "orthofinder_3_1_5_full": (PARITY + "orthofinder_3_1_5_full", "run/input", ".fasta"),
    "orthofinder_3_1_5_sequence_only": (PARITY + "orthofinder_3_1_5_sequence_only", None, None),
    "proteinortho_6_3_6": ("three_kingdoms/results/proteinortho", "run/input", ".fasta"),
    "sonicparanoid_2_0_9": ("three_kingdoms/results/sonicparanoid", "run/input", ".fasta"),
    "fastoma_0_3_5": (PARITY + "fastoma_0_3_5", "run/input/proteome", ".fa"),
    "orthomcl_1_4": (PARITY + "orthomcl_1_4", None, None),
}


def parse_sha256(text):
    rows = {}
    for line in text.splitlines():
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)", line)
        if not match or match[2] in rows:
            raise ValueError("Malformed or duplicate input checksum record")
        rows[match[2]] = match[1]
    if not rows:
        raise ValueError("Empty input checksum inventory")
    return rows


def classify(digest, canonical, raw):
    if digest == canonical:
        return "matches_staged"
    if digest == raw:
        return "matches_raw_not_staged"
    return "matches_neither"


def snapshot_rows(text):
    rows = {}
    for index, row in enumerate(csv.reader(text.splitlines(), delimiter="\t"), 1):
        if (len(row) != 5 or row[0] != str(index) or row[1] in rows
                or not re.fullmatch(r"[0-9a-f]{64}", row[2])
                or not row[3].isdigit() or not row[4].isdigit()):
            raise ValueError("Malformed SonicParanoid snapshot")
        rows[row[1]] = {"sha256": row[2], "proteins": int(row[3]), "reported_residues": int(row[4])}
    return rows


def audit(repo):
    source_path = repo / "benchmark_tools/results/three_kingdoms_sources_20260918.json"
    source = read_frozen(source_path, SOURCE_SHA)
    if source["status"] != "retained_three_kingdoms_lineage_and_reference_verified" or len(source["inputs"]) != 12:
        raise ValueError("Unexpected input source audit")
    expected = {r["code"]: r for r in source["inputs"]}
    records = [record(source_path)]
    for row in expected.values():
        for kind in ("raw", "staged"):
            check(row["files"][kind])
            records.append(row["files"][kind])
    methods = []
    for name, (directory, copy_dir, suffix) in METHODS.items():
        root = repo / directory
        evidence, copies = [], []
        sha_file = root / "input.sha256"
        manifest_status = "not_retained"
        if sha_file.exists():
            identity = record(sha_file)
            records.append(identity)
            evidence.append(identity)
            parsed = parse_sha256(sha_file.read_text())
            canonical = {r["files"]["staged"]["path"]: r["files"]["staged"]["sha256"] for r in expected.values()}
            if parsed != canonical:
                raise ValueError("Retained method input checksum manifest differs: " + name)
            manifest_status = "all_12_staged_hashes_match"
        if copy_dir:
            directory = root / copy_dir
            if {p.stem for p in directory.glob("*" + suffix)} != set(expected):
                raise ValueError("Unexpected copied input inventory")
            for code, row in sorted(expected.items()):
                path = directory / (code + suffix)
                if path.is_symlink():
                    raise ValueError("A symlink is not a retained independent input copy")
                item = record(path)
                records.append(item)
                copies.append({"code": code, "file": item,
                               "classification": classify(item["sha256"], row["files"]["staged"]["sha256"], row["files"]["raw"]["sha256"])})
        native_snapshot = None
        if name == "sonicparanoid_2_0_9":
            path = root / "run/sp_out/snapshot.tsv"
            item = record(path)
            records.append(item)
            evidence.append(item)
            native_snapshot = snapshot_rows(path.read_text())
            if set(native_snapshot) != {code + ".fasta" for code in expected}:
                raise ValueError("SonicParanoid snapshot species differ")
            for row in copies:
                native = native_snapshot[row["code"] + ".fasta"]
                if native["sha256"] != row["file"]["sha256"] or native["proteins"] != expected[row["code"]]["staged_proteins"]:
                    raise ValueError("SonicParanoid snapshot differs from retained input")
        if name == "orthohmm_high_sensitivity":
            path = root / "metrics.json"
            item = record(path)
            records.append(item)
            evidence.append(item)
            metrics = json.loads(path.read_text())
            if metrics["metadata"]["fasta_directory"] != str(repo / "three_kingdoms/input"):
                raise ValueError("Unexpected high-sensitivity input path")
        methods.append({"method": name, "input_hash_manifest": manifest_status, "copies": copies,
                        "native_snapshot": native_snapshot, "evidence": evidence,
                        "historical_consumption_of_staged_bytes_proven": False})
    for item in records:
        check(item)
    return {"status": "retained_method_input_evidence_reviewed", "source": record(__file__),
            "input_source_audit": record(source_path), "methods": methods,
            "publication_ready": False, "uniform_historical_input_equivalence_established": False,
            "limitations": ["Saved hash manifests and current copies are evidence, not immutable execution attestation.",
                "The high-sensitivity record identifies an input directory but does not retain per-file input hashes.",
                "Sequence-only OrthoFinder is a checkpoint diagnostic from the full run, not an independent native run.",
                "SonicParanoid snapshot verifies input file hashes; downstream effective sequence handling is not traced here.",
                "Input differences do not establish their impact on final predictions or reference scores.",
                "No old result is overwritten, rerun, rescored or admitted by this audit."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.repo.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(result["status"])
