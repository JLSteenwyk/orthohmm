"""Collect the full input-only historical annotation panel, preserving missingness."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import time

from Bio import SeqIO
from benchmark_tools.audit_unisave_fragment_source import acquire
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_blast_recovery_batch import save_status
from benchmark_tools.run_simulation_methods import read_frozen

INVENTORY_SHA = "912363276fa5456f11aeff9f398fce71aec8d377c8c8eda7b144105945e554f1"
PROTOCOL_SHA = "447714557c357e752152ab4d2c95e624d7250195a7f0a59a1a1995c19f7f8359"
HELPER_SHA = "0642beeee11199853518b328b15bdb12d297d3a1324bcd3327815ccbc7355680"


def extract(families, descriptors, fastas):
    wanted = {g for members in families.values() for g in members}
    if wanted != set(descriptors) or any(not m or len(set(m)) != len(m) for m in families.values()):
        raise ValueError("Invalid family/descriptor coverage")
    genes = {}
    for item in fastas:
        check(item)
        for entry in SeqIO.parse(item["path"], "fasta"):
            parts = entry.id.split("|")
            if len(parts) != 3 or parts[1] not in wanted:
                continue
            accession = parts[1]
            descriptor = descriptors[accession]
            if (accession in genes or entry.description != descriptor["description"]
                    or len(entry.seq) != descriptor["length"]
                    or entry.id != descriptor["input_id"]
                    or Path(item["path"]).name != descriptor["source_file"]):
                raise ValueError("Duplicate or changed corrected sequence descriptor")
            fields = {}
            for key in ("OX", "SV"):
                values = re.findall(rf"(?:^|\s){key}=(\d+)(?=\s|$)", entry.description)
                if len(values) != 1 or int(values[0]) <= 0:
                    raise ValueError("Missing or ambiguous taxon/sequence version")
                fields[key] = values[0]
            genes[accession] = dict(accession=accession, sequence_version=int(fields["SV"]),
                sequence_sha256=hashlib.sha256(str(entry.seq).encode("ascii")).hexdigest(),
                taxid=fields["OX"], length=len(entry.seq))
        check(item)
    if set(genes) != wanted:
        raise ValueError("Incomplete corrected sequence extraction")
    return genes


def prepare(root):
    results = root / "benchmark_tools/results"
    inventory_path = results / "corrected_swiss_sequence_strata_20260918.json"
    inventory = read_frozen(inventory_path, INVENTORY_SHA)
    protocol = record(results / "SWISS_HISTORICAL_FRAGMENT_PROTOCOL_20260923.md")
    helper = record(Path(__file__).with_name("audit_unisave_fragment_source.py"))
    if protocol["sha256"] != PROTOCOL_SHA or helper["sha256"] != HELPER_SHA:
        raise ValueError("Changed frozen acquisition protocol/helper")
    families = inventory["family_memberships"]
    if len(families) != 18 or len(inventory["genes"]) != 563 or len(inventory["fasta_inputs"]) != 78:
        raise ValueError("Changed panel size")
    genes = extract(families, inventory["genes"], inventory["fasta_inputs"])
    records = [record(inventory_path), protocol, helper, record(__file__), *inventory["fasta_inputs"]]
    for item in records:
        check(item)
    return dict(families=families, genes=genes, records=records)


def collect(preflight, output, delay=1.0):
    if delay < 1.0:
        raise ValueError("Require at least one second between accession requests")
    output.mkdir(parents=True, exist_ok=False)
    report = dict(status="collecting", preflight=preflight, entries={},
        job_id=os.environ.get("SLURM_JOB_ID"), prediction_statistics_evaluated=False,
        annotation_panel_admitted=False, publication_ready=False)
    save_status(output / "status.json", report)
    for accession, gene in sorted(preflight["genes"].items()):
        time.sleep(delay)
        try:
            result = acquire(accession, gene["sequence_version"], gene["sequence_sha256"],
                             gene["taxid"], output / accession)
            report["entries"][accession] = dict(status="sequence_matched", audit=record(output / accession / "audit.json"),
                selection_class=result["selection_class"])
        except Exception as error:
            report["entries"][accession] = dict(status="missing", error_type=type(error).__name__, error=str(error))
        save_status(output / "status.json", report)
    for item in preflight["records"]:
        check(item)
    report["status"] = "collection_complete_pending_independent_validation"
    report["matched"] = sum(r["status"] == "sequence_matched" for r in report["entries"].values())
    report["missing"] = len(report["entries"]) - report["matched"]
    save_status(output / "status.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    prepared = prepare(args.root.resolve())
    if args.check_only:
        print(json.dumps(dict(genes=len(prepared["genes"]), families=len(prepared["families"]))))
    elif args.output is None:
        parser.error("--output is required unless --check-only")
    else:
        collect(prepared, args.output.absolute())
